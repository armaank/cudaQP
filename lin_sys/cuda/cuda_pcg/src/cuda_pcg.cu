/* cuda pcg algorithm */

#include "cuda_pcg.h"
#include "csr_type.h"
#include "cuda_handler.h"
#include "cuda_malloc.h"
#include "cuda_lin_alg.h"
#include "cuda_wrapper.h"
#include "helper_cuda.h"

#ifdef __cplusplus
extern "C" {
extern CUDA_Handle_t *CUDA_handle;
}
#endif

__global__ void scalar_division_kernel(c_float *res, const c_float *num, const c_float *den)
{
    *res = (*num) / (*den);
}

/* computes:  d_y = (P + sigma*I + A'*R*A) * d_x */
static void mat_vec_prod(cudapcg_solver *s, c_float *d_y, const c_float  *d_x, c_int device)
{

    c_float *sigma;
    c_float H_0 = 0.0;
    c_float H_1  = 1.0;
    c_int n = s->n;
    c_int m = s->m;
    csr *P  = s->P;
    csr *A  = s->A;
    csr *At = s->At;

    if (device)
    {
        sigma = s->d_sigma;
    }
    else
    {
        sigma = s->h_sigma;
    }

    /* d_y = d_x */
    checkCudaErrors(cudaMemcpy(d_y, d_x, n * sizeof(c_float), cudaMemcpyDeviceToDevice));

    /* d_y *= sigma */
    checkCudaErrors(cublasTscal(CUDA_handle->cublasHandle, n, sigma, d_y, 1));

    /* d_y += P * d_x */
    checkCudaErrors(cusparseCsrmv(CUDA_handle->cusparseHandle, P->alg, P->m, P->n, P->nnz, &H_1, P->MatDescription, P->val, P->row_ptr, P->col_ind, d_x, &H_1, d_y, P->buffer));

    if (m == 0) return;

    if (!s->d_rho_vec)
    {
        /* d_z = rho * A * d_x */
        checkCudaErrors(cusparseCsrmv(CUDA_handle->cusparseHandle, A->alg, A->m, A->n, A->nnz, s->h_rho, A->MatDescription, A->val, A->row_ptr, A->col_ind, d_x, &H_0, s->d_z, A->buffer));
    }
    else
    {
        /* d_z = A * d_x */
        checkCudaErrors(cusparseCsrmv(CUDA_handle->cusparseHandle, A->alg, A->m, A->n, A->nnz, &H_1, A->MatDescription, A->val, A->row_ptr, A->col_ind, d_x, &H_0, s->d_z, A->buffer));

        /* d_z = diag(d_rho_vec) * dz */
        cuda_vec_ew_prod(s->d_z, s->d_z, s->d_rho_vec, m);
    }

    /* d_y += A' * d_z */
    checkCudaErrors(cusparseCsrmv(CUDA_handle->cusparseHandle, At->alg, At->m, At->n, At->nnz, &H_1, At->MatDescription, At->val, At->row_ptr, At->col_ind, s->d_z, &H_1, d_y, A->buffer));
}

/* pcg algorithm */
c_int cuda_pcg(cudapcg_solver *s, c_float eps, c_int max_niter)
{

    c_float *ptr_tmp;
    c_int niter = 0;
    c_int n = s->n;
    c_float H_m_1 = -1.0;

    /* set up problem */

    if (!s->warm_start)
    {
        /* d_x = 0 */
        checkCudaErrors(cudaMemset(s->d_x, 0, n * sizeof(c_float)));
    }

    /* d_p = 0 */
    checkCudaErrors(cudaMemset(s->d_p, 0, n * sizeof(c_float)));

    /* d_r = K * d_x */
    mat_vec_prod(s, s->d_r, s->d_x, 0);

    /* d_r -= d_rhs */
    checkCudaErrors(cublasTaxpy(CUDA_handle->cublasHandle, n, &H_m_1, s->d_rhs, 1, s->d_r, 1));

    /* h_r_norm = |d_r| */
    s->vector_norm(s->d_r, n, s->h_r_norm);

    /* need to change CUBLAS mode */
    cublasSetPointerMode(CUDA_handle->cublasHandle, CUBLAS_POINTER_MODE_DEVICE);

    if (s->precondition)
    {
        /* d_y = M \ d_r */
        cuda_vec_ew_prod(s->d_y, s->d_diag_precond_inv, s->d_r, n);
    }

    /* d_p = -d_y */
    checkCudaErrors(cublasTaxpy(CUDA_handle->cublasHandle, n, s->D_MINUS_ONE, s->d_y, 1, s->d_p, 1));

    /* rTy = d_r' * d_y */
    checkCudaErrors(cublasTdot(CUDA_handle->cublasHandle, n, s->d_y, 1, s->d_r, 1, s->rTy));

    /* synchronize for timing */
    cudaDeviceSynchronize();

    /* Run the PCG algorithm */
    while ( *(s->h_r_norm) > eps && niter < max_niter )
    {

        /* d_Kp = K * d_p */
        mat_vec_prod(s, s->d_Kp, s->d_p, 1);

        /* pKp = d_p' * d_Kp */
        checkCudaErrors(cublasTdot(CUDA_handle->cublasHandle, n, s->d_p, 1, s->d_Kp, 1, s->pKp));

        /* alpha = rTy / pKp */
        scalar_division_kernel<<<1,1>>>(s->alpha, s->rTy, s->pKp);

        /* d_x += alpha * d_p */
        checkCudaErrors(cublasTaxpy(CUDA_handle->cublasHandle, n, s->alpha, s->d_p, 1, s->d_x, 1));

        /* d_r += alpha * d_Kp */
        checkCudaErrors(cublasTaxpy(CUDA_handle->cublasHandle, n, s->alpha, s->d_Kp, 1, s->d_r, 1));

        if (s->precondition)
        {
            /* d_y = M \ d_r */
            cuda_vec_ew_prod(s->d_y, s->d_diag_precond_inv, s->d_r, n);
        }

        /* Swap pointers to rTy and rTy_prev */
        ptr_tmp = s->rTy_prev;
        s->rTy_prev = s->rTy;
        s->rTy = ptr_tmp;

        /* rTy = d_r' * d_y */
        checkCudaErrors(cublasTdot(CUDA_handle->cublasHandle, n, s->d_y, 1, s->d_r, 1, s->rTy));

        /* Update residual norm */
        s->vector_norm(s->d_r, n, s->d_r_norm);
        checkCudaErrors(cudaMemcpyAsync(s->h_r_norm, s->d_r_norm, sizeof(c_float), cudaMemcpyDeviceToHost));

        /* beta = rTy / rTy_prev */
        scalar_division_kernel<<<1,1>>>(s->beta, s->rTy, s->rTy_prev);

        /* d_p *= beta */
        checkCudaErrors(cublasTscal(CUDA_handle->cublasHandle, n, s->beta, s->d_p, 1));

        /* d_p -= d_y */
        checkCudaErrors(cublasTaxpy(CUDA_handle->cublasHandle, n, s->D_MINUS_ONE, s->d_y, 1, s->d_p, 1));

        cudaDeviceSynchronize();
        niter++;

    } /* End of the PCG algorithm */

    /* change CUBLAS pointer mode back */
    cublasSetPointerMode(CUDA_handle->cublasHandle, CUBLAS_POINTER_MODE_HOST);

    return niter;
}

/* update preconditioning  */
void cuda_pcg_update_precond(cudapcg_solver *s, c_int P_updated, c_int A_updated, c_int R_updated)
{

    void    *buffer;
    c_float *mem_tmp;
    c_int    n  = s->n;
    csr     *At = s->At;

    size_t buff_size = n * (sizeof(c_float) + sizeof(c_int));

    if (!P_updated && !A_updated && !R_updated) return;

    if (P_updated)
    {
        /* Update d_P_diag_val */
        checkCudaErrors(cusparseTgthr(CUDA_handle->cusparseHandle, n, s->P->val, s->d_P_diag_val, s->d_P_diag_ind, CUSPARSE_INDEX_BASE_ZERO));
    }

    if (A_updated || R_updated)
    {
        /* Allocate memory */
        cuda_malloc((void **) &mem_tmp, At->nnz * sizeof(c_float));
        cuda_malloc((void **) &buffer, buff_size);

        /* Update d_AtRA_diag_val */
        if (!s->d_rho_vec)
        {   /* R = rho*I  -->  A'*R*A = rho * A'*A */

            if (A_updated)
            {
                /* Update d_AtA_diag_val */
                cuda_vec_ew_prod(mem_tmp, At->val, At->val, At->nnz);
                cuda_vec_segmented_sum(mem_tmp, At->row_ind, s->d_AtA_diag_val, buffer, n, At->nnz);
            }

            /* d_AtRA_diag_val = rho * d_AtA_diag_val */
            cuda_vec_add_scaled(s->d_AtRA_diag_val, s->d_AtA_diag_val, NULL, *s->h_rho, 0.0, n);
        }
        else
        {   /* R = diag(d_rho_vec)  -->  A'*R*A = A' * diag(d_rho_vec) * A */
            cuda_mat_rmult_diag_new(At, mem_tmp, s->d_rho_vec);   /* mem_tmp = A' * R */
            cuda_vec_ew_prod(mem_tmp, mem_tmp, At->val, At->nnz);     /* mem_tmp = mem_tmp * A */
            cuda_vec_segmented_sum(mem_tmp, At->row_ind, s->d_AtRA_diag_val, buffer, n, At->nnz);
        }

        cuda_free((void **) &mem_tmp);
        cuda_free((void **) &buffer);
    }

    /* d_diag_precond = sigma */
    cuda_vec_set_sc(s->d_diag_precond, *s->h_sigma, n);

    /* d_diag_precond += d_P_diag_val + d_AtRA_diag_val */
    cuda_vec_add_scaled3(s->d_diag_precond, s->d_diag_precond, s->d_P_diag_val, s->d_AtRA_diag_val, 1.0, 1.0, 1.0, n);

    /* d_diag_precond_inv = 1 / d_diag_precond */
    cuda_vec_reciprocal(s->d_diag_precond_inv, s->d_diag_precond, n);
}
