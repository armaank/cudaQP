/* interface for pcg_solver with main osqp api */

#include "cuda_pcg_interface.h"
#include "cuda_pcg.h"

#include "cuda_lin_alg.h"
#include "cuda_malloc.h"

#include "glob_opts.h"

static c_float compute_tolerance(cudapcg_solver *s, c_int admm_iter)
{
    c_float eps;
    c_float rhs_norm;

    /* compute the norm of RHS of the linear system */
    s->vector_norm(s->d_rhs, s->n, &rhs_norm);

    if (s->polish) return c_max(rhs_norm * CUDA_PCG_POLISH_ACCURACY, CUDA_PCG_EPS_MIN);

    switch (s->eps_strategy)
    {

    /* SCS strategy */
    case SCS_STRATEGY:
        eps = rhs_norm * s->start_tol  / pow((admm_iter + 1), s->dec_rate);
        eps = c_max(eps, CUDA_PCG_EPS_MIN);
        break;

    /* residual strategy */
    case RESIDUAL_STRATEGY:
        if (admm_iter == 1)
        {
            /* In case rhs = 0.0 we don't want to set eps_prev to 0.0 */
            if (rhs_norm < CUDA_PCG_EPS_MIN)
                s->eps_prev = 1.0;
            else
                s->eps_prev = rhs_norm * s->reduction_factor;
            /* Return early since scaled_pri_res and scaled_dua_res are meaningless before the first ADMM iteration */
            return s->eps_prev;
        }

        if (s->zero_pcg_iter >= s->reduction_threshold) {
            s->reduction_factor /= 2;
            s->zero_pcg_iter = 0;
        }

        eps = s->reduction_factor * sqrt((*s->scaled_pri_res) * (*s->scaled_dua_res));
        eps = c_max(c_min(eps, s->eps_prev), CUDA_PCG_EPS_MIN);
        s->eps_prev = eps;
        break;
    }
    return eps;
}

/* d_rhs = d_b1 + A' * rho * d_b2 */
static void compute_rhs(cudapcg_solver *s, c_float *d_b)
{

    c_int n = s->n;
    c_int m = s->m;

    /* d_rhs = d_b1 */
    cuda_vec_copy_d2d(s->d_rhs, d_b, n);

    if (m == 0)
        return;

    /* d_z = d_b2 */
    cuda_vec_copy_d2d(s->d_z, d_b + n, m);

    if (!s->d_rho_vec)
    {
        /* d_z *= rho */
        cuda_vec_mult_sc(s->d_z, *s->h_rho, m);
    }
    else
    {
        /* d_z = diag(d_rho_vec) * d_z */
        cuda_vec_ew_prod(s->d_z, s->d_z, s->d_rho_vec, m);
    }

    /* d_rhs += A' * d_z */
    cuda_mat_Axpy(s->At, s->d_z, s->d_rhs, 1.0, 1.0);
}


/* api fcns as perscribed by the osqp documentation */
c_int init_linsys_solver_cudapcg(cudapcg_solver **sp,
                                 const OSQPMatrix *P,
                                 const OSQPMatrix *A,
                                 const OSQPVectorf *rho_vec,
                                 OSQPSettings *settings,
                                 c_float *scaled_pri_res,
                                 c_float *scaled_dua_res,
                                 c_int polish) {

    c_int n, m;
    c_float H_m_1 = -1.0;

    /* allocate linsys solver structure */
    cudapcg_solver *s = c_calloc(1, sizeof(cudapcg_solver));
    *sp = s;

    /* assign type and the number of threads */
    s->type = CUDA_PCG_SOLVER;
    s->nthreads = 1; // dummy, this changes for target device

    /* dimensions */
    n = OSQPMatrix_get_n(P);
    m = OSQPMatrix_get_m(A);
    s->n = n;
    s->m = m;

    /* pcg states */
    s->polish = polish;
    s->zero_pcg_iter = 0;

    /* default norm and tolerance strategy */
    s->eps_strategy = RESIDUAL_STRATEGY;
    s->norm = CUDA_PCG_NORM;
    s->precondition = CUDA_PCG_PRECONDITION;
    s->warm_start_pcg = CUDA_PCG_WARM_START;
    s->max_iter = (polish) ? CUDA_PCG_POLISH_MAX_ITER : CUDA_PCG_MAX_ITER;

    /* tolerance strategy parameters */
    s->start_tol = CUDA_PCG_START_TOL;
    s->dec_rate = CUDA_PCG_DECAY_RATE;
    s->reduction_threshold = CUDA_PCG_REDUCTION_THRESHOLD;
    s->reduction_factor = CUDA_PCG_REDUCTION_FACTOR;
    s->scaled_pri_res = scaled_pri_res;
    s->scaled_dua_res = scaled_dua_res;

    /* set pointers settings and data */
    s->A = A->S;
    s->At = A->At;
    s->P = P->S;
    s->d_P_diag_ind = P->d_P_diag_ind;
    if (rho_vec)
        s->d_rho_vec = rho_vec->d_val;
    if (!polish)
    {
        s->h_sigma = &settings->sigma;
        s->h_rho = &settings->rho;
    }
    else
    {
        s->h_sigma = &settings->delta;
        s->h_rho = (c_float*) c_malloc(sizeof(c_float));
        *s->h_rho = 1. / settings->delta;
    }

    /* allocate pcg iterates */
    cuda_calloc((void **) &s->d_x, n * sizeof(c_float));
    cuda_malloc((void **) &s->d_p, n * sizeof(c_float));
    cuda_malloc((void **) &s->d_Kp, n * sizeof(c_float));
    cuda_malloc((void **) &s->d_y, n * sizeof(c_float));
    cuda_malloc((void **) &s->d_r, n * sizeof(c_float));
    cuda_malloc((void **) &s->d_rhs, n * sizeof(c_float));

    if (m != 0) cuda_malloc((void **) &s->d_z, m * sizeof(c_float));

    /* allocate scalar in host memory that is page-locked and accessible to target */
    cuda_malloc_host((void **) &s->h_r_norm, sizeof(c_float));

    /* allocate target-side scalar values. This way scalars are packed in target memory */
    cuda_malloc((void **) &s->d_r_norm, 8 * sizeof(c_float));
    s->rTy = s->d_r_norm + 1;
    s->rTy_prev = s->d_r_norm + 2;
    s->alpha = s->d_r_norm + 3;
    s->beta = s->d_r_norm + 4;
    s->pKp = s->d_r_norm + 5;
    s->D_MINUS_ONE = s->d_r_norm + 6;
    s->d_sigma = s->d_r_norm + 7;
    cuda_vec_copy_h2d(s->D_MINUS_ONE, &H_m_1, 1);
    cuda_vec_copy_h2d(s->d_sigma, s->h_sigma, 1);

    /* allocate memory for PCG preconditioning */
    if (s->precondition)
    {
        cuda_malloc((void **) &s->d_P_diag_val, n * sizeof(c_float));
        cuda_malloc((void **) &s->d_AtRA_diag_val, n * sizeof(c_float));
        cuda_malloc((void **) &s->d_diag_precond, n * sizeof(c_float));
        cuda_malloc((void **) &s->d_diag_precond_inv, n * sizeof(c_float));
        if (!s->d_rho_vec) cuda_malloc((void **) &s->d_AtA_diag_val, n * sizeof(c_float));
    }

    /* Set the vector norm */
    switch (s->norm) {
    case 0:
        s->vector_norm = &cuda_vec_norm_inf;
        break;

    case 2:
        s->vector_norm = &cuda_vec_norm_2;
        break;
    }

    s->solve = &solve_linsys_cudapcg;
    s->warm_start = &warm_start_linsys_solver_cudapcg;
    s->free = &free_linsys_solver_cudapcg;
    s->update_matrices = &update_linsys_solver_matrices_cudapcg;
    s->update_rho_vec  = &update_linsys_solver_rho_vec_cudapcg;

    /* init PCG preconditioner */
    if (s->precondition) cuda_pcg_update_precond(s, 1, 1, 1);

    return 0;
}

/* main driver function for pcg */
c_int solve_linsys_cudapcg(cudapcg_solver *s, OSQPVectorf *b,c_int admm_iter)
{

    c_int pcg_iter;
    c_float eps;

    /* compute the RHS of the reduced KKT system */
    compute_rhs(s, b->d_val);

    /* compute the required solution precision */
    eps = compute_tolerance(s, admm_iter);

    /* solve the linear system with PCG */
    pcg_iter = cuda_pcg(s, eps, s->max_iter);

    /* copy the first part of the solution to b->d_val */
    cuda_vec_copy_d2d(b->d_val, s->d_x, s->n);

    /* solution polishing */
    if (!s->polish)
    {
        /* Compute d_z = A * d_x */
        if (s->m) cuda_mat_Axpy(s->A, s->d_x, b->d_val + s->n, 1.0, 0.0);
    }
    else
    {
        /* Compute yred = (A * d_x - b) / delta */
        cuda_mat_Axpy(s->A, s->d_x, b->d_val + s->n, 1.0, -1.0);
        cuda_vec_mult_sc(b->d_val + s->n, *s->h_rho, s->m);
    }

    if (pcg_iter == 0)
        s->zero_pcg_iter++;

    return 0;
}

/* remove me */
void warm_start_linsys_solver_cudapcg(cudapcg_solver *s, const OSQPVectorf *x)
{
    cuda_vec_copy_d2d(s->d_x, x->d_val, x->length);
}


void free_linsys_solver_cudapcg(cudapcg_solver *s)
{

    if (s)
    {
        if (s->polish)
            c_free(s->h_rho);

        /* PCG iterates */
        cuda_free((void **) &s->d_x);
        cuda_free((void **) &s->d_p);
        cuda_free((void **) &s->d_Kp);
        cuda_free((void **) &s->d_y);
        cuda_free((void **) &s->d_r);
        cuda_free((void **) &s->d_rhs);
        cuda_free((void **) &s->d_z);

        /* free host memory */
        cuda_free_host((void **) &s->h_r_norm);

        /* target-side scalar values */
        cuda_free((void **) &s->d_r_norm);

        /* pcg preconditioner */
        cuda_free((void **) &s->d_P_diag_val);
        cuda_free((void **) &s->d_AtA_diag_val);
        cuda_free((void **) &s->d_AtRA_diag_val);
        cuda_free((void **) &s->d_diag_precond);
        cuda_free((void **) &s->d_diag_precond_inv);

        c_free(s);
    }
}


c_int update_linsys_solver_matrices_cudapcg(cudapcg_solver *s, const OSQPMatrix *P, const OSQPMatrix *A)
{

    if (s->precondition) cuda_pcg_update_precond(s, 1, 1, 0);
    return 0;
}


c_int update_linsys_solver_rho_vec_cudapcg(cudapcg_solver *s, const OSQPVectorf *rho_vec, c_float rho_sc)
{

    if (s->precondition) cuda_pcg_update_precond(s, 0, 0, 1);
    return 0;
}

