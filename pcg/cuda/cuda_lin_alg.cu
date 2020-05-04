/**
 *  Copyright (c) 2019 ETH Zurich, Automatic Control Lab, Michel Schubiger, Goran Banjac.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

 #include "./include/cuda_lin_alg.h"
//  #include "cuda_configure.h"
//  #include "cuda_handler.h"
//  #include "cuda_malloc.h"
//  #include "cuda_wrapper.h"
//  #include "helper_cuda.h"    /* --> checkCudaErrors */
 
//  #include "csr_type.h"
// //  #include "glob_opts.h"

 
//  #include <thrust/reduce.h>
//  #include <thrust/execution_policy.h>

 
 
 /*******************************************************************************
  *                              GPU Kernels                                    *
  *******************************************************************************/
 
  __global__ void vec_set_sc_kernel(float *a,
                                    float  sc,
                                    int    n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < n; i += grid_size) {
     a[i] = sc;
   }
 }
 
 __global__ void vec_set_sc_cond_kernel(float     *a,
                                        const int *test,
                                        float      sc_if_neg,
                                        float      sc_if_zero,
                                        float      sc_if_pos,
                                        int        n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < n; i += grid_size) {
     if (test[i] == 0)      a[i] = sc_if_zero;
     else if (test[i] > 0)  a[i] = sc_if_pos;
     else                   a[i] = sc_if_neg;
   }
 }
 
 __global__ void vec_prod_pos_kernel(const float *a,
                                     const float *b,
                                     float       *res,
                                     int          n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   float res_kernel = 0.0;
 
   for(int i = idx; i < n; i += grid_size) {
     res_kernel += a[i] * c_max(b[i], 0.0);
   }
   atomicAdd(res, res_kernel);
 }
 
 __global__ void vec_prod_neg_kernel(const float *a,
                                     const float *b,
                                     float       *res,
                                     int          n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   float res_kernel = 0.0;
 
   for(int i = idx; i < n; i += grid_size) {
     res_kernel += a[i] * c_min(b[i], 0.0);
   }
   atomicAdd(res, res_kernel);
 }
 
 __global__ void vec_ew_prod_kernel(float       *c,
                                    const float *a,
                                    const float *b,
                                    int          n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < n; i += grid_size) {
//  #ifdef DFLOAT
    //  c[i] = __fmul_rn(a[i], b[i]);
//  #else
     c[i] = __dmul_rn(a[i], b[i]);
//  #endif
   }
 }
 
 __global__ void vec_leq_kernel(const float *l,
                                const float *u,
                                int         *res,
                                int          n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < n; i += grid_size) {
     if (l[i] > u[i]) atomicAnd(res, 0);
   }
 }
 
 __global__ void vec_bound_kernel(float       *x,
                                  const float *z,
                                  const float *l,
                                  const float *u,
                                  int          n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < n; i += grid_size) {
     x[i] = c_min(c_max(z[i], l[i]), u[i]);
   }
 }
 
 __global__ void vec_project_polar_reccone_kernel(float       *y,
                                                  const float *l,
                                                  const float *u,
                                                  float        infval,
                                                  int          n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < n; i += grid_size) {
     if (u[i] > +infval) {
       if (l[i] < -infval) {
         /* Both bounds infinite */
         y[i] = 0.0;
       }
       else {
         /* Only upper bound infinite */
         y[i] = c_min(y[i], 0.0);
       }
     }
     else if (l[i] < -infval) {
       /* Only lower bound infinite */
       y[i] = c_max(y[i], 0.0);
     }
   }
 }
 
 __global__ void vec_in_reccone_kernel(const float *y,
                                       const float *l,
                                       const float *u,
                                       float        infval,
                                       float        tol,
                                       int         *res,
                                       int          n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < n; i += grid_size) {
     if ( (u[i] < +infval && y[i] > +tol) ||
          (l[i] > -infval && y[i] < -tol) )
       atomicAnd(res, 0);
   }
 }
 
 __global__ void vec_reciprocal_kernel(float       *b,
                                       const float *a,
                                       int          n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < n; i += grid_size) {
//  #ifdef DFLOAT
    //  b[i] = __frcp_rn(a[i]);
//  #else
     b[i] = __drcp_rn(a[i]);
//  #endif
   }
 }
 
 __global__ void vec_sqrt_kernel(float *a,
                                 int    n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < n; i += grid_size) {
//  #ifdef DFLOAT
     a[i] = __fsqrt_rn(a[i]);
//  #else
     a[i] = __dsqrt_rn(a[i]);
//  #endif
   }
 }
 
 __global__ void vec_max_kernel(float       *c,
                                const float *a,
                                const float *b,
                                int          n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < n; i += grid_size) {
     c[i] = c_max(a[i], b[i]);
   }
 }
 
 __global__ void vec_min_kernel(float       *c,
                                const float *a,
                                const float *b,
                                int          n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < n; i += grid_size) {
     c[i] = c_min(a[i], b[i]);
   }
 }
 
 __global__ void vec_bounds_type_kernel(int         *iseq,
                                        const float *l,
                                        const float *u,
                                        float        infval,
                                        float        tol,
                                        int         *has_changed,
                                        int          n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < n; i += grid_size) {
     if (u[i] - l[i] < tol) {
       /* Equality constraints */
       if (iseq[i] != 1) {
         iseq[i] = 1;
         atomicOr(has_changed, 1);
       }
     }
     else if ( (l[i] < -infval) && (u[i] > infval) ) {
       /* Loose bounds */
       if (iseq[i] != -1) {
         iseq[i] = -1;
         atomicOr(has_changed, 1);
       }
     }
     else {
       /* Inequality constraints */
       if (iseq[i] != 0) {
         iseq[i] = 0;
         atomicOr(has_changed, 1);
       }
     }
   }
 }
 
 __global__ void vec_set_sc_if_lt_kernel(float       *x,
                                         const float *z,
                                         float        testval,
                                         float        newval,
                                         int          n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < n; i += grid_size) {
     x[i] = z[i] < testval ? newval : z[i];
   }
 }
 
 __global__ void vec_set_sc_if_gt_kernel(float       *x,
                                         const float *z,
                                         float        testval,
                                         float        newval,
                                         int          n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < n; i += grid_size) {
     x[i] = z[i] > testval ? newval : z[i];
   }
 }
 
 __global__ void mat_lmult_diag_kernel(const int   *row_ind,
                                       const float *diag,
                                       float       *data,
                                       int          nnz) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < nnz; i += grid_size) {
     int row = row_ind[i];
     data[i] *= diag[row];
   }
 }
 
 __global__ void mat_rmult_diag_kernel(const int   *col_ind,
                                       const float *diag,
                                       float       *data,
                                       int          nnz) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < nnz; i += grid_size) {
     int column = col_ind[i];
     data[i] *= diag[column];
   }
 }
 
 __global__ void mat_rmult_diag_new_kernel(const int   *col_ind,
                                           const float *diag,
                                           const float *data_in,
                                           float       *data_out,
                                           int          nnz) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < nnz; i += grid_size) {
     int column = col_ind[i];
     data_out[i] = data_in[i] * diag[column];
   }
 }
 
 __global__ void vec_abs_kernel(float *a,
                                int    n) {
 
   int i  = threadIdx.x + blockDim.x * blockIdx.x;
 
   if (i < n) {
//  #ifdef DFLOAT
     a[i] = fabsf(a[i]);
//  #else
     a[i] = fabs(a[i]);
//  #endif
   }
 }
 
 __global__ void scatter_kernel(float       *out,
                                const float *in,
                                const int   *ind,
                                int          n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < n; i += grid_size) {
     int j = ind[i];
     out[j] = in[i];
   }
 }
 
 /*
  * This code complements the cublasITamax routine which only returns the 
  * one-based index to the maximum absolute value in d_x. 
 */
 __global__ void abs_kernel(const int   *index_one_based,
                            const float *d_x,
                            float       *res) {
 
   /* cublasITamax returns one-based index */
   (*res) = abs(d_x[(*index_one_based)-1]);
 }
 
 
 /*******************************************************************************
  *                         Private functions                                   *
  *******************************************************************************/
 
 /*
  *  out[j] = in[i], where j = ind[i] for i in [0,n-1]
  */
 void scatter(float       *out,
              const float *in,
              const int   *ind,
              int          n) {
 
   int num_blocks = (n / THREADS_PER_BLOCK) + 1;
   scatter_kernel<<<num_blocks, THREADS_PER_BLOCK>>>(out, in, ind, n);
 }
 
 
 /*******************************************************************************
  *                          Thrust-related functions                           *
  *******************************************************************************/
 
 template<typename BinaryFunction>
 void Segmented_reduce(const int    *key_start,
                       int           number_of_keys,
                       int           num_segments,
                       const float  *values,
                       void           *buffer,
                       float        *result,
                       BinaryFunction  binary_op) {
  
   int num_nnz_rows;
 
  /*  Memory layout of buffer:
   *  [ m*sizeof(float) Bytes | m*sizeof(int) Bytes]
   *  where m = "number of rows"
   */
   float *intermediate_result = (float*) buffer; 
   int   *nnz_rows            = (int*) (&intermediate_result[num_segments]);
 
   thrust::pair<int*,float*> new_end;
   thrust::equal_to<int> binary_pred;
   
   new_end = thrust::reduce_by_key(thrust::device,
                                   key_start,
                                   key_start + number_of_keys,
                                   values,
                                   nnz_rows,
                                   intermediate_result,
                                   binary_pred,
                                   binary_op);
 
   num_nnz_rows = new_end.first - nnz_rows;
   checkCudaErrors(cudaMemset(result, 0, num_segments * sizeof(float)));
   scatter(result, intermediate_result, nnz_rows, num_nnz_rows);
 }
 
 template<typename T>
 struct abs_maximum {
   typedef T first_argument_type;
   typedef T second_argument_type;
   typedef T result_type;
   __host__ __device__ T operator()(const T &lhs, const T &rhs) const {return max(abs(lhs), abs(rhs));}
  };
 
 template void Segmented_reduce<abs_maximum<float>>(const int          *key_start,
                                                      int                 number_of_keys,
                                                      int                 number_of_segments,
                                                      const float        *values,
                                                      void                 *buffer,
                                                      float              *result,
                                                      abs_maximum<float>  binary_op);
 
 
 /*******************************************************************************
  *                           API Functions                                     *
  *******************************************************************************/
 
 void cuda_vec_copy_d2d(float       *d_y,
                        const float *d_x,
                        int          n) {
 
   checkCudaErrors(cudaMemcpy(d_y, d_x, n * sizeof(float), cudaMemcpyDeviceToDevice));
 }
 
 void cuda_vec_copy_h2d(float       *d_y,
                        const float *h_x,
                        int          n) {
 
   checkCudaErrors(cudaMemcpy(d_y, h_x, n * sizeof(float), cudaMemcpyHostToDevice));
 }
 
 void cuda_vec_copy_d2h(float       *h_y,
                        const float *d_x,
                        int          n) {
 
   checkCudaErrors(cudaMemcpy(h_y, d_x, n * sizeof(float), cudaMemcpyDeviceToHost));
 }
 
 void cuda_veint_copy_h2d(int       *d_y,
                            const int *h_x,
                            int        n) {
 
   checkCudaErrors(cudaMemcpy(d_y, h_x, n * sizeof(int), cudaMemcpyHostToDevice));
 }
 
 void cuda_veint_copy_d2h(int       *h_y,
                            const int *d_x,
                            int        n) {
 
   checkCudaErrors(cudaMemcpy(h_y, d_x, n * sizeof(int), cudaMemcpyDeviceToHost));
 }
 
 void cuda_vec_set_sc(float *d_a,
                      float  sc,
                      int    n) {
 
   int number_of_blocks = (n / THREADS_PER_BLOCK) + 1;
   vec_set_sc_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(d_a, sc, n);
 }
 
 void cuda_vec_set_sc_cond(float     *d_a,
                           const int *d_test,
                           float      sc_if_neg,
                           float      sc_if_zero,
                           float      sc_if_pos,
                           float      n) {
 
   int number_of_blocks = (n / THREADS_PER_BLOCK) + 1;
 
   vec_set_sc_cond_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(d_a, d_test, sc_if_neg, sc_if_zero, sc_if_pos, n);
 }
 
 void cuda_vec_mult_sc(float *d_a,
                       float  sc,
                       int    n) {
 
   checkCudaErrors(cublasTscal(CUDA_handle->cublasHandle, n, &sc, d_a, 1));
 }
 
 void cuda_vec_add_scaled(float       *d_x,
                          const float *d_a,
                          const float *d_b,
                          float        sca,
                          float        scb,
                          int          n) {
 
   if (d_x != d_a || sca != 1.0) {
     if (sca == 1.0) {
       /* d_x = d_a */
       checkCudaErrors(cudaMemcpy(d_x, d_a, n * sizeof(float), cudaMemcpyDeviceToDevice));
     }
     else if (d_x == d_a) {
       /* d_x *= sca */
       checkCudaErrors(cublasTscal(CUDA_handle->cublasHandle, n, &sca, d_x, 1));
     }
     else {
       /* d_x = 0 */
       checkCudaErrors(cudaMemset(d_x, 0, n * sizeof(float)));
 
       /* d_x += sca * d_a */
       checkCudaErrors(cublasTaxpy(CUDA_handle->cublasHandle, n, &sca, d_a, 1, d_x, 1));
     }
   }
 
   /* d_x += scb * d_b */
   checkCudaErrors(cublasTaxpy(CUDA_handle->cublasHandle, n, &scb, d_b, 1, d_x, 1));
 }
 
 void cuda_vec_add_scaled3(float       *d_x,
                           const float *d_a,
                           const float *d_b,
                           const float *d_c,
                           float        sca,
                           float        scb,
                           float        scc,
                           int          n) {
 
   if (d_x != d_a || sca != 1.0) {
     if (sca == 1.0) {
       /* d_x = d_a */
       checkCudaErrors(cudaMemcpy(d_x, d_a, n * sizeof(float), cudaMemcpyDeviceToDevice));
     }
     else if (d_x == d_a) {
       /* d_x *= sca */
       checkCudaErrors(cublasTscal(CUDA_handle->cublasHandle, n, &sca, d_x, 1));
     }
     else {
       /* d_x = 0 */
       checkCudaErrors(cudaMemset(d_x, 0, n * sizeof(float)));
 
       /* d_x += sca * d_a */
       checkCudaErrors(cublasTaxpy(CUDA_handle->cublasHandle, n, &sca, d_a, 1, d_x, 1));
     }
   }
 
   /* d_x += scb * d_b */
   checkCudaErrors(cublasTaxpy(CUDA_handle->cublasHandle, n, &scb, d_b, 1, d_x, 1));
 
   /* d_x += scc * d_c */
   checkCudaErrors(cublasTaxpy(CUDA_handle->cublasHandle, n, &scc, d_c, 1, d_x, 1));
 }
 
 void cuda_vec_norm_inf(const float *d_x,
                        int          n,
                        float       *h_res) {
 
   cublasPointerMode_t mode;
   checkCudaErrors(cublasGetPointerMode(CUDA_handle->cublasHandle, &mode));
 
   if (mode == CUBLAS_POINTER_MODE_DEVICE) {
     checkCudaErrors(cublasITamax(CUDA_handle->cublasHandle, n, d_x, 1, CUDA_handle->d_index));
     abs_kernel<<<1,1>>>(CUDA_handle->d_index, d_x, h_res);  /* d_res actually */
   }
   else {
     int idx;
     checkCudaErrors(cublasITamax(CUDA_handle->cublasHandle, n, d_x, 1, &idx));
     checkCudaErrors(cudaMemcpy(h_res, d_x + (idx-1), sizeof(float), cudaMemcpyDeviceToHost));
     (*h_res) = abs(*h_res);
   }
 }
 
 void cuda_vec_norm_1(const float *d_x,
                      int          n,
                      float       *h_res) {
 
   cublasTasum(CUDA_handle->cublasHandle, n, d_x, 1, h_res);
 }
 
 void cuda_vec_norm_2(const float *d_x,
                      int          n,
                      float       *h_res) {
 
   cublasTnrm2(CUDA_handle->cublasHandle, n, d_x, 1, h_res);
 }
 
 void cuda_vec_scaled_norm_inf(const float *d_S,
                               const float *d_v,
                               int          n,
                               float       *h_res) {
 
   float *d_v_scaled;
 
   cuda_malloc((void **) &d_v_scaled, n * sizeof(float));
 
   /* d_v_scaled = d_S * d_v */
   cuda_vec_ew_prod(d_v_scaled, d_S, d_v, n);
 
   /* (*h_res) = |d_v_scaled|_inf */
   cuda_vec_norm_inf(d_v_scaled, n, h_res);
 
   cuda_free((void **) &d_v_scaled);
 }
 
 void cuda_vec_diff_norm_inf(const float *d_a,
                             const float *d_b,
                             int          n,
                             float       *h_res) {
 
   float *d_diff;
 
   cuda_malloc((void **) &d_diff, n * sizeof(float));
 
   /* d_diff = d_a - d_b */
   cuda_vec_add_scaled(d_diff, d_a, d_b, 1.0, -1.0, n);
 
   /* (*h_res) = |d_diff|_inf */
   cuda_vec_norm_inf(d_diff, n, h_res);
 
   cuda_free((void **) &d_diff);
 }
 
 void cuda_vec_mean(const float *d_x,
                    int          n,
                    float       *h_res) {
 
   cublasTasum(CUDA_handle->cublasHandle, n, d_x, 1, h_res);
   (*h_res) /= n;
 }
 
 void cuda_vec_prod(const float *d_a,
                    const float *d_b,
                    int          n,
                    float       *h_res) {
 
   checkCudaErrors(cublasTdot(CUDA_handle->cublasHandle, n, d_a, 1, d_b, 1, h_res));
 }
 
 void cuda_vec_prod_signed(const float *d_a,
                           const float *d_b,
                           int          sign,
                           int          n,
                           float       *h_res) {
 
   float *d_res;
   int number_of_blocks = (n / THREADS_PER_BLOCK) + 1;
 
   cuda_calloc((void **) &d_res, sizeof(float));
 
   if (sign == 1) {
     vec_prod_pos_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(d_a, d_b, d_res, n);
     checkCudaErrors(cudaMemcpy(h_res, d_res, sizeof(float), cudaMemcpyDeviceToHost));
   }
   else if (sign == -1) {
     vec_prod_neg_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(d_a, d_b, d_res, n);
     checkCudaErrors(cudaMemcpy(h_res, d_res, sizeof(float), cudaMemcpyDeviceToHost));
   }
   else {
     checkCudaErrors(cublasTdot(CUDA_handle->cublasHandle, n, d_a, 1, d_b, 1, h_res));
   }
 
   cuda_free((void **) &d_res);
 }
 
 void cuda_vec_ew_prod(float       *d_c,
                       const float *d_a,
                       const float *d_b,
                       int          n) {
 
   int number_of_blocks = (n / THREADS_PER_BLOCK) + 1;
 
   vec_ew_prod_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(d_c, d_a, d_b, n);
 }
 
 void cuda_vec_leq(const float *d_l,
                    const float *d_u,
                    int          n,
                    int         *h_res) {
 
   int *d_res;
   int number_of_blocks = (n / THREADS_PER_BLOCK) + 1;
 
   cuda_malloc((void **) &d_res, sizeof(int));
 
   /* Initialize d_res to 1 */
   *h_res = 1;
   checkCudaErrors(cudaMemcpy(d_res, h_res, sizeof(int), cudaMemcpyHostToDevice));
 
   vec_leq_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(d_l, d_u, d_res, n);
 
   checkCudaErrors(cudaMemcpy(h_res, d_res, sizeof(int), cudaMemcpyDeviceToHost));
 
   cuda_free((void **) &d_res);
 }
 
 void cuda_vec_bound(float       *d_x,
                     const float *d_z,
                     const float *d_l,
                     const float *d_u,
                     int          n) {
 
   int number_of_blocks = (n / THREADS_PER_BLOCK) + 1;
 
   vec_bound_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(d_x, d_z, d_l, d_u, n);
 }
 
 void cuda_vec_project_polar_reccone(float       *d_y,
                                     const float *d_l,
                                     const float *d_u,
                                     float        infval,
                                     int          n) {
 
   int number_of_blocks = (n / THREADS_PER_BLOCK) + 1;
 
   vec_project_polar_reccone_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(d_y, d_l, d_u, infval, n);
 }
 
 void cuda_vec_in_reccone(const float *d_y,
                          const float *d_l,
                          const float *d_u,
                          float        infval,
                          float        tol,
                          int          n,
                          int         *h_res) {
 
   int *d_res;
   int number_of_blocks = (n / THREADS_PER_BLOCK) + 1;
 
   cuda_malloc((void **) &d_res, sizeof(int));
 
   /* Initialize d_res to 1 */
   *h_res = 1;
   checkCudaErrors(cudaMemcpy(d_res, h_res, sizeof(int), cudaMemcpyHostToDevice));
 
   vec_in_reccone_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(d_y, d_l, d_u, infval, tol, d_res, n);
 
   checkCudaErrors(cudaMemcpy(h_res, d_res, sizeof(int), cudaMemcpyDeviceToHost));
 
   cuda_free((void **) &d_res);
 }
 
 void cuda_vec_reciprocal(float       *d_b,
                          const float *d_a,
                          int          n) {
 
   int number_of_blocks = (n / THREADS_PER_BLOCK) + 1;
 
   vec_reciprocal_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(d_b, d_a, n);
 }
 
 void cuda_vec_sqrt(float *d_a,
                    int    n) {
 
   int number_of_blocks = (n / THREADS_PER_BLOCK) + 1;
 
   vec_sqrt_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(d_a, n);
 }
 
 void cuda_vec_max(float       *d_c,
                   const float *d_a,
                   const float *d_b,
                   int          n) {
 
   int number_of_blocks = (n / THREADS_PER_BLOCK) + 1;
 
   vec_max_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(d_c, d_a, d_b, n);
 }
 
 void cuda_vec_min(float       *d_c,
                   const float *d_a,
                   const float *d_b,
                   int          n) {
 
   int number_of_blocks = (n / THREADS_PER_BLOCK) + 1;
 
   vec_min_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(d_c, d_a, d_b, n);
 }
 
 void cuda_vec_bounds_type(int         *d_iseq,
                           const float *d_l,
                           const float *d_u,
                           float        infval,
                           float        tol,
                           int          n,
                           int         *h_has_changed) {
 
   int *d_has_changed;
   int number_of_blocks = (n / THREADS_PER_BLOCK) + 1;
 
   /* Initialize d_has_changed to zero */
   cuda_calloc((void **) &d_has_changed, sizeof(int));
 
   vec_bounds_type_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(d_iseq, d_l, d_u, infval, tol, d_has_changed, n);
 
   checkCudaErrors(cudaMemcpy(h_has_changed, d_has_changed, sizeof(int), cudaMemcpyDeviceToHost));
 
   cuda_free((void **) &d_has_changed);
 }
 
 void cuda_vec_set_sc_if_lt(float       *d_x,
                            const float *d_z,
                            float        testval,
                            float        newval,
                            int          n) {
 
   int number_of_blocks = (n / THREADS_PER_BLOCK) + 1;
 
   vec_set_sc_if_lt_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(d_x, d_z, testval, newval, n);
 }
 
 void cuda_vec_set_sc_if_gt(float       *d_x,
                            const float *d_z,
                            float        testval,
                            float        newval,
                            int          n) {
 
   int number_of_blocks = (n / THREADS_PER_BLOCK) + 1;
 
   vec_set_sc_if_gt_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(d_x, d_z, testval, newval, n);
 }
 
 void cuda_vec_segmented_sum(const float *d_values,
                             const int   *d_keys,
                             float       *d_res,
                             void          *d_buffer,
                             int          num_segments,
                             int          num_elements) {
 
   thrust::plus<float> binary_op;
   Segmented_reduce(d_keys, num_elements, num_segments, d_values, d_buffer, d_res, binary_op);
 }
 
 void cuda_mat_mult_sc(csr     *S,
                       csr     *At,
                       int    symmetric,
                       float  sc) {
 
   checkCudaErrors(cublasTscal(CUDA_handle->cublasHandle, S->nnz, &sc, S->val, 1));
 
   if (!symmetric) {
     /* Update At as well */
     checkCudaErrors(cublasTscal(CUDA_handle->cublasHandle, At->nnz, &sc, At->val, 1));
   }
 }
 
 void cuda_mat_lmult_diag(csr           *S,
                          csr           *At,
                          int          symmetric,
                          const float *d_diag) {
 
   int nnz = S->nnz;
   int number_of_blocks = (nnz / THREADS_PER_BLOCK) / ELEMENTS_PER_THREAD + 1;
 
   mat_lmult_diag_kernel<<<number_of_blocks,THREADS_PER_BLOCK>>>(S->row_ind, d_diag, S->val, nnz);
 
   if (!symmetric) {
     /* Multiply At from right */
     mat_rmult_diag_kernel<<<number_of_blocks,THREADS_PER_BLOCK>>>(At->col_ind, d_diag, At->val, nnz);
   }
 }
 
 void cuda_mat_rmult_diag(csr           *S,
                          csr           *At,
                          int          symmetric,
                          const float *d_diag) {
 
   int nnz = S->nnz;
   int number_of_blocks = (nnz / THREADS_PER_BLOCK) / ELEMENTS_PER_THREAD + 1;
 
   mat_rmult_diag_kernel<<<number_of_blocks,THREADS_PER_BLOCK>>>(S->col_ind, d_diag, S->val, nnz);
 
   if (!symmetric) {
     /* Multiply At from left */
     mat_lmult_diag_kernel<<<number_of_blocks,THREADS_PER_BLOCK>>>(At->row_ind, d_diag, At->val, nnz);
   }
 }
 
 void cuda_mat_rmult_diag_new(const csr     *S,
                              float       *d_buffer,
                              const float *d_diag) {
 
   int nnz = S->nnz;
   int number_of_blocks = (nnz / THREADS_PER_BLOCK) / ELEMENTS_PER_THREAD + 1;
 
   mat_rmult_diag_new_kernel<<<number_of_blocks,THREADS_PER_BLOCK>>>(S->col_ind, d_diag, S->val, d_buffer, nnz);
 }
 
 void cuda_mat_Axpy(const csr     *A,
                    const float *d_x,
                    float       *d_y,
                    float        alpha,
                    float        beta) {
 
   if (A->nnz == 0 || alpha == 0.0) {
     /* d_y = beta * d_y */
     cuda_vec_mult_sc(d_y, beta, A->m);
     return;
   }
 
   checkCudaErrors(cusparseCsrmv(CUDA_handle->cusparseHandle, A->alg, A->m, A->n, A->nnz, &alpha, A->MatDescription, A->val, A->row_ptr, A->col_ind, d_x, &beta, d_y, A->buffer));
 }
 
 void cuda_mat_quad_form(const csr     *P,
                         const float *d_x,
                         float       *h_res) {
 
   int n = P->n;
   float *d_Px;
 
   cuda_malloc((void **) &d_Px, n * sizeof(float));
 
   /* d_Px = P * x */
   cuda_mat_Axpy(P, d_x, d_Px, 1.0, 0.0);
 
   /* h_res = d_Px' * d_x */
   cuda_vec_prod(d_Px, d_x, n, h_res);
 
   /* h_res *= 0.5 */
   (*h_res) *= 0.5;
 
   cuda_free((void **) &d_Px);
 }
 
 void cuda_mat_row_norm_inf(const csr *S,
                            float   *d_res) {
 
   int nnz      = S->nnz;
   int num_rows = S->m;
 
   if (nnz == 0) return;
 
   abs_maximum<float> binary_op;
   void *d_buffer;
   cuda_malloc(&d_buffer, num_rows * (sizeof(float) + sizeof(int)));
 
   /* 
   *  For rows with only one element, the element itself is returned.
   *  Therefore, we have to take the absolute value to get the inf-norm.
   */
   Segmented_reduce(S->row_ind, nnz, num_rows, S->val, d_buffer, d_res, binary_op);
   vec_abs_kernel<<<num_rows/THREADS_PER_BLOCK+1,THREADS_PER_BLOCK>>>(d_res, num_rows);
 
   cuda_free(&d_buffer);
 }