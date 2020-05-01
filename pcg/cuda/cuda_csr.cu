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

 #include "cuda_csr.h"
 #include "cuda_configure.h"
 #include "cuda_handler.h"
 #include "cuda_malloc.h"
 #include "cuda_wrapper.h"
 #include "helper_cuda.h"    /* --> checkCudaErrors */
 
 #include "csr_type.h"
 #include "glob_opts.h"
 
 #include <thrust/scan.h>
 #include <thrust/execution_policy.h>
  
 /* This function is implemented in cuda_lin_alg.cu */
 extern void scatter(float *out, const float *in, const int *ind, int n);
 
 
 /*******************************************************************************
  *                            GPU Kernels                                      *
  *******************************************************************************/
 
  /*
  * Expand an upper triangular matrix given in COO format to a symmetric
  * matrix. Each entry is duplicated with its column- and row index switched.
  * In the case of a diagonal element we set the indices to a value  that is
  * larger than n to easily remove it later. This is done to keep the memory
  * patern one to one (MAP operation).
  * 
  * Additionally, it adds additional n diagonal elements to have a full 
  * diagonal.
  * 
  * The output arrays row_ind_out and col_ind_out have to be of size 2*nnz+n.
  */
 __global__ void fill_full_matrix_kernel(int       *row_ind_out,
                                         int       *col_ind_out,
                                         int       *nnz_on_diag,
                                         int       *has_non_zero_diag_element,
                                         const int *__restrict__ row_ind_in,
                                         const int *__restrict__ col_ind_in,
                                         int        nnz,
                                         int        n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < nnz; i += grid_size) {
     int row = row_ind_in[i];
     int column = col_ind_in[i];
 
     row_ind_out[i] = row;
     col_ind_out[i] = column;
 
     if (row == column) {
       has_non_zero_diag_element[row] = 1;
       row_ind_out[i + nnz] = column + n; /* dummy value for sorting and removal later on */
       col_ind_out[i + nnz] = row + n;
       atomicAdd(nnz_on_diag, 1);
     }
     else {
       row_ind_out[i + nnz] = column;
       col_ind_out[i + nnz] = row;
     }
   }
 }
 
 /**
  * Insert elements at structural zeros on the diagonal of the sparse matrix
  * specified by row and column index (COO format). To keep a one-to-one memory
  * patern we add n new elements to the matrix. In case where there already is a
  * diagonal element we add a dummy entry. The dummy entries will be removed later.
  */
 __global__ void add_diagonal_kernel(int       *row_ind,
                                     int       *col_ind,
                                     const int *has_non_zero_diag_element,
                                     int        n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int row = idx; row < n; row += grid_size) {
     if (has_non_zero_diag_element[row] == 0) {
       row_ind[row] = row; 
       col_ind[row] = row;
     }
     else {
       row_ind[row] = row + n; /* dummy value, for easy removal after sorting */
       col_ind[row] = row + n;
     }
   }
 }
 
 /*
  * Permutation in: (size n, range 2*nnz+n):
  * 
  * Gathers from the following array to create the full matrix :
  * 
  *       |P_lower->val|P_lower->val|zeros(n)|
  *
  *       
  * Permutation out: (size n, range new_range)
  * 
  * Gathers from the following array to create the full matrix :
  * 
  *          |P_lower->val|zeros(1)|
  *                             
  *          | x[i] mod new_range    if x[i] <  2 * new_range
  * x[i] ->  | new_range             if x[i] >= 2 * new_range   
  * 
  */
 __global__ void reduce_permutation_kernel(int *permutation,
                                           int  new_range,
                                           int  n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for(int i = idx; i < n; i += grid_size) {
     if (permutation[i] < 2 * new_range) {
       permutation[i] = permutation[i] % new_range;
     }
     else {
       permutation[i] = new_range; /* gets the 0 element at nnz+1 of the value array */
     }
   }
 }
 
 __global__ void get_diagonal_indices_kernel(int *row_ind,
                                             int *col_ind,
                                             int  nnz,
                                             int *diag_index) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_size = blockDim.x * gridDim.x;
 
   for (int index = idx; index < nnz; index += grid_size) {
     int row = row_ind[index];
     int column = col_ind[index];
 
     if (row == column) {
       diag_index[row] = index;
     }
   }
 }
 
 __global__ void predicate_generator_kernel(const int *row_ind,
                                            const int *row_predicate,
                                            int       *predicate,
                                            int        nnz) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
   int grid_stride = gridDim.x * blockDim.x;
 
   for(int i = idx; i < nnz; i += grid_stride) {
     int row = row_ind[i];
     predicate[i] = row_predicate[row];
   }
 }
 
 template<typename T>
 __global__ void compact(const T *data_in,
                         T       *data_out,
                         int   *predicate,
                         int   *scatter_addres,
                         int    n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
 
   if(idx < n) {
     if(predicate[idx]) {
       int write_ind = scatter_addres[idx] - 1;
       data_out[write_ind] = data_in[idx];
     }
   }
 }
 
 __global__ void compact_rows(const int *row_ind,
                              int       *data_out,
                              int       *new_row_number,
                              int       *predicate,
                              int       *scatter_addres,
                              int        n) {
 
   int idx = threadIdx.x + blockDim.x * blockIdx.x;
 
   if(idx < n) {
     if(predicate[idx]) {
       int write_ind = scatter_addres[idx] - 1;
       int row = row_ind[idx];
       data_out[write_ind] = new_row_number[row]-1;
     }
   }
 }
 
 __global__ void vector_init_abs_kernel(const int *a,
                                        int       *b,
                                        int        n) {
 
   int i  = threadIdx.x + blockDim.x * blockIdx.x;
 
   if (i < n) {
     b[i] = abs(a[i]);
   }
 }
 
 
 /*******************************************************************************
  *                         Private Functions                                   *
  *******************************************************************************/
 
  /*
  *  Update the size of buffer used for the merge path based
  *  sparse matrix-vector product (spmv).
  */
 void update_mp_buffer(csr *P) {
 
   size_t bufferSizeInBytes = 0;
   float alpha = 1.0;
 
   checkCudaErrors(cusparseCsrmv_bufferSize(CUDA_handle->cusparseHandle,
                                            P->alg, P->m, P->n, P->nnz,
                                            &alpha,
                                            P->MatDescription, P->val, P->row_ptr, P->col_ind,
                                            NULL,
                                            &alpha,
                                            NULL,
                                            &bufferSizeInBytes));
   
   if (bufferSizeInBytes > P->bufferSizeInBytes) {
     cuda_free((void **) &P->buffer);                            
     cuda_malloc((void **) &P->buffer, bufferSizeInBytes);
     P->bufferSizeInBytes = bufferSizeInBytes;
   }
 }
 
  /*
  *  Creates a CSR matrix with the specified dimension (m,n,nnz).
  *  
  *  If specified, it allocates proper amount of device memory
  *  allocate_on_device = 1: device memory for CSR
  *  allocate_on_device = 2: device memory for CSR (+ col_ind)  
  */
 csr* csr_alloc(int m,
                int n,
                int nnz,
                int allocate_on_device) {
 
   csr *dev_mat = (csr*) calloc(1, sizeof(csr));
 
   if (!dev_mat) return NULL;
 
   dev_mat->m   = m;
   dev_mat->n   = n;
   dev_mat->nnz = nnz;
       
//  #ifdef IS_WINDOWS
//    /* MERGE_PATH is not working properly on WINDOWS */
//    dev_mat->alg = CUSPARSE_ALG_NAIVE;
//  #else
   dev_mat->alg = CUSPARSE_ALG_MERGE_PATH;
//  #endif
 
   dev_mat->buffer = NULL;
   dev_mat->bufferSizeInBytes = 0;
 
   checkCudaErrors(cusparseCreateMatDescr(&dev_mat->MatDescription));
   cusparseSetMatType(dev_mat->MatDescription, CUSPARSE_MATRIX_TYPE_GENERAL);
   cusparseSetMatIndexBase(dev_mat->MatDescription, CUSPARSE_INDEX_BASE_ZERO);
 
   if (allocate_on_device > 0) {
     cuda_calloc((void **) &dev_mat->val, (dev_mat->nnz + 1) * sizeof(float));
     cuda_malloc((void **) &dev_mat->row_ptr, (dev_mat->m + 1) * sizeof(int)); 
     cuda_malloc((void **) &dev_mat->col_ind, dev_mat->nnz * sizeof(int));
 
     if (allocate_on_device > 1) {
       cuda_malloc((void **) &dev_mat->row_ind, dev_mat->nnz * sizeof(int));
     } 
   }
   return dev_mat;
 }
 
 /*
  *  Copy CSR matrix from host to device.
  *  The device memory should be pre-allocated.
  */
 void csr_copy_h2d(csr           *dev_mat,
                   const int   *h_row_ptr,
                   const int   *h_col_ind,
                   const float *h_val) {
 
   checkCudaErrors(cudaMemcpy(dev_mat->row_ptr, h_row_ptr, (dev_mat->m + 1) * sizeof(int), cudaMemcpyHostToDevice));
   checkCudaErrors(cudaMemcpy(dev_mat->col_ind, h_col_ind, dev_mat->nnz * sizeof(int), cudaMemcpyHostToDevice));
   checkCudaErrors(cudaMemcpy(dev_mat->val, h_val, dev_mat->nnz * sizeof(float), cudaMemcpyHostToDevice));
 }
 
 csr* csr_init(int          m,
               int          n,
               const int   *h_row_ptr,
               const int   *h_col_ind,
               const float *h_val) {
     
   csr *dev_mat = csr_alloc(m, n, h_row_ptr[m], 1);
   
   if (!dev_mat) return NULL;
   
   if (m == 0) return dev_mat;
 
   /* copy_matrix_to_device */
   csr_copy_h2d(dev_mat, h_row_ptr, h_col_ind, h_val);
   update_mp_buffer(dev_mat);
 
   return dev_mat;
 }
 
 /*
  *  Compress row indices from the COO format to the row pointer
  *  of the CSR format.
  */
 void compress_row_ind(csr *mat) {
 
   cuda_free((void** ) &mat->row_ptr);
   cuda_malloc((void** ) &mat->row_ptr, (mat->m + 1) * sizeof(float));
   checkCudaErrors(cusparseXcoo2csr(CUDA_handle->cusparseHandle, mat->row_ind, mat->nnz, mat->m, mat->row_ptr, CUSPARSE_INDEX_BASE_ZERO));
 }
 
 void csr_expand_row_ind(csr *mat) {
 
   if (!mat->row_ind) {
     cuda_malloc((void** ) &mat->row_ind, mat->nnz * sizeof(float));
     checkCudaErrors(cusparseXcsr2coo(CUDA_handle->cusparseHandle, mat->row_ptr, mat->nnz, mat->m, mat->row_ind, CUSPARSE_INDEX_BASE_ZERO));
   }
 }
 
 /*
  *  Sorts matrix in COO format by row. It returns a permutation
  *  vector that describes reordering of the elements.
  */
 int* coo_sort(csr *A) {
 
   int *A_to_At_permutation;
   char *pBuffer;
   size_t pBufferSizeInBytes;
 
   cuda_malloc((void **) &A_to_At_permutation, A->nnz * sizeof(int));
   checkCudaErrors(cusparseCreateIdentityPermutation(CUDA_handle->cusparseHandle, A->nnz, A_to_At_permutation));
 
   checkCudaErrors(cusparseXcoosort_bufferSizeExt(CUDA_handle->cusparseHandle, A->m, A->n, A->nnz, A->row_ind, A->col_ind, &pBufferSizeInBytes));
 
   cuda_malloc((void **) &pBuffer, pBufferSizeInBytes * sizeof(char));
 
   checkCudaErrors(cusparseXcoosortByRow(CUDA_handle->cusparseHandle, A->m, A->n, A->nnz, A->row_ind, A->col_ind, A_to_At_permutation, pBuffer));
 
   cuda_free((void **) &pBuffer);
 
   return A_to_At_permutation;
 }
 
 /*
  * Compute transpose of a matrix in COO format.
  */
 void coo_tranpose(csr* A) {
   int m = A->m;
   A->m = A->n;
   A->n = m;
 
   int *row_ind = A->row_ind;
   A->row_ind = A->col_ind;
   A->col_ind = row_ind;
 }
 
 /*
  *  values[i] = values[permutation[i]] for i in [0,n-1]
  */
 void permute_vector(float     *values,
                     const int *permutation,
                     int        n) {
 
   float *permuted_values;
   cuda_malloc((void **) &permuted_values, n * sizeof(float));
 
   checkCudaErrors(cusparseTgthr(CUDA_handle->cusparseHandle, n, values, permuted_values, permutation, CUSPARSE_INDEX_BASE_ZERO));
 
   checkCudaErrors(cudaMemcpy(values, permuted_values, n * sizeof(float), cudaMemcpyDeviceToDevice));
   cuda_free((void **) &permuted_values);
 }
 
 /*
  *  target[i] = source[permutation[i]] for i in [0,n-1]
  *  
  *  target and source cannot point to the same location
  */
 void permute_vector(float       *target,
                     const float *source,
                     const int   *permutation,
                     int          n) {
 
   checkCudaErrors(cusparseTgthr(CUDA_handle->cusparseHandle, n, source, target, permutation, CUSPARSE_INDEX_BASE_ZERO));
 }
 
 /*
  *  Copy the values and pointers form target to the source matrix.
  *  The device memory of source has to be freed first to avoid a
  *  memory leak in case it holds allocated memory.
  *  
  *  The MatrixDescription has to be destroyed first since it is a
  *  pointer hidded by a typedef.
  *  
  *  The pointers of source matrix are set to NULL to avoid
  *  accidental freeing of the associated memory blocks.
  */
 void copy_csr(csr* target,
               csr* source) {
 
   target->m                 = source->m;
   target->n                 = source->n;
   target->nnz               = source->nnz;
   target->bufferSizeInBytes = source->bufferSizeInBytes;
   target->alg               = source->alg;
 
   cusparseDestroyMatDescr(target->MatDescription);
   cuda_free((void **) &target->val);
   cuda_free((void **) &target->row_ind);
   cuda_free((void **) &target->row_ptr);
   cuda_free((void **) &target->col_ind);
   cuda_free((void **) &target->buffer);
 
   target->val            = source->val;
   target->row_ind        = source->row_ind;
   target->row_ptr        = source->row_ptr;
   target->col_ind        = source->col_ind;
   target->buffer         = source->buffer;
   target->MatDescription = source->MatDescription; 
 
   source->val            = NULL;
   source->row_ind        = NULL;
   source->row_ptr        = NULL;
   source->col_ind        = NULL;
   source->buffer         = NULL;
   source->MatDescription = NULL;
 }
 
 void csr_triu_to_full(csr    *P_triu,
                       int **P_triu_to_full_permutation,
                       int **P_diag_indices) {
 
   int number_of_blocks;
   int *has_non_zero_diag_element, *d_nnz_diag;
   int h_nnz_diag, Full_nnz, nnz_triu, n, nnz_max_Full;
   int offset;
 
   nnz_triu     = P_triu->nnz;
   n            = P_triu->n;
   nnz_max_Full = 2*nnz_triu + n;
 
   csr *Full_P = csr_alloc(n, n, nnz_max_Full, 2);
   cuda_calloc((void **) &has_non_zero_diag_element, n * sizeof(int));
   cuda_calloc((void **) &d_nnz_diag, sizeof(int));
 
   csr_expand_row_ind(P_triu);
 
   number_of_blocks = (nnz_triu / THREADS_PER_BLOCK) + 1;
   fill_full_matrix_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(Full_P->row_ind, Full_P->col_ind, d_nnz_diag, has_non_zero_diag_element, P_triu->row_ind, P_triu->col_ind, nnz_triu, n);
 
   offset = 2 * nnz_triu;
   number_of_blocks = (n / THREADS_PER_BLOCK) + 1;
   add_diagonal_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(Full_P->row_ind + offset, Full_P->col_ind + offset, has_non_zero_diag_element, n);
 
   /* The Full matrix now is of size (2n)x(2n)
     *                  [P 0]
     *                  [0 D]
     * where P is the desired full matrix and D is
     * a diagonal that contains dummy values
   */
   
   checkCudaErrors(cudaMemcpy(&h_nnz_diag, d_nnz_diag, sizeof(int), cudaMemcpyDeviceToHost));
 
   Full_nnz = (2 * (nnz_triu - h_nnz_diag)) + n;
   int *d_P = coo_sort(Full_P);
 
   number_of_blocks = (nnz_triu / THREADS_PER_BLOCK) + 1;
   reduce_permutation_kernel<<<number_of_blocks,THREADS_PER_BLOCK>>>(d_P, nnz_triu, Full_nnz);
 
   permute_vector(Full_P->val, P_triu->val, d_P, Full_nnz);
 
   cuda_malloc((void **) P_triu_to_full_permutation, Full_nnz * sizeof(int));
   checkCudaErrors(cudaMemcpy(*P_triu_to_full_permutation, d_P, Full_nnz * sizeof(int), cudaMemcpyDeviceToDevice));
   cuda_malloc((void **) P_diag_indices, n * sizeof(int));
 
   number_of_blocks = (Full_nnz / THREADS_PER_BLOCK) + 1;
   get_diagonal_indices_kernel<<<number_of_blocks, THREADS_PER_BLOCK>>>(Full_P->row_ind, Full_P->col_ind, Full_nnz, *P_diag_indices);
 
   Full_P->nnz = Full_nnz;
   compress_row_ind(Full_P);
   update_mp_buffer(Full_P); 
   copy_csr(P_triu, Full_P);
 
   cuda_mat_free(Full_P);
   cuda_free((void **) &d_P);
   cuda_free((void **) &d_nnz_diag);
   cuda_free((void **) &has_non_zero_diag_element);
 }
 
 /**
  * Matrix A is converted from CSC to CSR. The data in A is interpreted as
  * being in CSC format, even if it is in CSR.
  * This operation is equivalent to a transpose. We temporarily allocate space
  * for the new matrix since this operation cannot be done inplace.
  * Additionally, a gather indices vector is generated to perform the conversion
  * from A to A' faster during a matrix update.
  */
 void csr_transpose(csr    *A,
                    int **A_to_At_permutation) {
 
   (*A_to_At_permutation) = NULL;
 
   if (A->nnz == 0) {
     int tmp = A->n;
     A->n = A->m;
     A->m = tmp;
     return;
   }
 
   csr_expand_row_ind(A);
   coo_tranpose(A);
   (*A_to_At_permutation) = coo_sort(A);
   compress_row_ind(A);
 
   permute_vector(A->val, *A_to_At_permutation, A->nnz);
 
   update_mp_buffer(A);
 }
 
 
 /*******************************************************************************
  *                           API Functions                                     *
  *******************************************************************************/
 
 void cuda_mat_init_P(const csc  *mat,
                      csr       **P,
                      float   **d_P_triu_val,
                      int     **d_P_triu_to_full_ind,
                      int     **d_P_diag_ind) {
 
   int n   = mat->n;
   int nnz = mat->p[n];
   
   /* Initialize upper triangular part of P */
   *P = csr_init(n, n, mat->p, mat->i, mat->x);
 
   /* Convert P to a full matrix. Store indices of diagonal and triu elements. */
   csr_triu_to_full(*P, d_P_triu_to_full_ind, d_P_diag_ind);
   csr_expand_row_ind(*P);
 
   /* We need 0.0 at val[nzz] -> nnz+1 elements */
   cuda_calloc((void **) d_P_triu_val, (nnz+1) * sizeof(float));
 
   /* Store triu elements */
   checkCudaErrors(cudaMemcpy(*d_P_triu_val, mat->x, nnz * sizeof(float), cudaMemcpyHostToDevice));
 }
 
 void cuda_mat_init_A(const csc  *mat,
                      csr       **A,
                      csr       **At,
                      int     **d_A_to_At_ind) {
 
   int m = mat->m;
   int n = mat->n;
 
   /* Initializing At is easy since it is equal to A in CSC */
   *At = csr_init(n, m, mat->p, mat->i, mat->x);
   csr_expand_row_ind(*At);
 
   /* We need to take transpose of At to get A */
   *A = csr_init(n, m, mat->p, mat->i, mat->x);
   csr_transpose(*A, d_A_to_At_ind);
   csr_expand_row_ind(*A);
 }
 
 void cuda_mat_update_P(const float  *Px,
                        const int    *Px_idx,
                        int           Px_n,
                        csr           **P,
                        float        *d_P_triu_val,
                        int          *d_P_triu_to_full_ind,
                        int          *d_P_diag_ind,
                        int           P_triu_nnz) {
 
   if (!Px_idx) { /* Update whole P */
     float *d_P_val_new;
 
     /* Allocate memory */
     cuda_malloc((void **) &d_P_val_new, (P_triu_nnz + 1) * sizeof(float));
 
     /* Copy new values from host to device */
     checkCudaErrors(cudaMemcpy(d_P_val_new, Px, P_triu_nnz * sizeof(float), cudaMemcpyHostToDevice));
 
     checkCudaErrors(cusparseTgthr(CUDA_handle->cusparseHandle, (*P)->nnz, d_P_val_new, (*P)->val, d_P_triu_to_full_ind, CUSPARSE_INDEX_BASE_ZERO));
 
     cuda_free((void **) &d_P_val_new);
   }
   else { /* Update P partially */
     float *d_P_val_new;
     int   *d_P_ind_new;
 
     /* Allocate memory */
     cuda_malloc((void **) &d_P_val_new, Px_n * sizeof(float));
     cuda_malloc((void **) &d_P_ind_new, Px_n * sizeof(int));
 
     /* Copy new values and indices from host to device */
     checkCudaErrors(cudaMemcpy(d_P_val_new, Px,     Px_n * sizeof(float), cudaMemcpyHostToDevice));
     checkCudaErrors(cudaMemcpy(d_P_ind_new, Px_idx, Px_n * sizeof(int),   cudaMemcpyHostToDevice));
 
     /* Update d_P_triu_val */
     scatter(d_P_triu_val, d_P_val_new, d_P_ind_new, Px_n);
 
     /* Gather from d_P_triu_val to update full P */
     checkCudaErrors(cusparseTgthr(CUDA_handle->cusparseHandle, (*P)->nnz, d_P_triu_val, (*P)->val, d_P_triu_to_full_ind, CUSPARSE_INDEX_BASE_ZERO));
 
     cuda_free((void **) &d_P_val_new);
     cuda_free((void **) &d_P_ind_new);
   }
 }
 
 void cuda_mat_update_A(const float  *Ax,
                        const int    *Ax_idx,
                        int           Ax_n,
                        csr           **A,
                        csr           **At,
                        int          *d_A_to_At_ind) {
 
   int Annz     = (*A)->nnz;
   float *Aval  = (*A)->val;
   float *Atval = (*At)->val;
 
   if (!Ax_idx) { /* Update whole A */
     /* Updating At is easy since it is equal to A in CSC */
     checkCudaErrors(cudaMemcpy(Atval, Ax, Annz * sizeof(float), cudaMemcpyHostToDevice));
 
     /* Updating A requires transpose of A_new */
     checkCudaErrors(cusparseTgthr(CUDA_handle->cusparseHandle, Annz, Atval, Aval, d_A_to_At_ind, CUSPARSE_INDEX_BASE_ZERO));
   }
   else { /* Update A partially */
     float *d_At_val_new;
     int   *d_At_ind_new;
 
     /* Allocate memory */
     cuda_malloc((void **) &d_At_val_new, Ax_n * sizeof(float));
     cuda_malloc((void **) &d_At_ind_new, Ax_n * sizeof(int));
 
     /* Copy new values and indices from host to device */
     checkCudaErrors(cudaMemcpy(d_At_val_new, Ax,     Ax_n * sizeof(float), cudaMemcpyHostToDevice));
     checkCudaErrors(cudaMemcpy(d_At_ind_new, Ax_idx, Ax_n * sizeof(int),   cudaMemcpyHostToDevice));
 
     /* Update At first since it is equal to A in CSC */
     scatter(Atval, d_At_val_new, d_At_ind_new, Ax_n);
 
     cuda_free((void **) &d_At_val_new);
     cuda_free((void **) &d_At_ind_new);
 
     /* Gather from Atval to construct Aval */
     checkCudaErrors(cusparseTgthr(CUDA_handle->cusparseHandle, Annz, Atval, Aval, d_A_to_At_ind, CUSPARSE_INDEX_BASE_ZERO));
   }
 }
 
 void cuda_mat_free(csr *mat) {
   if (mat) {
     cuda_free((void **) &mat->val);
     cuda_free((void **) &mat->row_ptr);
     cuda_free((void **) &mat->col_ind);
     cuda_free((void **) &mat->buffer);
     cuda_free((void **) &mat->row_ind);
     cusparseDestroyMatDescr(mat->MatDescription);
     free(mat);
   }
 }
 
 void cuda_submat_byrows(const csr    *A,
                         const int  *d_rows,
                         csr         **Ared,
                         csr         **Aredt) {
 
   int new_m = 0;
 
   int n   = A->n;
   int m   = A->m;
   int nnz = A->nnz;
 
   int *d_predicate;
   int *d_compact_address;
   int *d_row_predicate;
   int *d_new_row_number;
 
   cuda_malloc((void **) &d_row_predicate,  m * sizeof(int));
   cuda_malloc((void **) &d_new_row_number, m * sizeof(int));
 
   cuda_malloc((void **) &d_predicate,       nnz * sizeof(int));
   cuda_malloc((void **) &d_compact_address, nnz * sizeof(int));
 
   // Copy rows array to device and set -1s to ones
   checkCudaErrors(cudaMemcpy(d_row_predicate, d_rows, m * sizeof(int), cudaMemcpyDeviceToDevice));
   vector_init_abs_kernel<<<(m/THREADS_PER_BLOCK) + 1,THREADS_PER_BLOCK>>>(d_row_predicate, d_row_predicate, m);
 
   // Calculate new row numbering and get new number of rows
   thrust::inclusive_scan(thrust::device, d_row_predicate, d_row_predicate + m, d_new_row_number);
   if (m) {
     checkCudaErrors(cudaMemcpy(&new_m, &d_new_row_number[m-1], sizeof(int), cudaMemcpyDeviceToHost));
   }
   else {
     (*Ared) = (csr *) calloc(1, sizeof(csr));
     (*Ared)->n = n;
 
     (*Aredt) = (csr *) calloc(1, sizeof(csr));
     (*Aredt)->m = n;
 
     return;
   }
 
   // Generate predicates per element from per row predicate
   predicate_generator_kernel<<<(nnz/THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK>>>(A->row_ind, d_row_predicate, d_predicate, nnz);
 
   // Get array offset for compacting and new nnz
   thrust::inclusive_scan(thrust::device, d_predicate, d_predicate + nnz, d_compact_address);
   int nnz_new;
   if (nnz) checkCudaErrors(cudaMemcpy(&nnz_new, &d_compact_address[nnz-1], sizeof(int), cudaMemcpyDeviceToHost));
 
   // allocate new matrix (2 -> allocate row indices as well)
   (*Ared) = csr_alloc(new_m, n, nnz_new, 2);
 
   // Compact arrays according to given predicates, special care has to be taken for the rows
   compact_rows<<<(nnz/THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK>>>(A->row_ind, (*Ared)->row_ind, d_new_row_number, d_predicate, d_compact_address, nnz);
   compact<<<(nnz/THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK>>>(A->col_ind, (*Ared)->col_ind, d_predicate, d_compact_address, nnz);
   compact<<<(nnz/THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK>>>(A->val, (*Ared)->val, d_predicate, d_compact_address, nnz);
 
   // Generate row pointer
   compress_row_ind(*Ared);
 
   // Update merge path buffer (CsrmvEx)
   update_mp_buffer(*Ared);
 
   // We first make a copy of Ared
   *Aredt = csr_alloc(new_m, n, nnz_new, 1);
   checkCudaErrors(cudaMemcpy((*Aredt)->val,     (*Ared)->val,     nnz_new   * sizeof(float), cudaMemcpyDeviceToDevice));
   checkCudaErrors(cudaMemcpy((*Aredt)->row_ptr, (*Ared)->row_ptr, (new_m+1) * sizeof(int),   cudaMemcpyDeviceToDevice));
   checkCudaErrors(cudaMemcpy((*Aredt)->col_ind, (*Ared)->col_ind, nnz_new   * sizeof(int),   cudaMemcpyDeviceToDevice));
 
   int *d_A_to_At_ind;
   csr_transpose(*Aredt, &d_A_to_At_ind);
 
   // Update merge path buffer (CsrmvEx)
   update_mp_buffer(*Aredt);
 
   cuda_free((void**)&d_A_to_At_ind);
   cuda_free((void**)&d_predicate);
   cuda_free((void**)&d_compact_address);
   cuda_free((void**)&d_row_predicate);
   cuda_free((void**)&d_new_row_number);
 }
 
 void cuda_mat_get_m(const csr *mat,
                     int     *m) {
 
   (*m) = mat->m;
 }
 
 void cuda_mat_get_n(const csr *mat,
                     int     *n) {
 
   (*n) = mat->n;
 }
 
 void cuda_mat_get_nnz(const csr *mat,
                       int     *nnz) {
 
   (*nnz) = mat->nnz;
 }