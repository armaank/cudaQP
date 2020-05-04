
#ifndef CUDA_CSR_H
# define CUDA_CSR_H

#include "cuda_configure.h"
 #include "cuda_handler.h"
 #include "cuda_malloc.h"
 #include "cuda_wrapper.h"
 #include "helper_cuda.h"    /* --> checkCudaErrors */
 
 #include "csr_type.h"
//  #include "glob_opts.h" 

 #include <thrust/scan.h>
 #include <thrust/execution_policy.h>

#include "../../../include/qptypes.h" // for csc 

// #include "algebra_types.h"
// #include "csc_type.h"
// todo: figure out what to include here

void cuda_mat_init_P(const csc  *mat,
                     csr       **P,
                     float   **d_P_triu_val,
                     int     **d_P_triu_to_full_ind,
                     int     **d_P_diag_ind);
                     
void cuda_mat_init_A(const csc  *mat,
                     csr       **A,
                     csr       **At,
                     int     **d_A_to_At_ind);

void cuda_mat_update_P(const float  *Px,
                       const int    *Px_idx,
                       int           Px_n,
                       csr           **P,
                       float        *d_P_triu_val,
                       int          *d_P_triu_to_full_ind,
                       int          *d_P_diag_ind,
                       int           P_triu_nnz);

void cuda_mat_update_A(const float  *Ax,
                       const int    *Ax_idx,
                       int           Ax_n,
                       csr           **A,
                       csr           **At,
                       int          *d_A_to_At_ind);

void cuda_mat_free(csr *mat);

void cuda_submat_byrows(const csr    *A,
                        const int  *d_rows,
                        csr         **Ared,
                        csr         **Aredt);

void cuda_mat_get_m(const csr *mat,
                    int     *m);

void cuda_mat_get_n(const csr *mat,
                    int     *n);

void cuda_mat_get_nnz(const csr *mat,
                      int     *nnz);



#endif /* ifndef CUDA_CSR_H */