#ifndef CUDA_PCG_H
# define CUDA_PCG_H
 
#include "cuda_pcg_interface.h"

#ifdef __cplusplus
extern "C" {
#endif

/* main pcg algo */
c_int cuda_pcg(cudapcg_solver *s, c_float eps,c_int max_iter);

/* update preconditioning */
void cuda_pcg_update_precond(cudapcg_solver *s, c_int P_updated, c_int A_updated,c_int R_updated);

#ifdef __cplusplus
}
#endif

#endif /* ifndef CUDA_PCG_H */