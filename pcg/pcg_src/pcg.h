#ifndef PCG_H
# define PCG_H
 
#include "../pcg_interface.h" // is this needed, could be ancillary
#include "../include/csr_type.h"
#include "../include/cuda_handler.h"
#include "../include/cuda_malloc.h"
#include "../include/cuda_lin_alg.h"
#include "../include/cuda_wrapper.h"
#include "../include/helper_cuda.h"    /* --> checkCudaErrors */


/**
 *  Preconditioned Conjugate Gradient (PCG) algorithm.
 *  Computes an approximate solution to the linear system
 *
 *       K * x = rhs
 *
 *  The solution is stored in s->d_x.
 *  The function returns the number of PCG iterations evaluated.
 */
int cuda_pcg_alg(cudapcg_solver *s,
                   float         eps,
                   int           max_iter);

void cuda_pcg_update_precond(cudapcg_solver *s,
                             int           P_updated,
                             int           A_updated,
                             int           R_updated);



#endif /* ifndef CUDA_PCG_H */
