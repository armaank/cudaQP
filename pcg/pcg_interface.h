#ifndef PCG_INTERFACE_H
#define PCG_INTERFACE_H

// edit includes
// #include "osqp.h" // remove
// #include "types.h"  // remove       /* OSQPMatrix and OSQPVector[fi] types */
// #include "algebra_types.h"        /* csr type */

// #include "cuda_pcg_constants.h"   /* enum linsys_solver_type */
#include "./pcg_src/pcg_params.h"
#include "../include/qp.h"
#include "../include/qptypes.h" // for csc
#include "./pcg_src/pcg.h"

#include "./cuda/cuda_lin_alg.h"
#include "./cuda/cuda_malloc.h"

/**
 * CUDA PCG solver structure
 */
typedef struct cudapcg_solver_ {

  enum linsys_solver_type type;

  /**
   * @name Functions
   * @{
   */
  int (*solve)(struct cudapcg_solver_ *self,
                //  OSQPVectorf            *b,
                 const float *b
                 int                   admm_iter);

  void (*warm_start)(struct cudapcg_solver_  *self,
                    //  const OSQPVectorf       *x)
                        const float* x);

  void (*free)(struct cudapcg_solver_ *self);

  int (*update_matrices)(struct cudapcg_solver_ *self,
                        //    const OSQPMatrix       *P,
                        //    const OSQPMatrix       *A);
                              const csc *P,
                              const csc *A);

  int (*update_rho_vec)(struct cudapcg_solver_  *self,
                        //   const OSQPVectorf       *rho_vec,
                        const float *rho_vec,
                          float                  rho_sc);

  /* threads count */
  int nthreads;

  /* Dimensions */
  int n;                  ///<  dimension of the linear system
  int m;                  ///<  number of rows in A

  /* States */
  int polish;
  int zero_pcg_iters;     ///<  state that counts zero PCG iterations

  /* Settings */
  enum pcg_eps_strategy eps_strategy;
  int                 norm;
  int                 precondition;
  int                 warm_start_pcg;
  int                 max_iter;

  /* SCS tolerance strategy parameters */
  float start_tol;
  float dec_rate;
  
  /* Residual tolerance strategy parameters */
  int    reduction_threshold;
  float  reduction_factor;
  float  eps_prev;
  float *scaled_pri_res;
  float *scaled_dua_res;

  /* Pointers to problem data and ADMM settings */
  csr     *A;
  csr     *At;
  csr     *P;
  int   *d_P_diag_ind;
  float *d_rho_vec;
  float *h_sigma;
  float *h_rho;

  /* PCG iterates */
  float *d_x;             ///<  current iterate solution
  float *d_p;             ///<  current search direction
  float *d_Kp;            ///<  holds K*p
  float *d_y;             ///<  solution of the preconditioner r = M*y
  float *d_r;             ///<  residual r = K*x - b
  float *d_rhs;           ///<  right-hand side of Kx = b
  float *d_z;             ///<  holds z = A*x for computing A'*z = A'*(A*x);

  /* Pointer to page-locked host memory */
  float *h_r_norm;

  /* PCG scalar values (in device memory) */
  float *d_r_norm;
  float *rTy;
  float *rTy_prev;
  float *alpha;
  float *beta;
  float *pKp;
  float *D_MINUS_ONE;     ///<  constant -1.0 in device memory
  float *d_sigma;

  /* PCG preconditioner */
  float *d_P_diag_val;
  float *d_AtA_diag_val;
  float *d_AtRA_diag_val;
  float *d_diag_precond;
  float *d_diag_precond_inv;

  /* Function pointer to handle different vector norms */
  void (*vector_norm)(const float *d_x,
                      int          n,
                      float       *res);

} cudapcg_solver;



int init_linsys_solver_cudapcg(cudapcg_solver    **sp,
                                //  const OSQPMatrix   *P,
                                //  const OSQPMatrix   *A,
                                 const csc *P, 
                                 const csc *A,
                                //  const OSQPVectorf  *rho_vec,
                                 const float *rho_vec, 
                                 qpParams       *params,
                                 float            *scaled_pri_res,
                                 float            *scaled_dua_res,
                                 int               polish);


int solve_linsys_cudapcg(cudapcg_solver *s,
                           float *b,
                        //    OSQPVectorf    *b,
                           int           admm_iter);

void warm_start_linsys_solver_cudapcg(cudapcg_solver    *s,
                                    //   const OSQPVectorf *x
                                          const float* x);

void free_linsys_solver_cudapcg(cudapcg_solver *s);

int update_linsys_solver_matrices_cudapcg(cudapcg_solver   *s,
                                            // const OSQPMatrix *P,
                                            // const OSQPMatrix *A);
                                            const csc *P,
                                            const csc *A);

int update_linsys_solver_rho_vec_cudapcg(cudapcg_solver    *s,
                                        //    const OSQPVectorf *rho_vec,
                                           const float *rho_vec
                                           float            rho_sc);



#endif
