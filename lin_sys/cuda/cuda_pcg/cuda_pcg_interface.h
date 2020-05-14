#ifndef CUDA_PCG_INTERFACE_H
#define CUDA_PCG_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "osqp.h"
#include "types.h"
#include "algebra_types.h"

#include "cuda_pcg_constants.h"

/* cuda pcg solver structure */
typedef struct cudapcg_solver_ {

    enum linsys_solver_type type;

    c_int (*solve)(struct cudapcg_solver_ *self, OSQPVectorf *b, c_int admm_iter);

    void (*warm_start)(struct cudapcg_solver_ *self, const OSQPVectorf *x);

    void (*free)(struct cudapcg_solver_ *self);

    c_int (*update_matrices)(struct cudapcg_solver_ *self, const OSQPMatrix *P, const OSQPMatrix *A);

    c_int (*update_rho_vec)(struct cudapcg_solver_ *self, const OSQPVectorf *rho_vec, c_float rho_sc);

    /* threads count */
    c_int nthreads;

    /* dimensions */
    c_int n;
    c_int m;

    /* states */
    c_int polish;
    c_int zero_pcg_iters;

    /* settings */
    enum pcg_eps_strategy eps_strategy;
    c_int                 norm;
    c_int                 precondition;
    c_int                 warm_start_pcg;
    c_int                 max_iter;

    /* scs tolerance strategy parameters */
    c_float start_tol;
    c_float dec_rate;

    /* residual tolerance strategy parameters */
    c_int    reduction_threshold;
    c_float  reduction_factor;
    c_float  eps_prev;
    c_float *scaled_pri_res;
    c_float *scaled_dua_res;

    /* pointer to problem data and ADMM settings */
    csr     *A;
    csr     *At;
    csr     *P;
    c_int   *d_P_diag_ind;
    c_float *d_rho_vec;
    c_float *h_sigma;
    c_float *h_rho;

    /* pcg iterates */
    c_float *d_x;
    c_float *d_p;
    c_float *d_Kp;
    c_float *d_y;
    c_float *d_r;
    c_float *d_rhs;
    c_float *d_z;

    /* pointer to host memory */
    c_float *h_r_norm;

    /* pcg scalar values (in target memory) */
    c_float *d_r_norm;
    c_float *rTy;
    c_float *rTy_prev;
    c_float *alpha;
    c_float *beta;
    c_float *pKp;
    c_float *D_MINUS_ONE;
    c_float *d_sigma;

    /* pcg preconditioner */
    c_float *d_P_diag_val;
    c_float *d_AtA_diag_val;
    c_float *d_AtRA_diag_val;
    c_float *d_diag_precond;
    c_float *d_diag_precond_inv;

    /* function pointer to handle different vector norms */
    void (*vector_norm)(const c_float *d_x, c_int n, c_float *res);

} cudapcg_solver;

c_int init_linsys_solver_cudapcg(cudapcg_solver **sp,
                                 const OSQPMatrix *P,
                                 const OSQPMatrix *A,
                                 const OSQPVectorf *rho_vec,
                                 OSQPSettings *settings,
                                 c_float *scaled_pri_res,
                                 c_float *scaled_dua_res,
                                 c_int polish);


c_int solve_linsys_cudapcg(cudapcg_solver *s, OSQPVectorf *b, c_int admm_iter);

void warm_start_linsys_solver_cudapcg(cudapcg_solver *s, const OSQPVectorf *x);

void free_linsys_solver_cudapcg(cudapcg_solver *s);

c_int update_linsys_solver_matrices_cudapcg(cudapcg_solver *s, const OSQPMatrix *P, const OSQPMatrix *A);

c_int update_linsys_solver_rho_vec_cudapcg(cudapcg_solver *s, const OSQPVectorf *rho_vec, c_float rho_sc);


#ifdef __cplusplus
}
#endif

#endif /* ifndef OSQP_API_TYPES_H */
