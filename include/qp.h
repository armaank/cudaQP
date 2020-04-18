/* main header files for qp */
#ifndef QP_H
# define QP_H

# include "ops.h"

/* todo: add timer */

/* problem scaling matrices stored as vectors */
typedef struct {
    float  c;    // cost function scaling
    float *D;    // primal variable scaling
    float *E;    // dual variable scaling
    float  cinv; // cost function rescaling
    float *Dinv; // primal variable rescaling
    float *Einv; // dual variable rescaling
} qpScaling;

/* solution structure */
typedef struct {
    float *x; // primal solution
    float *y; // Lagrange multiplier associated to \f$l <= Ax <= u\f$
} qpSolution;

/* solver information */
typedef struct {
    int iter; // number of iterations
    char status[32]; // status string, e.g. 'solved'
    int status_val; // status defined in constants.h

    int status_polish; // polish status: successful (1), unperformed (0), (-1) unsuccessful

    float obj_val;     // primal objective
    float pri_res;     // norm of primal residual
    float dua_res;     // norm of dual residual

    float setup_time;  // time taken for setup phase (seconds)
    float solve_time;  // time taken for solve phase (seconds)
    float update_time; // time taken for update phase (seconds)
    float polish_time; // time taken for polish phase (seconds)
    float run_time;    // total time  (seconds)

    int   rho_updates;  // number of rho updates
    float rho_estimate; // best rho estimate so far from residuals
} qpInfo;

/* main data struct */
typedef struct {
    int n; // number of variables n
    int m; // number of constraints m
    csc *P; // the upper triangular part of the quadratic cost matrix P in csc format (size n x n).
    csc *A; // linear constraints matrix A in csc format (size m x n)
    float *q; // dense array for linear part of cost function (size n)
    float *l; // dense array for lower bound (size m)
    float *u; // dense array for upper bound (size m)
} qpData;


/* hyperparams */
typedef struct {
    float rho; // ADMM step rho
    float sigma; // ADMM step sigma
    int scaling; // heuristic data scaling iterations; if 0, then disabled.

    int adaptive_rho; // boolean, is rho step size adaptive?
    int adaptive_rho_interval; // number of iterations between rho adaptations; if 0, then it is automatic
    float adaptive_rho_tolerance; // tolerance X for adapting rho. The new rho has to be X times larger or 1/X times smaller than the current one to trigger a new factorization.
    float adaptive_rho_fraction; // interval for adapting rho (fraction of the setup time)

    int max_iter; // maximum number of iterations
    float eps_abs; // absolute convergence tolerance
    float eps_rel; // relative convergence tolerance
    float eps_prim_inf; // primal infeasibility tolerance
    float eps_dual_inf; // dual infeasibility tolerance
    float alpha; // relaxation parameter

    float delta; // regularization parameter for polishing
    int polish; // boolean, polish ADMM solution
    int polish_refine_iter; // number of iterative refinement steps in polishing

    int verbose; // boolean, write out progress

    int scaled_termination; // boolean, use scaled termination criteria
    int check_termination; // integer, check termination interval; if 0, then termination checking is disabled
    int warm_start; // boolean, warm start

    float time_limit; // maximum number of seconds allowed to solve the problem; if 0, then disabled
} qpHyperparams



