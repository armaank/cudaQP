/* main header files for qp */
#ifndef QP_H
# define QP_H

#include "ops.h"
#include "timer.h"
#include "linalg.h"


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
} qpHyperparams;

/* solution polishing */
typedef struct {
    csc *Ared; // active rows of A
    int n_low; // number of lower-active rows
    int n_upp; // number of upper-active rows
    int *A_to_Alow; // Maps indices in A to indices in Alow
    int *A_to_Aupp; // Maps indices in A to indices in Aupp
    int *Alow_to_A; // Maps indices in Alow to indices in A
    int *Aupp_to_A; // Maps indices in Aupp to indices in A
    float *x; // optimal x-solution obtained by polish
    float *z; // optimal z-solution obtained by polish
    float *y; // optimal y-solution obtained by polish
    float obj_val; // objective value at polished solution
    float pri_res; // primal residual at polished solution
    float dua_res; // dual residual at polished solution
} qpSolnPolish;

/* problem instance */
typedef struct {

    qpSolnPolish *pol; // solution polishing

    float *rho_vec; // vector of rho values
    float *rho_inv_vec; // vector of inv rho values
    float constr_type; // type of convex constraint: loose (-1), equality (1), inequality (0)
    float *x;        // Iterate x
    float *y;        // Iterate y
    float *z;        // Iterate z
    float *xz_tilde; // Iterate xz_tilde

    float *x_prev;   // Previous x

    float *z_prev;   // Previous z

    float *Ax;  // scaled A * x
    float *Px;  // scaled P * x
    float *Aty; // scaled A' * y

    float *delta_y;   // difference between consecutive dual iterates
    float *Atdelta_y; // A' * delta_y

    float *delta_x;  // difference between consecutive primal iterates
    float *Pdelta_x; // P * delta_x
    float *Adelta_x; // A * delta_x

    float *D_temp;   // temporary primal variable scaling vectors
    float *D_temp_A; // temporary primal variable scaling vectors storing norms of A columns
    float *E_temp;   // temporary constraints scaling vectors storing norms of A' columns

    qpHyperparams *hparams; // problem settings
    qpScaling *scaling;  // scaling vectors
    qpSolution *solution; // problem solution
    qpInfo *info;     // solver information

    qp_timer *timer;       // timer object

    // flag indicating whether the solve function has been run before
    int first_run;

    // flag indicating whether the update_time should be cleared
    int clear_update_time;

    // flag indicating that osqp_update_rho is called from osqp_solve function
    int rho_update_from_solve;

} qpInstance;

float compute_rho_estimate(qpInstance *problem);

/* Adapt rho value based on current unscaled primal/dual residuals */
int adapt_rho(qpInstance *problem);

/* Set values of rho vector based on constraint types */
void set_rho_vec(qpInstance *problem);

/* Update values of rho vector based on updated constraints */
int update_rho_vec(qpInstance *problem);

/* Swap float vector pointers */
void swap_vectors(float **a, float **b);

/* Cold start workspace variables xz and y */
void cold_start(qpInstance *problem);

/* Update x_tilde and z_tilde variable (first ADMM step) */
void update_xz_tilde(qpInstance *problem);

/*  Update x (second ADMM step) */
void update_x(qpInstance *problem);

/* Update z (third ADMM step) */
void update_z(qpInstance *problem);

/* Update y variable (fourth ADMM step) */
void update_y(qpInstance *problem);

/* Compute objective function from data at value x */
float compute_obj_val(qpInstance *problem, float *x);

/* Check whether QP has solution */
int has_solution(qpInfo *info);

/* Store the QP solution */
void store_solution(qpInstance *problem);

/* Update solver information */
void update_info(qpInstance *problem, int iter, int compute_objective, int polish);

/* Reset solver information (after problem updates) */
void reset_info(qpInfo *info);

/* update solver status (value and string) */
void update_status(qpInfo *info,int status_val);

/* check if termination conditions are satisfied */
int check_termination(qpInstance *problem, int approximate);

/* Validate problem data */
int validate_data(const qpData *data);

/* Validate problem settings */
int validate_settings(const qpHyperparams *settings);

/* Solution polish: Solve equality constrained QP with assumed active */
int polish(qpInstance *problem)

/* Define Projections onto set C involved in the ADMM algorithm */
void project(qpInstance *problem, float *z);

/*  Ensure z satisfies box constraints and y is is normal cone of z */
void project_normalcone(qpInstance *problem, float *z, float *y);

int scale_data(qpInstance *problem);

int unscale_data(qpInstance *problem);

int unscale_solution(qpInstance *problem);