#ifndef CONSTNATS_H
#define CONSTANTS_H

# define RHO (0.1)
# define SIGMA (1E-06)
# define MAX_ITER (4000)
# define EPS_ABS (1E-3)
# define EPS_REL (1E-3)
# define EPS_PRIM_INF (1E-4)
# define EPS_DUAL_INF (1E-4)
# define ALPHA (1.6)

# define LINSYS_SOLVER (LDL_SOLVER)

# define RHO_MIN (1e-06)
# define RHO_MAX (1e06)
# define RHO_EQ_OVER_RHO_INEQ (1e03)
# define RHO_TOL (1e-04) // tolerance for detecting if an inequality is set to equality
# define MIN_SCALING (1e-04) // minimum scaling value
# define MAX_SCALING (1e+04) // maximum scaling value
# define qpNAN ((float)0x7fc00000UL)  // not a number
# define qpINFTY ((float)1e30) // infinity
# define ADAPTIVE_RHO (1)
# define ADAPTIVE_RHO_INTERVAL (0)
# define ADAPTIVE_RHO_FRACTION (0.4) // fraction of setup time after which we update rho
# define ADAPTIVE_RHO_MULTIPLE_TERMINATION (4) // multiple of check_termination after which we update rho (if PROFILING disabled)
# define ADAPTIVE_RHO_FIXED (100) // number of iterations after which we update rho if termination_check  and PROFILING are disabled
# define ADAPTIVE_RHO_TOLERANCE (5) // tolerance for adopting new rho; minimum ratio between new rho and the current one

/******************
* Solver Status  *
******************/
# define qp_DUAL_INFEASIBLE_INACCURATE (4)
# define qp_PRIMAL_INFEASIBLE_INACCURATE (3)
# define qp_SOLVED_INACCURATE (2)
# define qp_SOLVED (1)
# define qp_MAX_ITER_REACHED (-2)
# define qp_PRIMAL_INFEASIBLE (-3)    /* primal infeasible  */
# define qp_DUAL_INFEASIBLE (-4)      /* dual infeasible */
# define qp_SIGINT (-5)               /* interrupted by user */
#  define qp_TIME_LIMIT_REACHED (-6)
# define qp_NON_CVX (-7)              /* problem non convex */
# define qp_UNSOLVED (-10)            /* Unsolved. Only setup function has been called */


/*************************
* Linear System Solvers *
*************************/
enum linsys_solver_type { LDL_SOLVER };
extern const char * LINSYS_SOLVER_NAME[];


/******************
* Solver Errors  *
******************/
enum qp_error_type {
    QP_DATA_VALIDATION_ERROR = 1,  /* Start errors from 1 */
    QP_SETTINGS_VALIDATION_ERROR,
    QP_LINSYS_SOLVER_LOAD_ERROR,
    QP_LINSYS_SOLVER_INIT_ERROR,
    QP_NONCVX_ERROR,
    QP_MEM_ALLOC_ERROR,
    QP_WORKSPACE_NOT_INIT_ERROR,
};
extern const char * QP_ERROR_MESSAGE[];

#endif CONSTANTS_H
