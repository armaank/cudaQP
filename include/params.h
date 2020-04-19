/* params */

#define RHO (0.1)
#define SIGMA (1E-06)
#define MAX_ITER (4000)
#define EPS_ABS (1E-3)
#define EPS_REL (1E-3)
#define EPS_PRIM_INF (1E-4)
#define EPS_DUAL_INF (1E-4)
#define ALPHA (1.6)

#define RHO_MIN (1e-06)
#define RHO_MAX (1e06)
#define RHO_EQ_OVER_RHO_INEQ (1e03)
#define RHO_TOL (1e-04) ///< tolerance for detecting if an inequality is set to equality

#define DELTA (1E-6)
#define POLISH (0)
#define POLISH_REFINE_ITER (3)

#define SCALED_TERMINATION (0)
#define CHECK_TERMINATION (25)
#define WARM_START (1)
#define SCALING (10)

#define INFTY ((float)1e30)   
#define qpNAN ((float)0x7fc00000UL)

#define MIN_SCALING (1e-04) ///< minimum scaling value
#define MAX_SCALING (1e+04) ///< maximum scaling value

#define ADAPTIVE_RHO (1)
#define ADAPTIVE_RHO_INTERVAL (0)
#define ADAPTIVE_RHO_FRACTION (0.4)           ///< fraction of setup time after which we update rho
#define ADAPTIVE_RHO_MULTIPLE_TERMINATION (4) ///< multiple of check_termination after which we update rho (if PROFILING disabled)
#define ADAPTIVE_RHO_FIXED (100)              ///< number of iterations after which we update rho if termination_check  and PROFILING are disabled
#define ADAPTIVE_RHO_TOLERANCE (5)            ///< tolerance for adopting new rho; minimum ratio between new rho and the current one

# define QP_DUAL_INFEASIBLE_INACCURATE (4)
# define QP_PRIMAL_INFEASIBLE_INACCURATE (3)
# define QP_SOLVED_INACCURATE (2)
# define QP_SOLVED (1)
# define QP_MAX_ITER_REACHED (-2)
# define QP_PRIMAL_INFEASIBLE (-3)    /* primal infeasible  */
# define QP_DUAL_INFEASIBLE (-4)      /* dual infeasible */
# define QP_SIGINT (-5)               /* interrupted by user */
# ifdef PROFILING
#  define QP_TIME_LIMIT_REACHED (-6)
# endif // ifdef PROFILING
# define QP_NON_CVX (-7)              /* problem non convex */
# define QP_UNSOLVED (-10)            /* Unsolved. Only setup function has been called */
