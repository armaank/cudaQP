# define RHO (0.1)
# define SIGMA (1E-06)
# define MAX_ITER (4000)
# define EPS_ABS (1E-3)
# define EPS_REL (1E-3)
# define EPS_PRIM_INF (1E-4)
# define EPS_DUAL_INF (1E-4)
# define ALPHA (1.6)
# define RHO_MIN (1e-06)
# define RHO_MAX (1e06)
# define RHO_EQ_OVER_RHO_INEQ (1e03)
# define RHO_TOL (1e-04) // tolerance for detecting if an inequality is set to equality
# define MIN_SCALING (1e-04) // minimum scaling value
# define MAX_SCALING (1e+04) // maximum scaling value
# define QPNAN ((c_float)0x7fc00000UL)  // not a number
# define QPINFTY ((c_float)1e30) // infinity
# define ADAPTIVE_RHO (1)
# define ADAPTIVE_RHO_INTERVAL (0)
# define ADAPTIVE_RHO_FRACTION (0.4) // fraction of setup time after which we update rho
# define ADAPTIVE_RHO_MULTIPLE_TERMINATION (4) // multiple of check_termination after which we update rho (if PROFILING disabled)
# define ADAPTIVE_RHO_FIXED (100) // number of iterations after which we update rho if termination_check  and PROFILING are disabled
# define ADAPTIVE_RHO_TOLERANCE (5) // tolerance for adopting new rho; minimum ratio between new rho and the current one
