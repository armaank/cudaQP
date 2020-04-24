#include "constants.h"

typedef struct {
  int    nzmax; ///< maximum number of entries
  int    m;     ///< number of rows
  int    n;     ///< number of columns
  int   *p;     ///< column pointers (size n+1); col indices (size nzmax) start from 0 when using triplet format (direct KKT matrix formation)
  int   *i;     ///< row indices, size nzmax starting from 0
  float *x;     ///< numerical values, size nzmax
  int    nz;    ///< number of entries in triplet matrix, -1 for csc
} csc;


typedef struct linsys_solver LinSysSolver;
typedef struct QP_TIMER qpTimer;

typedef struct {
  float  c;    ///< cost function scaling
  float *D;    ///< primal variable scaling
  float *E;    ///< dual variable scaling
  float  cinv; ///< cost function rescaling
  float *Dinv; ///< primal variable rescaling
  float *Einv; ///< dual variable rescaling
} qpScaling;

typedef struct {
  float *x; ///< primal solution
  float *y; ///< Lagrange multiplier associated to \f$l <= Ax <= u\f$
} qpSolution;

typedef struct {
  int iter;          ///< number of iterations taken
  char  status[32];    ///< status string, e.g. 'solved'
  int status_val;    ///< status as int, defined in constants.h

  int status_polish; ///< polish status: successful (1), unperformed (0), (-1) unsuccessful

  float obj_val;     ///< primal objective
  float pri_res;     ///< norm of primal residual
  float dua_res;     ///< norm of dual residual

  float setup_time;  ///< time taken for setup phase (seconds)
  float solve_time;  ///< time taken for solve phase (seconds)
  float update_time; ///< time taken for update phase (seconds)
  float polish_time; ///< time taken for polish phase (seconds)
  float run_time;    ///< total time  (seconds)

  int   rho_updates;  ///< number of rho updates
  float rho_estimate; ///< best rho estimate so far from residuals
} qpInfo;

typedef struct {
  int    n; ///< number of variables n
  int    m; ///< number of constraints m
  csc     *P; ///< the upper triangular part of the quadratic cost matrix P in csc format (size n x n).
  csc     *A; ///< linear constraints matrix A in csc format (size m x n)
  float *q; ///< dense array for linear part of cost function (size n)
  float *l; ///< dense array for lower bound (size m)
  float *u; ///< dense array for upper bound (size m)
} qpData;


typedef struct {
  csc *Ared;          ///< active rows of A
  ///<    Ared = vstack[Alow, Aupp]
  int    n_low;     ///< number of lower-active rows
  int    n_upp;     ///< number of upper-active rows
  int   *A_to_Alow; ///< Maps indices in A to indices in Alow
  int   *A_to_Aupp; ///< Maps indices in A to indices in Aupp
  int   *Alow_to_A; ///< Maps indices in Alow to indices in A
  int   *Aupp_to_A; ///< Maps indices in Aupp to indices in A
  float *x;         ///< optimal x-solution obtained by polish
  float *z;         ///< optimal z-solution obtained by polish
  float *y;         ///< optimal y-solution obtained by polish
  float  obj_val;   ///< objective value at polished solution
  float  pri_res;   ///< primal residual at polished solution
  float  dua_res;   ///< dual residual at polished solution
} qpPolish;



typedef struct {
  float rho;                    ///< ADMM step rho
  float sigma;                  ///< ADMM step sigma
  int   scaling;                ///< heuristic data scaling iterations; if 0, then disabled.

  int   adaptive_rho;           ///< boolean, is rho step size adaptive?
  int   adaptive_rho_interval;  ///< number of iterations between rho adaptations; if 0, then it is automatic
  float adaptive_rho_tolerance; ///< tolerance X for adapting rho. The new rho has to be X times larger or 1/X times smaller than the current one to trigger a new factorization.
  float adaptive_rho_fraction;  ///< interval for adapting rho (fraction of the setup time)


  int                   max_iter;      ///< maximum number of iterations
  float                 eps_abs;       ///< absolute convergence tolerance
  float                 eps_rel;       ///< relative convergence tolerance
  float                 eps_prim_inf;  ///< primal infeasibility tolerance
  float                 eps_dual_inf;  ///< dual infeasibility tolerance
  float                 alpha;         ///< relaxation parameter
  
  enum linsys_solver_type linsys_solver; ///< linear system solver to use

  float delta;                         ///< regularization parameter for polishing
  int   polish;                        ///< boolean, polish ADMM solution
  int   polish_refine_iter;            ///< number of iterative refinement steps in polishing

  int verbose;                         ///< boolean, write out progress

  int scaled_termination;              ///< boolean, use scaled termination criteria
  int check_termination;               ///< integer, check termination interval; if 0, then termination checking is disabled
  int warm_start;                      ///< boolean, warm start

  float time_limit;                    ///< maximum number of seconds allowed to solve the problem; if 0, then disabled
} qpParams;


typedef struct {
  /// Problem data to work on (possibly scaled)
  qpData *data;

  /// Linear System solver structure
  LinSysSolver *linsys_solver;

  /// Polish structure
  qpPolish *pol;

  /**
   * @name Vector used to store a vectorized rho parameter
   * @{
   */
  float *rho_vec;     ///< vector of rho values
  float *rho_inv_vec; ///< vector of inv rho values

  /** @} */

  int *constr_type; ///< Type of constraints: loose (-1), equality (1), inequality (0)

  /**
   * @name Iterates
   * @{
   */
  float *x;        ///< Iterate x
  float *y;        ///< Iterate y
  float *z;        ///< Iterate z
  float *xz_tilde; ///< Iterate xz_tilde

  float *x_prev;   ///< Previous x

  /**< NB: Used also as workspace vector for dual residual */
  float *z_prev;   ///< Previous z

  /**< NB: Used also as workspace vector for primal residual */

  /**
   * @name Primal and dual residuals workspace variables
   *
   * Needed for residuals computation, tolerances computation,
   * approximate tolerances computation and adapting rho
   * @{
   */
  float *Ax;  ///< scaled A * x
  float *Px;  ///< scaled P * x
  float *Aty; ///< scaled A' * y

  /** @} */

  /**
   * @name Primal infeasibility variables
   * @{
   */
  float *delta_y;   ///< difference between consecutive dual iterates
  float *Atdelta_y; ///< A' * delta_y

  /** @} */

  /**
   * @name Dual infeasibility variables
   * @{
   */
  float *delta_x;  ///< difference between consecutive primal iterates
  float *Pdelta_x; ///< P * delta_x
  float *Adelta_x; ///< A * delta_x

  /** @} */

  /**
   * @name Temporary vectors used in scaling
   * @{
   */

  float *D_temp;   ///< temporary primal variable scaling vectors
  float *D_temp_A; ///< temporary primal variable scaling vectors storing norms of A columns
  float *E_temp;   ///< temporary constraints scaling vectors storing norms of A' columns


  /** @} */

  qpParams *params; ///< problem settings
  qpScaling  *scaling;  ///< scaling vectors
  qpSolution *solution; ///< problem solution
  qpInfo     *info;     ///< solver information

  qpTimer *timer;       ///< timer object

  /// flag indicating whether the solve function has been run before
  int first_run;

  /// flag indicating whether the update_time should be cleared
  int clear_update_time;

  /// flag indicating that osqp_update_rho is called from osqp_solve function
  int rho_update_from_solve;

  int summary_printed; ///< Has last summary been printed? (true/false)

} qpWorkspace;

struct linsys_solver {
  enum linsys_solver_type type;                 ///< linear system solver type functions
  int (*solve)(LinSysSolver *self,
                 float      *b);              ///< solve linear system

  void (*free)(LinSysSolver *self);             ///< free linear system solver (only in desktop version)

  int (*update_matrices)(LinSysSolver *s,
                           const csc *P,            ///< update matrices P
                           const csc *A);           //   and A in the solver

  int (*update_rho_vec)(LinSysSolver  *s,
                          const float *rho_vec);  ///< Update rho_vec

  int nthreads; ///< number of threads active
};
