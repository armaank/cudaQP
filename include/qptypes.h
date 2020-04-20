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
  // enum linsys_solver_type linsys_solver; ///< linear system solver to use

  float delta;                         ///< regularization parameter for polishing
  int   polish;                        ///< boolean, polish ADMM solution
  int   polish_refine_iter;            ///< number of iterative refinement steps in polishing

  int verbose;                         ///< boolean, write out progress

  int scaled_termination;              ///< boolean, use scaled termination criteria
  int check_termination;               ///< integer, check termination interval; if 0, then termination checking is disabled
  int warm_start;                      ///< boolean, warm start

  float time_limit;                    ///< maximum number of seconds allowed to solve the problem; if 0, then disabled
} qpParams;
