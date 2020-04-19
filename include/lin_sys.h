#ifndef LIN_SYS_H
# define LIN_SYS_H

/* KKT linear system definition and solution */
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
/**
 * Load linear system solver shared library
 * @param	linsys_solver  Linear system solver
 * @return Zero on success, nonzero on failure.
 */
int load_linsys_solver(enum linsys_solver_type linsys_solver);


/**
 * Unload linear system solver shared library
 * @param	linsys_solver  Linear system solver
 * @return Zero on success, nonzero on failure.
 */
int unload_linsys_solver(enum linsys_solver_type linsys_solver);


// NB: Only the upper triangular part of P is stuffed!

/**
 * Initialize linear system solver structure
 * @param   s             Pointer to linear system solver structure
 * @param   P             Cost function matrix
 * @param	A             Constraint matrix
 * @param	sigma         Algorithm parameter
 * @param	rho_vec       Algorithm parameter
 * @param	linsys_solver Linear system solver
 * @param	polish        0/1 depending whether we are allocating for
 *polishing or not
 * @return                Exitflag for error (0 if no errors)
 */
int init_linsys_solver(LinSysSolver          **s,
                         const csc              *P,
                         const csc              *A,
                         float                 sigma,
                         const float          *rho_vec,
                         enum linsys_solver_type linsys_solver,
                         int                   polish);


#endif // ifndef LIN_SYS_H