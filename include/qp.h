#ifndef QP_H
# define QP_H


#include "auxil.h"
# include "linsys.h"
// # include "polish.h"
# include "ctrl.h"
#include "scaling.h"
#include "qptypes.h" // might remove this? 
// # include "util.h" // Needed for qp_set_default_params functions


// #  include "cs.h"

/********************
* Main Solver API  *
********************/

/**
 * @name Main solver API
 * @{
 */

/**
 * Set default params from constants.h file
 * assumes params already allocated in memory
 * @param params params structure
 */
void qp_set_default_params(qpParams *params);



/**
 * Initialize qp solver allocating memory.
 *
 * All the inputs must be already allocated in memory before calling.
 *
 * It performs:
 * - data and params validation
 * - problem data scaling
 * - automatic parameters tuning (if enabled)
 * - setup linear system solver:
 *      - direct solver: KKT matrix factorization is performed here
 *      - indirect solver: KKT matrix preconditioning is performed here
 *
 * NB: This is the only function that allocates dynamic memory and is not used
 *during code generation
 *
 * @param  workp        Solver workspace pointer
 * @param  data         Problem data
 * @param  params     Solver params
 * @return              Exitflag for errors (0 if no errors)
 */
int qp_setup(qpWorkspace** workp, const qpData* data, const qpParams* params);


/**
 * Solve quadratic program
 *
 * The final solver information is stored in the a work->info  structure
 *
 * The solution is stored in the  a work->solution  structure
 *
 * If the problem is primal infeasible, the certificate is stored
 * in a work->delta_y
 *
 * If the problem is dual infeasible, the certificate is stored in a
 * work->delta_x
 *
 * @param  work Workspace allocated
 * @return      Exitflag for errors
 */
int qp_solve(qpWorkspace *work);



/**
 * Cleanup workspace by deallocating memory
 *
 * This function is not used in code generation
 * @param  work Workspace
 * @return      Exitflag for errors
 */
int qp_cleanup(qpWorkspace *work);


/** @} */


/********************************************
* Sublevel API                             *
*                                          *
* Edit data without performing setup again *
********************************************/

/**
 * @name Sublevel API
 * @{
 */

/**
 * Update linear cost in the problem
 * @param  work  Workspace
 * @param  q_new New linear cost
 * @return       Exitflag for errors and warnings
 */
int qp_update_lin_cost(qpWorkspace *work,
                       const float *q_new);


/**
 * Update lower and upper bounds in the problem constraints
 * @param  work   Workspace
 * @param  l_new New lower bound
 * @param  u_new New upper bound
 * @return        Exitflag: 1 if new lower bound is not <= than new upper bound
 */
int qp_update_bounds(qpWorkspace *work,
                     const float *l_new,
                     const float *u_new);


/**
 * Update lower bound in the problem constraints
 * @param  work   Workspace
 * @param  l_new New lower bound
 * @return        Exitflag: 1 if new lower bound is not <= than upper bound
 */
int qp_update_lower_bound(qpWorkspace *work,
                          const float *l_new);


/**
 * Update upper bound in the problem constraints
 * @param  work   Workspace
 * @param  u_new New upper bound
 * @return        Exitflag: 1 if new upper bound is not >= than lower bound
 */
int qp_update_upper_bound(qpWorkspace *work,
                          const float *u_new);


/**
 * Warm start primal and dual variables
 * @param  work Workspace structure
 * @param  x    Primal variable
 * @param  y    Dual variable
 * @return      Exitflag
 */
int qp_warm_start(qpWorkspace *work,
                  const float *x,
                  const float *y);


/**
 * Warm start primal variable
 * @param  work Workspace structure
 * @param  x    Primal variable
 * @return      Exitflag
 */
int qp_warm_start_x(qpWorkspace *work,
                    const float *x);


/**
 * Warm start dual variable
 * @param  work Workspace structure
 * @param  y    Dual variable
 * @return      Exitflag
 */
int qp_warm_start_y(qpWorkspace *work,
                    const float *y);


/**
 * Update elements of matrix P (upper triangular)
 * without changing sparsity structure.
 *
 *
 *  If Px_new_idx is qp_NULL, Px_new is assumed to be as long as P->x
 *  and the whole P->x is replaced.
 *
 * @param  work       Workspace structure
 * @param  Px_new     Vector of new elements in P->x (upper triangular)
 * @param  Px_new_idx Index mapping new elements to positions in P->x
 * @param  P_new_n    Number of new elements to be changed
 * @return            output flag:  0: OK
 *                                  1: P_new_n > nnzP
 *                                 <0: error in the update
 */
int qp_update_P(qpWorkspace *work,
                const float *Px_new,
                const int   *Px_new_idx,
                int          P_new_n);


/**
 * Update elements of matrix A without changing sparsity structure.
 *
 *
 *  If Ax_new_idx is qp_NULL, Ax_new is assumed to be as long as A->x
 *  and the whole A->x is replaced.
 *
 * @param  work       Workspace structure
 * @param  Ax_new     Vector of new elements in A->x
 * @param  Ax_new_idx Index mapping new elements to positions in A->x
 * @param  A_new_n    Number of new elements to be changed
 * @return            output flag:  0: OK
 *                                  1: A_new_n > nnzA
 *                                 <0: error in the update
 */
int qp_update_A(qpWorkspace *work,
                const float *Ax_new,
                const int   *Ax_new_idx,
                int          A_new_n);


/**
 * Update elements of matrix P (upper triangular) and elements of matrix A
 * without changing sparsity structure.
 *
 *
 *  If Px_new_idx is qp_NULL, Px_new is assumed to be as long as P->x
 *  and the whole P->x is replaced.
 *
 *  If Ax_new_idx is qp_NULL, Ax_new is assumed to be as long as A->x
 *  and the whole A->x is replaced.
 *
 * @param  work       Workspace structure
 * @param  Px_new     Vector of new elements in P->x (upper triangular)
 * @param  Px_new_idx Index mapping new elements to positions in P->x
 * @param  P_new_n    Number of new elements to be changed
 * @param  Ax_new     Vector of new elements in A->x
 * @param  Ax_new_idx Index mapping new elements to positions in A->x
 * @param  A_new_n    Number of new elements to be changed
 * @return            output flag:  0: OK
 *                                  1: P_new_n > nnzP
 *                                  2: A_new_n > nnzA
 *                                 <0: error in the update
 */
int qp_update_P_A(qpWorkspace *work,
                  const float *Px_new,
                  const int   *Px_new_idx,
                  int          P_new_n,
                  const float *Ax_new,
                  const int   *Ax_new_idx,
                  int          A_new_n);

/**
 * Update rho. Limit it between RHO_MIN and RHO_MAX.
 * @param  work         Workspace
 * @param  rho_new      New rho setting
 * @return              Exitflag
 */
int qp_update_rho(qpWorkspace *work,
                  float        rho_new);


/** @} */


/**
 * @name Update params
 * @{
 */

/**
 * Update max_iter setting
 * @param  work         Workspace
 * @param  max_iter_new New max_iter setting
 * @return              Exitflag
 */
int qp_update_max_iter(qpWorkspace *work,
                       int          max_iter_new);


/**
 * Update absolute tolernace value
 * @param  work        Workspace
 * @param  eps_abs_new New absolute tolerance value
 * @return             Exitflag
 */
int qp_update_eps_abs(qpWorkspace *work,
                      float        eps_abs_new);


/**
 * Update relative tolernace value
 * @param  work        Workspace
 * @param  eps_rel_new New relative tolerance value
 * @return             Exitflag
 */
int qp_update_eps_rel(qpWorkspace *work,
                      float        eps_rel_new);


/**
 * Update primal infeasibility tolerance
 * @param  work          Workspace
 * @param  eps_prim_inf_new  New primal infeasibility tolerance
 * @return               Exitflag
 */
int qp_update_eps_prim_inf(qpWorkspace *work,
                           float        eps_prim_inf_new);


/**
 * Update dual infeasibility tolerance
 * @param  work          Workspace
 * @param  eps_dual_inf_new  New dual infeasibility tolerance
 * @return               Exitflag
 */
int qp_update_eps_dual_inf(qpWorkspace *work,
                           float        eps_dual_inf_new);


/**
 * Update relaxation parameter alpha
 * @param  work  Workspace
 * @param  alpha_new New relaxation parameter value
 * @return       Exitflag
 */
int qp_update_alpha(qpWorkspace *work,
                    float        alpha_new);


/**
 * Update warm_start setting
 * @param  work           Workspace
 * @param  warm_start_new New warm_start setting
 * @return                Exitflag
 */
int qp_update_warm_start(qpWorkspace *work,
                         int          warm_start_new);


/**
 * Update scaled_termination setting
 * @param  work                 Workspace
 * @param  scaled_termination_new  New scaled_termination setting
 * @return                      Exitflag
 */
int qp_update_scaled_termination(qpWorkspace *work,
                                 int          scaled_termination_new);

/**
 * Update check_termination setting
 * @param  work                   Workspace
 * @param  check_termination_new  New check_termination setting
 * @return                        Exitflag
 */
int qp_update_check_termination(qpWorkspace *work,
                                int          check_termination_new);


/**
 * Update regularization parameter in polish
 * @param  work      Workspace
 * @param  delta_new New regularization parameter
 * @return           Exitflag
 */
int qp_update_delta(qpWorkspace *work,
                    float        delta_new);


/**
 * Update polish setting
 * @param  work          Workspace
 * @param  polish_new New polish setting
 * @return               Exitflag
 */
int qp_update_polish(qpWorkspace *work,
                     int          polish_new);


/**
 * Update number of iterative refinement steps in polish
 * @param  work                Workspace
 * @param  polish_refine_iter_new New iterative reginement steps
 * @return                     Exitflag
 */
int qp_update_polish_refine_iter(qpWorkspace *work,
                                 int          polish_refine_iter_new);


/**
 * Update verbose setting
 * @param  work        Workspace
 * @param  verbose_new New verbose setting
 * @return             Exitflag
 */
int qp_update_verbose(qpWorkspace *work,
                      int          verbose_new);


/**
 * Update time_limit setting
 * @param  work            Workspace
 * @param  time_limit_new  New time_limit setting
 * @return                 Exitflag
 */
int qp_update_time_limit(qpWorkspace *work,
                         float        time_limit_new);


#endif //QP_H


