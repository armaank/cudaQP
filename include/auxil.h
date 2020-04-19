#ifndef AUXIL_H
# define AUXIL_H
# include "types.h"


/***********************************************************
* Auxiliary functions needed to compute ADMM iterations * *
***********************************************************/

/**
 * Compute rho estimate from residuals
 * @param work Workspace
 * @return     rho estimate
 */
float compute_rho_estimate(OSQPWorkspace *work);

/**
 * Adapt rho value based on current unscaled primal/dual residuals
 * @param work Workspace
 * @return     Exitflag
 */
int   adapt_rho(OSQPWorkspace *work);

/**
 * Set values of rho vector based on constraint types
 * @param work Workspace
 */
void    set_rho_vec(OSQPWorkspace *work);

/**
 * Update values of rho vector based on updated constraints.
 * If the constraints change, update the linear systems solver.
 *
 * @param work Workspace
 * @return     Exitflag
 */
int   update_rho_vec(OSQPWorkspace *work);


/**
 * Swap float vector pointers
 * @param a first vector
 * @param b second vector
 */
void swap_vectors(float **a,
                  float **b);


/**
 * Cold start workspace variables xz and y
 * @param work Workspace
 */
void cold_start(OSQPWorkspace *work);


/**
 * Update x_tilde and z_tilde variable (first ADMM step)
 * @param work [description]
 */
void update_xz_tilde(OSQPWorkspace *work);


/**
 * Update x (second ADMM step)
 * Update also delta_x (For for dual infeasibility)
 * @param work Workspace
 */
void update_x(OSQPWorkspace *work);


/**
 * Update z (third ADMM step)
 * @param work Workspace
 */
void update_z(OSQPWorkspace *work);


/**
 * Update y variable (fourth ADMM step)
 * Update also delta_y to check for primal infeasibility
 * @param work Workspace
 */
void update_y(OSQPWorkspace *work);


/**
 * Compute objective function from data at value x
 * @param  work OSQPWorkspace structure
 * @param  x    Value x
 * @return      Objective function value
 */
float compute_obj_val(OSQPWorkspace *work,
                        float       *x);

/**
 * Check whether QP has solution
 * @param info OSQPInfo
 */
int has_solution(OSQPInfo *info);

/**
 * Store the QP solution
 * @param work Workspace
 */
void store_solution(OSQPWorkspace *work);


/**
 * Update solver information
 * @param work               Workspace
 * @param iter               Iteration number
 * @param compute_objective  Boolean (if compute the objective or not)
 * @param polish             Boolean (if called from polish)
 */
void update_info(OSQPWorkspace *work,
                 int          iter,
                 int          compute_objective,
                 int          polish);


/**
 * Reset solver information (after problem updates)
 * @param info               Information structure
 */
void reset_info(OSQPInfo *info);


/**
 * Update solver status (value and string)
 * @param info OSQPInfo
 * @param status_val new status value
 */
void update_status(OSQPInfo *info,
                   int     status_val);


/**
 * Check if termination conditions are satisfied
 * If the boolean flag is ON, it checks for approximate conditions (10 x larger
 * tolerances than the ones set)
 *
 * @param  work        Workspace
 * @param  approximate Boolean
 * @return      Residuals check
 */
int check_termination(OSQPWorkspace *work,
                        int          approximate);



/**
 * Validate problem data
 * @param  data OSQPData to be validated
 * @return      Exitflag to check
 */
int validate_data(const OSQPData *data);


/**
 * Validate problem settings
 * @param  settings OSQPSettings to be validated
 * @return          Exitflag to check
 */
int validate_settings(const OSQPSettings *settings);


#endif // ifndef AUXIL_H
