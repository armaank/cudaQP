#ifndef AUXIL_H
# define AUXIL_H

// # include "qptypes.h"
// #include "linalg.h"
// #include "constants.h"
#include "scaling.h"
#include "proj.h"
// #include "qp.h"
// #include "timer.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/**
 * Compute rho estimate from residuals
 * @param work Workspace
 * @return     rho estimate
 */
float compute_rho_estimate(qpWorkspace *work);

/**
 * Adapt rho value based on current unscaled primal/dual residuals
 * @param work Workspace
 * @return     Exitflag
 */
int   adapt_rho(qpWorkspace *work);

/**
 * Set values of rho vector based on constraint types
 * @param work Workspace
 */
void    set_rho_vec(qpWorkspace *work);

/**
 * Update values of rho vector based on updated constraints.
 * If the constraints change, update the linear systems solver.
 *
 * @param work Workspace
 * @return     Exitflag
 */
int   update_rho_vec(qpWorkspace *work);


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
void cold_start(qpWorkspace *work);


/**
 * Update x_tilde and z_tilde variable (first ADMM step)
 * @param work [description]
 */
void update_xz_tilde(qpWorkspace *work);


/**
 * Update x (second ADMM step)
 * Update also delta_x (For for dual infeasibility)
 * @param work Workspace
 */
void update_x(qpWorkspace *work);


/**
 * Update z (third ADMM step)
 * @param work Workspace
 */
void update_z(qpWorkspace *work);


/**
 * Update y variable (fourth ADMM step)
 * Update also delta_y to check for primal infeasibility
 * @param work Workspace
 */
void update_y(qpWorkspace *work);


/**
 * Compute objective function from data at value x
 * @param  work qpWorkspace structure
 * @param  x    Value x
 * @return      Objective function value
 */
float compute_obj_val(qpWorkspace *work,
                      float       *x);

/**
 * Check whether QP has solution
 * @param info qpInfo
 */
int has_solution(qpInfo *info);

/**
 * Store the QP solution
 * @param work Workspace
 */
void store_solution(qpWorkspace *work);


/**
 * Update solver information
 * @param work               Workspace
 * @param iter               Iteration number
 * @param compute_objective  Boolean (if compute the objective or not)
 * @param polish             Boolean (if called from polish)
 */
void update_info(qpWorkspace *work,
                 int          iter,
                 int          compute_objective,
                 int          polish);


/**
 * Reset solver information (after problem updates)
 * @param info               Information structure
 */
void reset_info(qpInfo *info);


/**
 * Update solver status (value and string)
 * @param info qpInfo
 * @param status_val new status value
 */
void update_status(qpInfo *info,
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
int check_termination(qpWorkspace *work,
                      int          approximate);

/**
 * Validate problem data
 * @param  data qpData to be validated
 * @return      Exitflag to check
 */
int validate_data(const qpData *data);


/**
 * Validate problem settings
 * @param  settings qpSettings to be validated
 * @return          Exitflag to check
 */
int validate_params(const qpParams *params);

#endif //AUXIL_H
