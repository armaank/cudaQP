#ifndef SCALING_H
#define SCALING_H

// Functions to scale problem data
// # include "qptypes.h"
# include "linalg.h"
// # include "constants.h"

// Enable data scaling if EMBEDDED is disabled or if EMBEDDED == 2

/**
 * Scale problem matrices
 * @param  work Workspace
 * @return      exitflag
 */
int scale_data(qpWorkspace *work);


/**
 * Unscale problem matrices
 * @param  work Workspace
 * @return      exitflag
 */
int unscale_data(qpWorkspace *work);


/**
 * Unscale solution
 * @param  work Workspace
 * @return      exitflag
 */
int unscale_solution(qpWorkspace *work);

#endif SCALING_H
