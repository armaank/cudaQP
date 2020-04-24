// Functions to scale problem data
# include "qptypes.h"
# include "lin_alg.h"
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
