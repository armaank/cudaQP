/* Solution polish based on assuming the active set */
#ifndef POLISH_H
# define POLISH_H



# include "types.h"

/**
 * Solution polish: Solve equality constrained QP with assumed active
 *constraints
 * @param  work Workspace
 * @return      Exitflag
 */
c_int polish(OSQPWorkspace *work);


#endif // ifndef POLISH_H
