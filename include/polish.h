/* Solution polish based on assuming the active set */


# include "qptypes.h"
#include "auxil.h"
#include "lin_sys.h"
#include "kkt.h"
#include "proj.h"
#include "lin_alg.h"
#include "timer.h"

/**
 * Solution polish: Solve equality constrained QP with assumed active
 *constraints
 * @param  work Workspace
 * @return      Exitflag
 */
int polish(qpWorkspace *work);

