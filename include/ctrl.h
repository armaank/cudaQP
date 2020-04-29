#ifndef CTRL_H
#define CTRL_H
/*
 * Interface for OSQP signal handling.
 */



// # include "glob_opts.h"


/* Use sigaction for signal handling on non-Windows machines */
#  include <signal.h>


/* METHODS are the same for both */

/**
 * Start listener for ctrl-c interrupts
 */
void qp_start_interrupt_listener(void);

/**
 * End listener for ctrl-c interrupts
 */
void qp_end_interrupt_listener(void);

/**
 * Check if the solver has been interrupted
 * @return  Boolean indicating if the solver has been interrupted
 */
int qp_is_interrupted(void);

#endif //CTRL_H
