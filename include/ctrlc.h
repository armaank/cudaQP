/*
 * Interface for OSQP signal handling.
 */

#ifndef CTRLC_H
# define CTRLC_H

# include "glob_opts.h"
#  include <signal.h>


/* METHODS are the same for both */

/**
 * Start listener for ctrl-c interrupts
 */
void osqp_start_interrupt_listener(void);

/**
 * End listener for ctrl-c interrupts
 */
void osqp_end_interrupt_listener(void);

/**
 * Check if the solver has been interrupted
 * @return  Boolean indicating if the solver has been interrupted
 */
int osqp_is_interrupted(void);



#endif /* END IFDEF CTRLC */
