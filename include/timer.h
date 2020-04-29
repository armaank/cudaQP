#ifndef TIMER_H
#define TIMER_H

#include <time.h>
#include <sys/time.h>
#include "qptypes.h"

struct QP_TIMER {
    struct timespec tic;
    struct timespec toc;
};

/* start timer */
void qp_tic(qpTimer *t);

/* report time */
float qp_toc(qpTimer *t);

#endif TIMER_H
