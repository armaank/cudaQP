#ifndef TIME_H
#define TIME_H

#include <time.h>
#include <sys/time.h>

struct qp_timer {
    struct timespec tic;
    struct timespec toc;
}

/* start timer */
void qp_tic(qp_timer *t);

/* report time */
float qp_toc(qp_timer *t);


#endif TIME_H