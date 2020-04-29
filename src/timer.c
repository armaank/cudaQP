// #include "../include/timer.h"
#include "../include/qptypes.h"

void qp_tic(qpTimer *t)
{
    clock_gettime(CLOCK_MONOTONIC, &t->tic);
}

/* return time passed since last call to tic on this timer */
float qp_toc(qpTimer *t)
{
    struct timespec temp;

    clock_gettime(CLOCK_MONOTONIC, &t->toc);

    if ((t->toc.tv_nsec - t->tic.tv_nsec) < 0) {
        temp.tv_sec  = t->toc.tv_sec - t->tic.tv_sec - 1;
        temp.tv_nsec = 1e9 + t->toc.tv_nsec - t->tic.tv_nsec;
    } else {
        temp.tv_sec  = t->toc.tv_sec - t->tic.tv_sec;
        temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
    }
    return (float)temp.tv_sec + (float)temp.tv_nsec / 1e9;
}
