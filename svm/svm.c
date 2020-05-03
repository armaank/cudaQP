/* svm */

#include <stdio.h>
#include <stdlib.h>
#include "../include/qp.h" 

// problem data, needs to be changed manually arg
#include "./n100m1000_1/data.h"

#define ERR 1e-4

int main(int argc, char **argv) {

    
    // Exitflag
    int exitflag = 0;

    // Workspace structures
    qpWorkspace *work;
    qpParams *params = (qpParams *)malloc(sizeof(qpParams));
    qpData *data = (qpData *)malloc(sizeof(qpData));
    
    data = generate_problem();

    // Define Solver params as default
    qp_set_default_params(params);
    params->max_iter   = 2000;
    params->alpha      = 1.6;
    params->polish     = 1;
    params->scaling    = 0;
    params->verbose    = 1;
    params->warm_start = 0;

    // Setup workspace
    exitflag = qp_setup(&work, data, params);

    // Solve Problem
    qp_solve(work);

    // Cleanup
    free(params);

    // Clean workspace
    qp_cleanup(work);

    return 0;
}




