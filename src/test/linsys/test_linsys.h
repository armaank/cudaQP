#include <stdio.h>
#include <stdlib.h>

#include "../../../include/unittest.h"
#include "../../../include/linsys.h"
#include "../qptest.h"
#include "data.h"


#define TESTS_TOL 1e-4 // Define tests tolerance


static const char* test_solveKKT() {
    int m, exitflag = 0;
    float *rho_vec;
    LinSysSolver *s;  // Private structure to form KKT factorization
    qpParams *params = (qpParams *)malloc(sizeof(qpParams)); // Settings
    solve_linsys_sols_data *data = generate_problem_solve_linsys_sols_data();

    // Settings
    params->rho   = data->test_solve_KKT_rho;
    params->sigma = data->test_solve_KKT_sigma;

    // Set rho_vec
    m       = data->test_solve_KKT_A->m;
    rho_vec = (float*) calloc(m, sizeof(float));
    vec_add_scalar(rho_vec, params->rho, m);

    // Form and factorize KKT matrix
    exitflag = init_linsys_solver(&s, data->test_solve_KKT_Pu, data->test_solve_KKT_A,
                                  params->sigma, rho_vec, LINSYS_SOLVER, 0);

    // Solve  KKT x = b via LDL given factorization
    s->solve(s, data->test_solve_KKT_rhs);

    ut_assert(
        "Linear systems solve tests: error in forming and solving KKT system!",
        vec_norm_inf_diff(data->test_solve_KKT_rhs, data->test_solve_KKT_x,
                          data->test_solve_KKT_m + data->test_solve_KKT_n) < TESTS_TOL);


    // Cleanup
    s->free(s);
    free(params);
    free(rho_vec);
    clean_problem_solve_linsys_sols_data(data);

    return 0;
}


static const char* test_solve_linsys()
{
    ut_run_test(test_solveKKT);
    return 0;
}