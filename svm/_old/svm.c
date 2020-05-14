/* svm */
#include "osqp.h"
#include <stdlib.h>
#include "data.h"

int main(void) {
  c_int exitflag, tmp_int;
  c_float tmp_float;
  csc *tmp_mat, *P_tmp;

  // Problem settings
  OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

  // Structures
  OSQPSolver   *solver; // Solver
  OSQPTestData *data;      // Data
  basic_qp_sols_data *sols_data;

  // Populate data
  data = generate_problem_basic_qp();
  sols_data = generate_problem_basic_qp_sols_data();

  // Define Solver settings as default
  osqp_set_default_settings(settings);
  settings->max_iter   = 2000;
  settings->alpha      = 1.6;
  settings->polish     = 1;
  settings->scaling    = 0;
  settings->verbose    = 1;
  settings->warm_start = 0;

  // Setup solver
  exitflag = osqp_setup(&solver, data->P, data->q,
                        data->A, data->l, data->u,
                        data->m, data->n, settings);
  // Solve Problem
  osqp_solve(solver);
    // Clean solver
  osqp_cleanup(solver);

  // Cleanup data
  clean_problem_basic_qp(data);
  clean_problem_basic_qp_sols_data(sols_data);

  // Cleanup
  c_free(settings);
  return 0;


}
