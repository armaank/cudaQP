#include "osqp.h"    // OSQP API
#include "minunit.h" // Basic testing script header

#include "svm/data.h"

static const char* test_svm_solve()
{
  c_int exitflag;

  // Problem settings
  OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

  // Structures
  OSQPSolver *solver; // Workspace
  OSQPTestData *data;      // Data

  // Populate data
  data = generate_problem_svm();

  // Define Solver settings as default
  osqp_set_default_settings(settings);
  settings->alpha   = 1.6;
  settings->rho     = 0.1;
  settings->polish  = 1;
  settings->verbose = 1;

  // Setup workspace
  exitflag = osqp_setup(&solver, data->P, data->q,
                        data->A, data->l, data->u,
                        data->m, data->n, settings);

  // Solve problem 
  osqp_solve(solver);

  // Clean workspace
  osqp_cleanup(solver);

  // Cleanup settings and data
  c_free(settings);
  clean_problem_svm(data);
  //clean_problem_basic_svm_sols_data(sols_data);

  return 0;
}

static const char* test_svm()
{
  mu_run_test(test_svm_solve);

  return 0;
}
