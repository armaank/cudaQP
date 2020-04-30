// #include "qp.h"    // qp API
// #include "auxil.h"   // Needed for cold_start()
// #include "cs.h"      // CSC data structure
// #include "util.h"    // Utilities for testing
// #include "minunit.h" // Basic testing script header

// #include "basic_qp/data.h"
#include <stdio.h>
#include <stdlib.h>
#include "../../../include/qp.h"
#include "../../../include/unittest.h"
#include "../qptest.h"
#include "data.h"

#define TESTS_TOL 1e-4 // Define tests tolerance

static const char* test_basic_qp_solve()
{
  int exitflag, tmp_int;
  float tmp_float;
  csc *tmp_mat, *P_tmp;

  // Problem params
  qpParams *params = (qpParams *)malloc(sizeof(qpParams));

  // Structures
  qpWorkspace *work; // Workspace
  qpData *data;      // Data
  basic_qp_sols_data *sols_data;

  // Populate data
  data = generate_problem_basic_qp();
  sols_data = generate_problem_basic_qp_sols_data();

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

  // Setup correct
  ut_assert("Basic QP test solve: Setup error!", exitflag == 0);


  // Solve Problem
  qp_solve(work);

  // Compare solver statuses
  ut_assert("Basic QP test solve: Error in solver status!",
	    work->info->status_val == sols_data->status_test);

  // Compare primal solutions
  ut_assert("Basic QP test solve: Error in primal solution!",
	    vec_norm_inf_diff(work->solution->x, sols_data->x_test,
			      data->n) < TESTS_TOL);

  // Compare dual solutions
  ut_assert("Basic QP test solve: Error in dual solution!",
	    vec_norm_inf_diff(work->solution->y, sols_data->y_test,
			      data->m) < TESTS_TOL);


  // Compare objective values
  ut_assert("Basic QP test solve: Error in objective value!",
	    c_absval(work->info->obj_val - sols_data->obj_value_test) <
	    TESTS_TOL);

  // Try to set wrong params
  ut_assert("Basic QP test solve: Wrong value of rho not caught!",
	    qp_update_rho(work, -0.1) == 1);

  ut_assert("Basic QP test solve: Wrong value of max_iter not caught!",
	    qp_update_max_iter(work, -1) == 1);

  ut_assert("Basic QP test solve: Wrong value of eps_abs not caught!",
	    qp_update_eps_abs(work, -1.) == 1);

  ut_assert("Basic QP test solve: Wrong value of eps_rel not caught!",
	    qp_update_eps_rel(work, -1.) == 1);

  ut_assert("Basic QP test solve: Wrong value of eps_prim_inf not caught!",
	    qp_update_eps_prim_inf(work, -0.1) == 1);

  ut_assert("Basic QP test solve: Wrong value of eps_dual_inf not caught!",
	    qp_update_eps_dual_inf(work, -0.1) == 1);

  ut_assert("Basic QP test solve: Wrong value of alpha not caught!",
	    qp_update_alpha(work, 2.0) == 1);

  ut_assert("Basic QP test solve: Wrong value of warm_start not caught!",
	    qp_update_warm_start(work, -1) == 1);

  ut_assert("Basic QP test solve: Wrong value of scaled_termination not caught!",
	    qp_update_scaled_termination(work, 2) == 1);

  ut_assert("Basic QP test solve: Wrong value of check_termination not caught!",
	    qp_update_check_termination(work, -1) == 1);

  ut_assert("Basic QP test solve: Wrong value of delta not caught!",
	    qp_update_delta(work, 0.) == 1);

  ut_assert("Basic QP test solve: Wrong value of polish not caught!",
	    qp_update_polish(work, 2) == 1);

  ut_assert("Basic QP test solve: Wrong value of polish_refine_iter not caught!",
	    qp_update_polish_refine_iter(work, -1) == 1);

  ut_assert("Basic QP test solve: Wrong value of verbose not caught!",
	    qp_update_verbose(work, 2) == 1);

  // Clean workspace
  qp_cleanup(work);

  /* =============================
       SETUP WITH WRONG SETTINGS
     ============================= */

  // Setup workspace with empty params
  exitflag = qp_setup(&work, data, 0);
  ut_assert("Basic QP test solve: Setup should result in error due to empty params",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);

  // Setup workspace with a wrong number of scaling iterations
  tmp_int = params->scaling;
  params->scaling = -1;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to a negative number of scaling iterations",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->scaling = tmp_int;

  // Setup workspace with wrong params->adaptive_rho
  tmp_int = params->adaptive_rho;
  params->adaptive_rho = 2;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to non-boolean params->adaptive_rho",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->adaptive_rho = tmp_int;

  // Setup workspace with wrong params->adaptive_rho_interval
  tmp_int = params->adaptive_rho_interval;
  params->adaptive_rho_interval = -1;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to negative params->adaptive_rho_interval",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->adaptive_rho_interval = tmp_int;

  // Setup workspace with wrong params->adaptive_rho_fraction
  tmp_float = params->adaptive_rho_fraction;
  params->adaptive_rho_fraction = -1.5;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to non-positive params->adaptive_rho_fraction",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->adaptive_rho_fraction = tmp_float;

  // Setup workspace with wrong params->adaptive_rho_tolerance
  tmp_float = params->adaptive_rho_tolerance;
  params->adaptive_rho_tolerance = 0.5;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to wrong params->adaptive_rho_tolerance",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->adaptive_rho_tolerance = tmp_float;

  // Setup workspace with wrong params->polish_refine_iter
  tmp_int = params->polish_refine_iter;
  params->polish_refine_iter = -3;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to negative params->polish_refine_iter",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->polish_refine_iter = tmp_int;

  // Setup workspace with wrong params->rho
  tmp_float = params->rho;
  params->rho = 0.0;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to non-positive params->rho",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->rho = tmp_float;

  // Setup workspace with wrong params->sigma
  tmp_float = params->sigma;
  params->sigma = -0.1;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to non-positive params->sigma",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->sigma = tmp_float;

  // Setup workspace with wrong params->delta
  tmp_float = params->delta;
  params->delta = -1.1;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to non-positive params->delta",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->delta = tmp_float;

  // Setup workspace with wrong params->max_iter
  tmp_int = params->max_iter;
  params->max_iter = 0;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to non-positive params->max_iter",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->max_iter = tmp_int;

  // Setup workspace with wrong params->eps_abs
  tmp_float = params->eps_abs;
  params->eps_abs = -1.1;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to negative params->eps_abs",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->eps_abs = tmp_float;

  // Setup workspace with wrong params->eps_rel
  tmp_float = params->eps_rel;
  params->eps_rel = -0.1;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to negative params->eps_rel",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->eps_rel = tmp_float;

  // Setup workspace with wrong params->eps_prim_inf
  tmp_float = params->eps_prim_inf;
  params->eps_prim_inf = -0.1;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to non-positive params->eps_prim_inf",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->eps_prim_inf = tmp_float;

  // Setup workspace with wrong params->eps_dual_inf
  tmp_float = params->eps_dual_inf;
  params->eps_dual_inf = 0.0;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to non-positive params->eps_dual_inf",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->eps_dual_inf = tmp_float;

  // Setup workspace with wrong params->alpha
  tmp_float = params->alpha;
  params->alpha = 2.0;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to wrong params->alpha",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->alpha = tmp_float;

  // Setup workspace with wrong params->linsys_solver
  tmp_int = params->linsys_solver;
  params->linsys_solver = 5;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to wrong params->linsys_solver",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->linsys_solver = tmp_int;

  // Setup workspace with wrong params->verbose
  tmp_int = params->verbose;
  params->verbose = 2;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to non-boolean params->verbose",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->verbose = tmp_int;

  // Setup workspace with wrong params->scaled_termination
  tmp_int = params->scaled_termination;
  params->scaled_termination = 2;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to non-boolean params->scaled_termination",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->scaled_termination = tmp_int;

  // Setup workspace with wrong params->check_termination
  tmp_int = params->check_termination;
  params->check_termination = -1;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to non-boolean params->check_termination",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->check_termination = tmp_int;

  // Setup workspace with wrong params->warm_start
  tmp_int = params->warm_start;
  params->warm_start = 5;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to non-boolean params->warm_start",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->warm_start = tmp_int;

  // Setup workspace with wrong params->time_limit
  tmp_float = params->time_limit;
  params->time_limit = -0.2;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to wrong params->time_limit",
            exitflag == QP_SETTINGS_VALIDATION_ERROR);
  params->time_limit = tmp_float;


  /* =========================
       SETUP WITH WRONG DATA
     ========================= */

  // Setup workspace with empty data
  exitflag = qp_setup(&work, 0, params);
  ut_assert("Basic QP test solve: Setup should result in error due to empty data",
            exitflag == QP_DATA_VALIDATION_ERROR);

  // Setup workspace with wrong data->m
  tmp_int = data->m;
  data->m = data->m - 1;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to wrong data->m",
            exitflag == QP_DATA_VALIDATION_ERROR);
  data->m = tmp_int;

  // Setup workspace with wrong data->n
  tmp_int = data->n;
  data->n = data->n + 1;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to wrong data->n",
            exitflag == QP_DATA_VALIDATION_ERROR);

  // Setup workspace with zero data->n
  data->n = 0;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to zero data->n",
            exitflag == QP_DATA_VALIDATION_ERROR);
  data->n = tmp_int;

  // Setup workspace with wrong P->m
  tmp_int = data->P->m;
  data->P->m = data->n + 1;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to wrong P->m",
            exitflag == QP_DATA_VALIDATION_ERROR);
  data->P->m = tmp_int;

  // Setup workspace with wrong P->n
  tmp_int = data->P->n;
  data->P->n = data->n + 1;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to wrong P->n",
            exitflag == QP_DATA_VALIDATION_ERROR);
  data->P->n = tmp_int;

  // Setup workspace with non-upper-triangular P
  tmp_mat = data->P;

  // Construct non-upper-triangular P
  P_tmp = (csc*) malloc(sizeof(csc));
  P_tmp->m = 2;
  P_tmp->n = 2;
  P_tmp->nz = -1;
  P_tmp->nzmax = 4;
  P_tmp->x = (float*) malloc(4 * sizeof(float));
  P_tmp->x[0] = 4.0;
  P_tmp->x[1] = 1.0;
  P_tmp->x[2] = 1.0;
  P_tmp->x[3] = 2.0;
  P_tmp->i = (int*) malloc(4 * sizeof(int));
  P_tmp->i[0] = 0;
  P_tmp->i[1] = 1;
  P_tmp->i[2] = 0;
  P_tmp->i[3] = 1;
  P_tmp->p = (int*) malloc((2 + 1) * sizeof(int));
  P_tmp->p[0] = 0;
  P_tmp->p[1] = 2;
  P_tmp->p[2] = 4;

  data->P = P_tmp;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to non-triu structure of P",
            exitflag == QP_DATA_VALIDATION_ERROR);
  data->P = tmp_mat;

  // Setup workspace with non-consistent bounds
  data->l[0] = data->u[0] + 1.0;
  exitflag = qp_setup(&work, data, params);
  ut_assert("Basic QP test solve: Setup should result in error due to non-consistent bounds",
            exitflag == QP_DATA_VALIDATION_ERROR);


  // Cleanup data
  clean_problem_basic_qp(data);
  clean_problem_basic_qp_sols_data(sols_data);

  // Cleanup
  free(params);
  free(P_tmp->x);
  free(P_tmp->i);
  free(P_tmp->p);
  free(P_tmp);

  return 0;
}

static const char* test_basic_qp_update()
{
  int exitflag;

  // Problem params
  qpParams *params = (qpParams *)malloc(sizeof(qpParams));

  // Structures
  qpWorkspace *work; // Workspace
  qpData *data;      // Data
  basic_qp_sols_data *sols_data;

  // Populate data
  data = generate_problem_basic_qp();
  sols_data = generate_problem_basic_qp_sols_data();


  // Define Solver params as default
  qp_set_default_params(params);
  params->max_iter   = 200;
  params->alpha      = 1.6;
  params->polish     = 1;
  params->scaling    = 0;
  params->verbose    = 1;
  params->warm_start = 0;

  // Setup workspace
  exitflag = qp_setup(&work, data, params);

  // Setup correct
  ut_assert("Basic QP test update: Setup error!", exitflag == 0);


  // ====================================================================
  //  Update data
  // ====================================================================

  // Update linear cost
  qp_update_lin_cost(work, sols_data->q_new);
  ut_assert("Basic QP test update: Error in updating linear cost!",
            vec_norm_inf_diff(work->data->q, sols_data->q_new,
                              data->n) < TESTS_TOL);

  // UPDATE BOUND
  // Try to update with non-consistent values
  ut_assert("Basic QP test update: Error in bounds update ordering not caught!",
            qp_update_bounds(work, sols_data->u_new, sols_data->l_new) == 1);

  // Now update with correct values
  ut_assert("Basic QP test update: Error in bounds update ordering!",
            qp_update_bounds(work, sols_data->l_new, sols_data->u_new) == 0);

  ut_assert("Basic QP test update: Error in bounds update, lower bound!",
            vec_norm_inf_diff(work->data->l, sols_data->l_new,
                              data->m) < TESTS_TOL);

  ut_assert("Basic QP test update: Error in bounds update, upper bound!",
            vec_norm_inf_diff(work->data->u, sols_data->u_new,
                              data->m) < TESTS_TOL);

  // Return original values
  qp_update_bounds(work, data->l, data->u);


  // UPDATE LOWER BOUND
  // Try to update with non-consistent values
  ut_assert(
    "Basic QP test update: Error in lower bound update ordering not caught!",
    qp_update_lower_bound(work, sols_data->u_new) == 1);

  // Now update with correct values
  ut_assert("Basic QP test update: Error in lower bound update ordering!",
            qp_update_lower_bound(work, sols_data->l_new) == 0);

  ut_assert("Basic QP test update: Error in updating lower bound!",
            vec_norm_inf_diff(work->data->l, sols_data->l_new,
                              data->m) < TESTS_TOL);

  // Return original values
  qp_update_lower_bound(work, data->l);


  // UPDATE UPPER BOUND
  // Try to update with non-consistent values
  ut_assert(
    "Basic QP test update: Error in upper bound update: ordering not caught!",
    qp_update_upper_bound(work, sols_data->l_new) == 1);

  // Now update with correct values
  ut_assert("Basic QP test update: Error in upper bound update: ordering!",
            qp_update_upper_bound(work, sols_data->u_new) == 0);

  ut_assert("Basic QP test update: Error in updating upper bound!",
            vec_norm_inf_diff(work->data->u, sols_data->u_new,
                              data->m) < TESTS_TOL);


  // Clean workspace
  qp_cleanup(work);


  // Cleanup data
  clean_problem_basic_qp(data);
  clean_problem_basic_qp_sols_data(sols_data);

  // Cleanup
  free(params);

  return 0;
}

static const char* test_basic_qp_check_termination()
{
  int exitflag;

  // Problem params
  qpParams *params = (qpParams *)malloc(sizeof(qpParams));

  // Structures
  qpWorkspace *work; // Workspace
  qpData *data;      // Data
  basic_qp_sols_data *sols_data;

  // Populate data
  data = generate_problem_basic_qp();
  sols_data = generate_problem_basic_qp_sols_data();


  // Define Solver params as default
  qp_set_default_params(params);
  params->max_iter          = 200;
  params->alpha             = 1.6;
  params->polish            = 0;
  params->scaling           = 0;
  params->verbose           = 1;
  params->check_termination = 0;
  params->warm_start        = 0;

  // Setup workspace
  exitflag = qp_setup(&work, data, params);

  // Setup correct
  ut_assert("Basic QP test solve: Setup error!", exitflag == 0);

  // Solve Problem
  qp_solve(work);

  // Check if iter == max_iter
  ut_assert(
    "Basic QP test check termination: Error in number of iterations taken!",
    work->info->iter == work->params->max_iter);

  // Compare solver statuses
  ut_assert("Basic QP test check termination: Error in solver status!",
            work->info->status_val == sols_data->status_test);

  // Compare primal solutions
  ut_assert("Basic QP test check termination: Error in primal solution!",
            vec_norm_inf_diff(work->solution->x, sols_data->x_test,
                              data->n) < TESTS_TOL);

  // Compare dual solutions
  // print_vec(work->solution->y, data->m, "y_sol");
  // print_vec(sols_data->y_test, data->m, "y_test");
  ut_assert("Basic QP test check termination: Error in dual solution!",
            vec_norm_inf_diff(work->solution->y, sols_data->y_test,
                              data->m) < TESTS_TOL);

  // Compare objective values
  ut_assert("Basic QP test check termination: Error in objective value!",
            c_absval(work->info->obj_val - sols_data->obj_value_test) <
            TESTS_TOL);

  // Clean workspace
  qp_cleanup(work);

  // Cleanup data
  clean_problem_basic_qp(data);
  clean_problem_basic_qp_sols_data(sols_data);

  // Cleanup
  free(params);

  return 0;
}

static const char* test_basic_qp_update_rho()
{
  int extiflag;

  // Problem params
  qpParams *params = (qpParams *)malloc(sizeof(qpParams));

  // Structures
  qpWorkspace *work; // Workspace
  qpData *data;      // Data
  basic_qp_sols_data *sols_data;

  // Exitflag
  int exitflag;

  // rho to use
  float rho;

  // Define number of iterations to compare
  int n_iter_new_solver, n_iter_update_rho;

  // Populate data
  data = generate_problem_basic_qp();
  sols_data = generate_problem_basic_qp_sols_data();


  // Define Solver params as default
  rho = 0.7;
  qp_set_default_params(params);
  params->rho               = rho;
  params->adaptive_rho      = 0; // Disable adaptive rho for this test
  params->eps_abs           = 5e-05;
  params->eps_rel           = 5e-05;
  params->check_termination = 1;

  // Setup workspace
  exitflag = qp_setup(&work, data, params);

  // Setup correct
  ut_assert("Basic QP test update rho: Setup error!", exitflag == 0);

  // Solve Problem
  qp_solve(work);

  // Store number of iterations
  n_iter_new_solver = work->info->iter;

  // Compare solver statuses
  ut_assert("Update rho test solve: Error in solver status!",
            work->info->status_val == sols_data->status_test);

  // Compare primal solutions
  ut_assert("Update rho test solve: Error in primal solution!",
            vec_norm_inf_diff(work->solution->x, sols_data->x_test,
                              data->n)/vec_norm_inf(sols_data->x_test, data->n) < TESTS_TOL);

  // Compare dual solutions
  ut_assert("Update rho test solve: Error in dual solution!",
            vec_norm_inf_diff(work->solution->y, sols_data->y_test,
                              data->m)/vec_norm_inf(sols_data->y_test, data->m) < TESTS_TOL);

  // Compare objective values
  ut_assert("Update rho test solve: Error in objective value!",
            c_absval(work->info->obj_val - sols_data->obj_value_test) <
            TESTS_TOL);

  // Clean workspace
  qp_cleanup(work);


  // Create new problem with different rho and update it
  qp_set_default_params(params);
  params->rho               = 0.1;
  params->adaptive_rho      = 0;
  params->check_termination = 1;
  params->eps_abs           = 5e-05;
  params->eps_rel           = 5e-05;

  // Setup workspace
  exitflag = qp_setup(&work, data, params);

  // Setup correct
  ut_assert("Basic QP test update rho: Setup error!", exitflag == 0);

  // Update rho
  exitflag = qp_update_rho(work, rho);
  ut_assert("Basic QP test update rho: Error update rho!", exitflag == 0);

  // Solve Problem
  qp_solve(work);

  // Compare solver statuses
  ut_assert("Basic QP test update rho: Error in solver status!",
            work->info->status_val == sols_data->status_test);

  // Compare primal solutions
  ut_assert("Basic QP test update rho: Error in primal solution!",
            vec_norm_inf_diff(work->solution->x, sols_data->x_test,
                              data->n)/vec_norm_inf(sols_data->x_test, data->n) < TESTS_TOL);

  // Compare dual solutions
  ut_assert("Basic QP test update rho: Error in dual solution!",
            vec_norm_inf_diff(work->solution->y, sols_data->y_test,
                              data->m)/vec_norm_inf(sols_data->y_test, data->m)< TESTS_TOL);

  // Compare objective values
  ut_assert("Basic QP test update rho: Error in objective value!",
            c_absval(work->info->obj_val - sols_data->obj_value_test) <
            TESTS_TOL);

  // Get number of iterations
  n_iter_update_rho = work->info->iter;

  // Assert same number of iterations
  ut_assert("Basic QP test update rho: Error in number of iterations!",
            n_iter_new_solver == n_iter_update_rho);

  // Cleanup solver
  qp_cleanup(work);

  // Cleanup data
  clean_problem_basic_qp(data);
  clean_problem_basic_qp_sols_data(sols_data);

  // Cleanup
  free(params);

  return 0;
}

static const char* test_basic_qp_time_limit()
{
  int exitflag;

  // Problem params
  qpParams *params = (qpParams *)malloc(sizeof(qpParams));

  // Structures
  qpWorkspace *work; // Workspace
  qpData *data;      // Data
  basic_qp_sols_data *sols_data;

  // Populate data
  data = generate_problem_basic_qp();
  sols_data = generate_problem_basic_qp_sols_data();

  // Define Solver params as default
  qp_set_default_params(params);
  params->rho = 20;
  params->adaptive_rho = 0;

  // Check default time limit
  ut_assert("Basic QP test time limit: Default not correct", params->time_limit == 0);

  // Setup workspace
  exitflag = qp_setup(&work, data, params);

  // Setup correct
  ut_assert("Basic QP test time limit: Setup error!", exitflag == 0);

  // Solve Problem
  qp_solve(work);

  // Compare solver statuses
  ut_assert("Basic QP test time limit: Error in no time limit solver status!",
	    work->info->status_val == sols_data->status_test);

  // Update time limit
  qp_update_time_limit(work, 1e-5);
  qp_update_eps_rel(work, 1e-09);
  qp_update_eps_abs(work, 1e-09);
// # else
//   // Not printing makes the code run a lot faster, so we need to make it work harder
//   // to fail by time limit exceeded
//   qp_update_time_limit(work, 1e-7);
//   qp_update_eps_rel(work, 1e-12);
//   qp_update_eps_abs(work, 1e-12);
// # endif
  qp_update_max_iter(work, (int)2e9);
  qp_update_check_termination(work, 0);

  // Solve Problem
  cold_start(work);
  qp_solve(work);

  // Compare solver statuses
  ut_assert("Basic QP test time limit: Error in timed out solver status!",
	    work->info->status_val == qp_TIME_LIMIT_REACHED);

  // Cleanup solver
  qp_cleanup(work);

  // Cleanup data
  clean_problem_basic_qp(data);
  clean_problem_basic_qp_sols_data(sols_data);

  // Cleanup
  free(params);

  return 0;
}


static const char* test_basic_qp_warm_start()
{
  int exitflag, iter;

  // Cold started variables
  float x0[2] = { 0.0, 0.0, };
  float y0[4] = { 0.0, 0.0, 0.0, 0.0, };

  // Optimal solution
  float xopt[2] = { 0.3, 0.7, };
  float yopt[4] = {-2.9, 0.0, 0.2, 0.0, };

  // Problem params
  qpParams *params = (qpParams *)malloc(sizeof(qpParams));

  // Structures
  qpWorkspace *work; // Workspace
  qpData *data;      // Data
  basic_qp_sols_data *sols_data;

  // Populate data
  data = generate_problem_basic_qp();
  sols_data = generate_problem_basic_qp_sols_data();

  // Define Solver params as default
  qp_set_default_params(params);
  params->check_termination = 1;

  // Setup workspace
  exitflag = qp_setup(&work, data, params);

  // Solve Problem
  qp_solve(work);
  iter = work->info->iter;

  // Cold start and solve again
  qp_warm_start(work, x0, y0);
  qp_solve(work);

  // Check if the number of iterations is the same
  ut_assert("Basic QP test warm start: Cold start error!", work->info->iter == iter);

  // Warm start from the solution and solve again
  qp_warm_start_x(work, xopt);
  qp_warm_start_y(work, yopt);
  qp_solve(work);

  // Check that the number of iterations equals 1
  ut_assert("Basic QP test warm start: Warm start error!", work->info->iter == 1);

  // Cleanup solver
  qp_cleanup(work);

  // Cleanup data
  clean_problem_basic_qp(data);
  clean_problem_basic_qp_sols_data(sols_data);

  // Cleanup
  free(params);

  return 0;
}


static const char* test_basic_qp()
{
  ut_run_test(test_basic_qp_solve);
// #ifdef ENABLE_MKL_PARDISO
//   ut_run_test(test_basic_qp_solve_pardiso);
// #endif
  ut_run_test(test_basic_qp_update);
  ut_run_test(test_basic_qp_check_termination);
  ut_run_test(test_basic_qp_update_rho);
// #ifdef PROFILING
//   ut_run_test(test_basic_qp_time_limit);
// #endif
  ut_run_test(test_basic_qp_warm_start);

  return 0;
}
