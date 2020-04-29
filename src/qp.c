#include "../include/qp.h"
// #include "auxil.h"
// #include "util.h"
// #include "scaling.h"
// #include "glob_opts.h"
// #include "error.h"


// #ifndef EMBEDDED
// # include "polish.h"
// #endif /* ifndef EMBEDDED */

// #ifdef CTRLC
// # include "ctrlc.h"
// #endif /* ifdef CTRLC */

// #ifndef EMBEDDED
// # include "lin_sys.h"
// #endif /* ifndef EMBEDDED */

/**********************
* Main API Functions *
**********************/
static void print_line(void) {

  int HEADER_LINE_LEN = 65;
  char  the_line[HEADER_LINE_LEN + 1];
  int i;

  for (i = 0; i < HEADER_LINE_LEN; ++i) the_line[i] = '-';
  the_line[HEADER_LINE_LEN] = '\0';
  printf("%s\n", the_line);
}

void print_footer(qpInfo *info, int polish) {
  printf("\n"); // Add space after iterations

  printf("status:               %s\n", info->status);

  if (polish && (info->status_val == qp_SOLVED)) {
    if (info->status_polish == 1) {
      printf("solution polish:      successful\n");
    } else if (info->status_polish < 0) {
      printf("solution polish:      unsuccessful\n");
    }
  }

  printf("number of iterations: %i\n", (int)info->iter);

  if ((info->status_val == qp_SOLVED) ||
      (info->status_val == qp_SOLVED_INACCURATE)) {
    printf("optimal objective:    %.4f\n", info->obj_val);
  }

  printf("run time:             %.2es\n", info->run_time);

  printf("optimal rho estimate: %.2e\n", info->rho_estimate);
  printf("\n");
}
void print_header(void) {
  printf("iter   ");


  // Main information
  printf("objective    pri res    dua res    rho");
  printf("        time");
  printf("\n");
}
void print_setup_header(const qpWorkspace *work) {
  qpData *data;
  qpParams *params;
  int nnz; // Number of nonzeros in the problem

  data     = work->data;
  params = work->params;

  // Number of nonzeros
  nnz = data->P->p[data->P->n] + data->A->p[data->A->n];

  print_line();
  printf("           OSQP v%d  -  Operator Splitting QP Solver\n"
          "              (c) Bartolomeo Stellato,  Goran Banjac\n"
          "        University of Oxford  -  Stanford University 2019\n",
          26);
  print_line();

  // Print variables and constraints
  printf("problem:  ");
  printf("variables n = %i, constraints m = %i\n          ",
                                    (int)data->n,
          (int)data->m);
  printf("nnz(P) + nnz(A) = %i\n", (int)nnz);

  // Print Settings
  printf("params: ");
  printf("linear system solver = %s",
          LINSYS_SOLVER_NAME[params->linsys_solver]);

  if (work->linsys_solver->nthreads != 1) {
    printf(" (%d threads)", (int)work->linsys_solver->nthreads);
  }
  printf(",\n          ");

  printf("eps_abs = %.1e, eps_rel = %.1e,\n          ",
          params->eps_abs, params->eps_rel);
  printf("eps_prim_inf = %.1e, eps_dual_inf = %.1e,\n          ",
          params->eps_prim_inf, params->eps_dual_inf);
  printf("rho = %.2e ", params->rho);

  if (params->adaptive_rho) printf("(adaptive)");
  printf(",\n          ");
  printf("sigma = %.2e, alpha = %.2f, ",
          params->sigma, params->alpha);
  printf("max_iter = %i\n", (int)params->max_iter);

  if (params->check_termination) printf(
      "          check_termination: on (interval %i),\n",
      (int)params->check_termination);
  else printf("          check_termination: off,\n");

  if (params->time_limit) printf("          time_limit: %.2e sec,\n",
                                    params->time_limit);

  if (params->scaling) printf("          scaling: on, ");
  else printf("          scaling: off, ");

  if (params->scaled_termination) printf("scaled_termination: on\n");
  else printf("scaled_termination: off\n");

  if (params->warm_start) printf("          warm start: on, ");
  else printf("          warm start: off, ");

  if (params->polish) printf("polish: on, ");
  else printf("polish: off, ");

  if (params->time_limit) printf("time_limit: %.2e sec\n", params->time_limit);
  else printf("time_limit: off\n");

  printf("\n");
}


void print_summary(qpWorkspace *work) {
  qpInfo *info;

  info = work->info;

  printf("%4i",     (int)info->iter);
  printf(" %12.4e", info->obj_val);
  printf("  %9.2e", info->pri_res);
  printf("  %9.2e", info->dua_res);
  printf("  %9.2e", work->params->rho);

  if (work->first_run) {
    // total time: setup + solve
    printf("  %9.2es", info->setup_time + info->solve_time);
  } else {
    // total time: update + solve
    printf("  %9.2es", info->update_time + info->solve_time);
  }
  printf("\n");

  work->summary_printed = 1; // Summary has been printed
}



void qp_set_default_params(qpParams *params) {

    params->rho           = (float)RHO;            /* ADMM step */
    params->sigma         = (float)SIGMA;          /* ADMM step */
    params->scaling = SCALING;                       /* heuristic problem scaling */

    params->adaptive_rho           = ADAPTIVE_RHO;
    params->adaptive_rho_interval  = ADAPTIVE_RHO_INTERVAL;
    params->adaptive_rho_tolerance = (float)ADAPTIVE_RHO_TOLERANCE;

    params->adaptive_rho_fraction = (float)ADAPTIVE_RHO_FRACTION;


    params->max_iter      = MAX_ITER;                /* maximum iterations to
                                                        take */
    params->eps_abs       = (float)EPS_ABS;        /* absolute convergence
                                                        tolerance */
    params->eps_rel       = (float)EPS_REL;        /* relative convergence
                                                        tolerance */
    params->eps_prim_inf  = (float)EPS_PRIM_INF;   /* primal infeasibility
                                                        tolerance */
    params->eps_dual_inf  = (float)EPS_DUAL_INF;   /* dual infeasibility
                                                        tolerance */
    params->alpha         = (float)ALPHA;          /* relaxation parameter */
    params->linsys_solver = LINSYS_SOLVER;           /* relaxation parameter */

    params->delta              = DELTA;              /* regularization parameter
                                                        for polish */
    params->polish             = POLISH;             /* ADMM solution polish: 1
                                                      */
    params->polish_refine_iter = POLISH_REFINE_ITER; /* iterative refinement
                                                        steps in polish */
    params->verbose            = VERBOSE;            /* print output */

    params->scaled_termination = SCALED_TERMINATION; /* Evaluate scaled
                                                        termination criteria*/
    params->check_termination  = CHECK_TERMINATION;  /* Interval for evaluating
                                                        termination criteria */
    params->warm_start         = WARM_START;         /* warm starting */

    params->time_limit = TIME_LIMIT;
}

qpParams* copy_params(const qpParams *params) {
  qpParams *new = malloc(sizeof(qpParams));

  if (!new) return 0;

  // Copy params
  // NB. Copying them explicitly because memcpy is not
  // defined when PRINTING is disabled (appears in string.h)
  new->rho = params->rho;
  new->sigma = params->sigma;
  new->scaling = params->scaling;

  new->adaptive_rho = params->adaptive_rho;
  new->adaptive_rho_interval = params->adaptive_rho_interval;
  new->adaptive_rho_tolerance = params->adaptive_rho_tolerance;
  new->adaptive_rho_fraction = params->adaptive_rho_fraction;

  new->max_iter = params->max_iter;
  new->eps_abs = params->eps_abs;
  new->eps_rel = params->eps_rel;
  new->eps_prim_inf = params->eps_prim_inf;
  new->eps_dual_inf = params->eps_dual_inf;
  new->alpha = params->alpha;
  new->linsys_solver = params->linsys_solver;
  new->delta = params->delta;
  new->polish = params->polish;
  new->polish_refine_iter = params->polish_refine_iter;
  new->verbose = params->verbose;
  new->scaled_termination = params->scaled_termination;
  new->check_termination = params->check_termination;
  new->warm_start = params->warm_start;
  new->time_limit = params->time_limit;

  return new;
}


int qp_setup(qpWorkspace** workp, const qpData *data, const qpParams *params) {
    int exitflag;

    qpWorkspace * work;

    // Validate data
    if (validate_data(data)) return qp_error(QP_DATA_VALIDATION_ERROR);

    // Validate params
    if (validate_params(params)) return qp_error(QP_SETTINGS_VALIDATION_ERROR);

    // Allocate empty workspace
    work = calloc(1, sizeof(qpWorkspace));
    if (!(work)) return qp_error(QP_MEM_ALLOC_ERROR);
    *workp = work;

    // Start and allocate directly timer
    work->timer = malloc(sizeof(qpTimer));
    if (!(work->timer)) return qp_error(QP_MEM_ALLOC_ERROR);
    qp_tic(work->timer);

    // Copy problem data into workspace
    work->data = malloc(sizeof(qpData));
    if (!(work->data)) return qp_error(QP_MEM_ALLOC_ERROR);
    work->data->n = data->n;
    work->data->m = data->m;

    // Cost function
    work->data->P = copy_csc_mat(data->P);
    work->data->q = vec_copy(data->q, data->n);
    if (!(work->data->P) || !(work->data->q)) return qp_error(QP_MEM_ALLOC_ERROR);

    // Constraints
    work->data->A = copy_csc_mat(data->A);
    if (!(work->data->A)) return qp_error(QP_MEM_ALLOC_ERROR);
    work->data->l = vec_copy(data->l, data->m);
    work->data->u = vec_copy(data->u, data->m);
    if ( data->m && (!(work->data->l) || !(work->data->u)) )
        return qp_error(QP_MEM_ALLOC_ERROR);

    // Vectorized rho parameter
    work->rho_vec     = malloc(data->m * sizeof(float));
    work->rho_inv_vec = malloc(data->m * sizeof(float));
    if ( data->m && (!(work->rho_vec) || !(work->rho_inv_vec)) )
        return qp_error(QP_MEM_ALLOC_ERROR);

    // Type of constraints
    work->constr_type = calloc(data->m, sizeof(int));
    if (data->m && !(work->constr_type)) return qp_error(QP_MEM_ALLOC_ERROR);

    // Allocate internal solver variables (ADMM steps)
    work->x        = calloc(data->n, sizeof(float));
    work->z        = calloc(data->m, sizeof(float));
    work->xz_tilde = calloc(data->n + data->m, sizeof(float));
    work->x_prev   = calloc(data->n, sizeof(float));
    work->z_prev   = calloc(data->m, sizeof(float));
    work->y        = calloc(data->m, sizeof(float));
    if (!(work->x) || !(work->xz_tilde) || !(work->x_prev))
        return qp_error(QP_MEM_ALLOC_ERROR);
    if ( data->m && (!(work->z) || !(work->z_prev) || !(work->y)) )
        return qp_error(QP_MEM_ALLOC_ERROR);

    // Initialize variables x, y, z to 0
    cold_start(work);

    // Primal and dual residuals variables
    work->Ax  = calloc(data->m, sizeof(float));
    work->Px  = calloc(data->n, sizeof(float));
    work->Aty = calloc(data->n, sizeof(float));

    // Primal infeasibility variables
    work->delta_y   = calloc(data->m, sizeof(float));
    work->Atdelta_y = calloc(data->n, sizeof(float));

    // Dual infeasibility variables
    work->delta_x  = calloc(data->n, sizeof(float));
    work->Pdelta_x = calloc(data->n, sizeof(float));
    work->Adelta_x = calloc(data->m, sizeof(float));

    if (!(work->Px) || !(work->Aty) || !(work->Atdelta_y) ||
            !(work->delta_x) || !(work->Pdelta_x))
        return qp_error(QP_MEM_ALLOC_ERROR);
    if ( data->m && (!(work->Ax) || !(work->delta_y) || !(work->Adelta_x)) )
        return qp_error(QP_MEM_ALLOC_ERROR);

    // Copy params
    work->params = copy_params(params);
    if (!(work->params)) return qp_error(QP_MEM_ALLOC_ERROR);

    // Perform scaling
    if (params->scaling) {
        // Allocate scaling structure
        work->scaling = malloc(sizeof(qpScaling));
        if (!(work->scaling)) return qp_error(QP_MEM_ALLOC_ERROR);
        work->scaling->D    = malloc(data->n * sizeof(float));
        work->scaling->Dinv = malloc(data->n * sizeof(float));
        work->scaling->E    = malloc(data->m * sizeof(float));
        work->scaling->Einv = malloc(data->m * sizeof(float));
        if (!(work->scaling->D) || !(work->scaling->Dinv))
            return qp_error(QP_MEM_ALLOC_ERROR);
        if ( data->m && (!(work->scaling->E) || !(work->scaling->Einv)) )
            return qp_error(QP_MEM_ALLOC_ERROR);


        // Allocate workspace variables used in scaling
        work->D_temp   = malloc(data->n * sizeof(float));
        work->D_temp_A = malloc(data->n * sizeof(float));
        work->E_temp   = malloc(data->m * sizeof(float));
        // if (!(work->D_temp) || !(work->D_temp_A) || !(work->E_temp))
        //   return qp_error(QP_MEM_ALLOC_ERROR);
        if (!(work->D_temp) || !(work->D_temp_A)) return qp_error(QP_MEM_ALLOC_ERROR);
        if (data->m && !(work->E_temp))           return qp_error(QP_MEM_ALLOC_ERROR);

        // Scale data
        scale_data(work);
    } else {
        work->scaling  = 0;
        work->D_temp   = 0;
        work->D_temp_A = 0;
        work->E_temp   = 0;
    }

    // Set type of constraints
    set_rho_vec(work);

    // Load linear system solver
    if (load_linsys_solver(work->params->linsys_solver)) return qp_error(QP_LINSYS_SOLVER_LOAD_ERROR);

    // Initialize linear system solver structure
    exitflag = init_linsys_solver(&(work->linsys_solver), work->data->P, work->data->A,
                                  work->params->sigma, work->rho_vec,
                                  work->params->linsys_solver, 0);

    if (exitflag) {
        return qp_error(exitflag);
    }

    // Initialize active constraints structure
    work->pol = malloc(sizeof(qpPolish));
    if (!(work->pol)) return qp_error(QP_MEM_ALLOC_ERROR);
    work->pol->Alow_to_A = malloc(data->m * sizeof(int));
    work->pol->Aupp_to_A = malloc(data->m * sizeof(int));
    work->pol->A_to_Alow = malloc(data->m * sizeof(int));
    work->pol->A_to_Aupp = malloc(data->m * sizeof(int));
    work->pol->x         = malloc(data->n * sizeof(float));
    work->pol->z         = malloc(data->m * sizeof(float));
    work->pol->y         = malloc(data->m * sizeof(float));
    if (!(work->pol->x)) return qp_error(QP_MEM_ALLOC_ERROR);
    if ( data->m && (!(work->pol->Alow_to_A) || !(work->pol->Aupp_to_A) ||
                     !(work->pol->A_to_Alow) || !(work->pol->A_to_Aupp) ||
                     !(work->pol->z) || !(work->pol->y)) )
        return qp_error(QP_MEM_ALLOC_ERROR);

    // Allocate solution
    work->solution = calloc(1, sizeof(qpSolution));
    if (!(work->solution)) return qp_error(QP_MEM_ALLOC_ERROR);
    work->solution->x = calloc(1, data->n * sizeof(float));
    work->solution->y = calloc(1, data->m * sizeof(float));
    if (!(work->solution->x))            return qp_error(QP_MEM_ALLOC_ERROR);
    if (data->m && !(work->solution->y)) return qp_error(QP_MEM_ALLOC_ERROR);

    // Allocate and initialize information
    work->info = calloc(1, sizeof(qpInfo));
    if (!(work->info)) return qp_error(QP_MEM_ALLOC_ERROR);
    work->info->status_polish = 0;              // Polishing not performed
    update_status(work->info, qp_UNSOLVED);

    work->info->solve_time  = 0.0;                   // Solve time to zero
    work->info->update_time = 0.0;                   // Update time to zero
    work->info->polish_time = 0.0;                   // Polish time to zero
    work->info->run_time    = 0.0;                   // Total run time to zero
    work->info->setup_time  = qp_toc(work->timer); // Update timer information

    work->first_run         = 1;
    work->clear_update_time = 0;
    work->rho_update_from_solve = 0;

    work->info->rho_updates  = 0;                    // Rho updates set to 0
    work->info->rho_estimate = work->params->rho;  // Best rho estimate

    // Print header

    if (work->params->verbose) print_setup_header(work);
    work->summary_printed = 0; // Initialize last summary  to not printed



    // If adaptive rho and automatic interval, but profiling disabled, we need to
    // set the interval to a default value

    if (work->params->adaptive_rho && !work->params->adaptive_rho_interval) {
        if (work->params->check_termination) {
            // If check_termination is enabled, we set it to a multiple of the check
            // termination interval
            work->params->adaptive_rho_interval = ADAPTIVE_RHO_MULTIPLE_TERMINATION *
                                                  work->params->check_termination;
        } else {
            // If check_termination is disabled we set it to a predefined fix number
            work->params->adaptive_rho_interval = ADAPTIVE_RHO_FIXED;
        }
    }


    // Return exit flag
    return 0;
}




int qp_solve(qpWorkspace *work) {

    int exitflag;
    int iter;
    int compute_cost_function; // Boolean: compute the cost function in the loop or not
    int can_check_termination; // Boolean: check termination or not


    float temp_run_time;       // Temporary variable to store current run time



    int can_print;             // Boolean whether you can print


    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);


    if (work->clear_update_time == 1)
        work->info->update_time = 0.0;
    work->rho_update_from_solve = 1;


    // Initialize variables
    exitflag              = 0;
    can_check_termination = 0;

    can_print = work->params->verbose;


    compute_cost_function = work->params->verbose; // Compute cost function only
    // if verbose is on

    compute_cost_function = 0;                       // Never compute cost
    // function during the
    // iterations if no printing
    // enabled





    qp_tic(work->timer); // Start timer





    if (work->params->verbose) {
        // Print Header for every column
        print_header();
        // printf("header");
    }




    // initialize Ctrl-C support
    qp_start_interrupt_listener();


    // Initialize variables (cold start or warm start depending on params)
    if (!work->params->warm_start) cold_start(work);  // If not warm start ->
    // set x, z, y to zero

    // Main ADMM algorithm
    for (iter = 1; iter <= work->params->max_iter; iter++) {
        // Update x_prev, z_prev (preallocated, no malloc)
        swap_vectors(&(work->x), &(work->x_prev));
        swap_vectors(&(work->z), &(work->z_prev));

        /* ADMM STEPS */
        /* Compute \tilde{x}^{k+1}, \tilde{z}^{k+1} */
        update_xz_tilde(work);

        /* Compute x^{k+1} */
        update_x(work);

        /* Compute z^{k+1} */
        update_z(work);

        /* Compute y^{k+1} */
        update_y(work);

        /* End of ADMM Steps */



        // Check the interrupt signal
        if (qp_is_interrupted()) {
            update_status(work->info, qp_SIGINT);

            printf("Solver interrupted\n");

            exitflag = 1;
            goto exit;
        }


        // Check if solver time_limit is enabled. In case, check if the current
        // run time is more than the time_limit option.
        if (work->first_run) {
            temp_run_time = work->info->setup_time + qp_toc(work->timer);
        }
        else {
            temp_run_time = work->info->update_time + qp_toc(work->timer);
        }

        if (work->params->time_limit &&
                (temp_run_time >= work->params->time_limit)) {
            update_status(work->info, qp_TIME_LIMIT_REACHED);

            if (work->params->verbose) printf("run time limit reached\n");
            can_print = 0;  // Not printing at this iteration

            break;
        }



        // Can we check for termination ?
        can_check_termination = work->params->check_termination &&
                                (iter % work->params->check_termination == 0);



        // Can we print ?
        can_print = work->params->verbose &&
                    ((iter % PRINT_INTERVAL == 0) || (iter == 1));

        if (can_check_termination || can_print) { // Update status in either of
            // these cases
            // Update information
            update_info(work, iter, compute_cost_function, 0);

            if (can_print) {
                // Print summary
                print_summary(work);
                // printf("work");
            }

            if (can_check_termination) {
                // Check algorithm termination
                if (check_termination(work, 0)) {
                    // Terminate algorithm
                    break;
                }
            }
        }


        if (can_check_termination) {
            // Update information and compute also objective value
            update_info(work, iter, compute_cost_function, 0);

            // Check algorithm termination
            if (check_termination(work, 0)) {
                // Terminate algorithm
                break;
            }
        }





        // If adaptive rho with automatic interval, check if the solve time is a
        // certain fraction
        // of the setup time.
        if (work->params->adaptive_rho && !work->params->adaptive_rho_interval) {
            // Check time
            if (qp_toc(work->timer) >
                    work->params->adaptive_rho_fraction * work->info->setup_time) {
                // Enough time has passed. We now get the number of iterations between
                // the updates.
                if (work->params->check_termination) {
                    // If check_termination is enabled, we round the number of iterations
                    // between
                    // rho updates to the closest multiple of check_termination
                    work->params->adaptive_rho_interval = (int)c_roundmultiple(iter,
                                                          work->params->check_termination);
                } else {
                    // If check_termination is disabled, we round the number of iterations
                    // between
                    // updates to the closest multiple of the default check_termination
                    // interval.
                    work->params->adaptive_rho_interval = (int)c_roundmultiple(iter,
                                                          CHECK_TERMINATION);
                }

                // Make sure the interval is not 0 and at least check_termination times
                work->params->adaptive_rho_interval = c_max(
                        work->params->adaptive_rho_interval,
                        work->params->check_termination);
            } // If time condition is met
        }   // If adaptive rho enabled and interval set to auto


        // Adapt rho
        if (work->params->adaptive_rho &&
                work->params->adaptive_rho_interval &&
                (iter % work->params->adaptive_rho_interval == 0)) {
            // Update info with the residuals if it hasn't been done before


            if (!can_check_termination && !can_print) {
                // Information has not been computed neither for termination or printing
                // reasons
                update_info(work, iter, compute_cost_function, 0);
            }


            if (!can_check_termination) {
                // Information has not been computed before for termination check
                update_info(work, iter, compute_cost_function, 0);
            }


            // Actually update rho
            if (adapt_rho(work)) {

                printf("Failed rho update");

                exitflag = 1;
                goto exit;
            }
        }


    }        // End of ADMM for loop


    // Update information and check termination condition if it hasn't been done
    // during last iteration (max_iter reached or check_termination disabled)
    if (!can_check_termination) {
        /* Update information */


        if (!can_print) {
            // Update info only if it hasn't been updated before for printing
            // reasons
            update_info(work, iter - 1, compute_cost_function, 0);
        }


        // If no printing is enabled, update info directly
        update_info(work, iter - 1, compute_cost_function, 0);




        /* Print summary */
        if (work->params->verbose && !work->summary_printed) print_summary(work);


        /* Check whether a termination criterion is triggered */
        check_termination(work, 0);
    }

    // Compute objective value in case it was not
    // computed during the iterations
    if (!compute_cost_function && has_solution(work->info)) {
        work->info->obj_val = compute_obj_val(work, work->x);
    }



    /* Print summary for last iteration */
    if (work->params->verbose && !work->summary_printed) {
        print_summary(work);
        // printf("work");
    }


    /* if max iterations reached, change status accordingly */
    if (work->info->status_val == qp_UNSOLVED) {
        if (!check_termination(work, 1)) { // Try to check for approximate
            update_status(work->info, qp_MAX_ITER_REACHED);
        }
    }


    /* if time-limit reached check termination and update status accordingly */
    if (work->info->status_val == qp_TIME_LIMIT_REACHED) {
        if (!check_termination(work, 1)) { // Try for approximate solutions
            update_status(work->info, qp_TIME_LIMIT_REACHED); /* Change update status back to qp_TIME_LIMIT_REACHED */
        }
    }




    /* Update rho estimate */
    work->info->rho_estimate = compute_rho_estimate(work);


    /* Update solve time */

    work->info->solve_time = qp_toc(work->timer);




    // Polish the obtained solution
    if (work->params->polish && (work->info->status_val == qp_SOLVED))
        polish(work);



    /* Update total time */
    if (work->first_run) {
        // total time: setup + solve + polish
        work->info->run_time = work->info->setup_time +
                               work->info->solve_time +
                               work->info->polish_time;
    } else {
        // total time: update + solve + polish
        work->info->run_time = work->info->update_time +
                               work->info->solve_time +
                               work->info->polish_time;
    }

    // Indicate that the solve function has already been executed
    if (work->first_run) work->first_run = 0;

    // Indicate that the update_time should be set to zero
    work->clear_update_time = 1;

    // Indicate that qp_update_rho is not called from qp_solve
    work->rho_update_from_solve = 0;



    /* Print final footer */
    if (work->params->verbose) print_footer(work->info, work->params->polish);


    // Store solution
    store_solution(work);


// Define exit flag for quitting function

exit:



    // Restore previous signal handler
    qp_end_interrupt_listener();


    return exitflag;
}




int qp_cleanup(qpWorkspace *work) {
    int exitflag = 0;

    if (work) { // If workspace has been allocated
        // Free Data
        if (work->data) {
            if (work->data->P) csc_spfree(work->data->P);
            if (work->data->A) csc_spfree(work->data->A);
            if (work->data->q) free(work->data->q);
            if (work->data->l) free(work->data->l);
            if (work->data->u) free(work->data->u);
            free(work->data);
        }

        // Free scaling variables
        if (work->scaling) {
            if (work->scaling->D)    free(work->scaling->D);
            if (work->scaling->Dinv) free(work->scaling->Dinv);
            if (work->scaling->E)    free(work->scaling->E);
            if (work->scaling->Einv) free(work->scaling->Einv);
            free(work->scaling);
        }

        // Free temp workspace variables for scaling
        if (work->D_temp)   free(work->D_temp);
        if (work->D_temp_A) free(work->D_temp_A);
        if (work->E_temp)   free(work->E_temp);

        // Free linear system solver structure
        if (work->linsys_solver) {
            if (work->linsys_solver->free) {
                work->linsys_solver->free(work->linsys_solver);
            }
        }

        // Unload linear system solver after free
        if (work->params) {
            exitflag = unload_linsys_solver(work->params->linsys_solver);
        }


        // Free active constraints structure
        if (work->pol) {
            if (work->pol->Alow_to_A) free(work->pol->Alow_to_A);
            if (work->pol->Aupp_to_A) free(work->pol->Aupp_to_A);
            if (work->pol->A_to_Alow) free(work->pol->A_to_Alow);
            if (work->pol->A_to_Aupp) free(work->pol->A_to_Aupp);
            if (work->pol->x)         free(work->pol->x);
            if (work->pol->z)         free(work->pol->z);
            if (work->pol->y)         free(work->pol->y);
            free(work->pol);
        }


        // Free other Variables
        if (work->rho_vec)     free(work->rho_vec);
        if (work->rho_inv_vec) free(work->rho_inv_vec);

        if (work->constr_type) free(work->constr_type);

        if (work->x)           free(work->x);
        if (work->z)           free(work->z);
        if (work->xz_tilde)    free(work->xz_tilde);
        if (work->x_prev)      free(work->x_prev);
        if (work->z_prev)      free(work->z_prev);
        if (work->y)           free(work->y);
        if (work->Ax)          free(work->Ax);
        if (work->Px)          free(work->Px);
        if (work->Aty)         free(work->Aty);
        if (work->delta_y)     free(work->delta_y);
        if (work->Atdelta_y)   free(work->Atdelta_y);
        if (work->delta_x)     free(work->delta_x);
        if (work->Pdelta_x)    free(work->Pdelta_x);
        if (work->Adelta_x)    free(work->Adelta_x);

        // Free Settings
        if (work->params) free(work->params);

        // Free solution
        if (work->solution) {
            if (work->solution->x) free(work->solution->x);
            if (work->solution->y) free(work->solution->y);
            free(work->solution);
        }

        // Free information
        if (work->info) free(work->info);


        // Free timer
        if (work->timer) free(work->timer);


        // Free work
        free(work);
    }

    return exitflag;
}




/************************
* Update problem data  *
************************/
int qp_update_lin_cost(qpWorkspace *work, const float *q_new) {

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);


    if (work->clear_update_time == 1) {
        work->clear_update_time = 0;
        work->info->update_time = 0.0;
    }
    qp_tic(work->timer); // Start timer


    // Replace q by the new vector
    prea_vec_copy(q_new, work->data->q, work->data->n);

    // Scaling
    if (work->params->scaling) {
        vec_ew_prod(work->scaling->D, work->data->q, work->data->q, work->data->n);
        vec_mult_scalar(work->data->q, work->scaling->c, work->data->n);
    }

    // Reset solver information
    reset_info(work->info);


    work->info->update_time += qp_toc(work->timer);


    return 0;
}

int qp_update_bounds(qpWorkspace *work,
                     const float *l_new,
                     const float *u_new) {
    int i, exitflag = 0;

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);


    if (work->clear_update_time == 1) {
        work->clear_update_time = 0;
        work->info->update_time = 0.0;
    }
    qp_tic(work->timer); // Start timer


    // Check if lower bound is smaller than upper bound
    for (i = 0; i < work->data->m; i++) {
        if (l_new[i] > u_new[i]) {

            printf("lower bound must be lower than or equal to upper bound");

            return 1;
        }
    }

    // Replace l and u by the new vectors
    prea_vec_copy(l_new, work->data->l, work->data->m);
    prea_vec_copy(u_new, work->data->u, work->data->m);

    // Scaling
    if (work->params->scaling) {
        vec_ew_prod(work->scaling->E, work->data->l, work->data->l, work->data->m);
        vec_ew_prod(work->scaling->E, work->data->u, work->data->u, work->data->m);
    }

    // Reset solver information
    reset_info(work->info);


    // Update rho_vec and refactor if constraints type changes
    exitflag = update_rho_vec(work);



    work->info->update_time += qp_toc(work->timer);


    return exitflag;
}

int qp_update_lower_bound(qpWorkspace *work, const float *l_new) {
    int i, exitflag = 0;

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);


    if (work->clear_update_time == 1) {
        work->clear_update_time = 0;
        work->info->update_time = 0.0;
    }
    qp_tic(work->timer); // Start timer


    // Replace l by the new vector
    prea_vec_copy(l_new, work->data->l, work->data->m);

    // Scaling
    if (work->params->scaling) {
        vec_ew_prod(work->scaling->E, work->data->l, work->data->l, work->data->m);
    }

    // Check if lower bound is smaller than upper bound
    for (i = 0; i < work->data->m; i++) {
        if (work->data->l[i] > work->data->u[i]) {

            printf("upper bound must be greater than or equal to lower bound");

            return 1;
        }
    }

    // Reset solver information
    reset_info(work->info);


    // Update rho_vec and refactor if constraints type changes
    exitflag = update_rho_vec(work);



    work->info->update_time += qp_toc(work->timer);


    return exitflag;
}

int qp_update_upper_bound(qpWorkspace *work, const float *u_new) {
    int i, exitflag = 0;

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);


    if (work->clear_update_time == 1) {
        work->clear_update_time = 0;
        work->info->update_time = 0.0;
    }
    qp_tic(work->timer); // Start timer


    // Replace u by the new vector
    prea_vec_copy(u_new, work->data->u, work->data->m);

    // Scaling
    if (work->params->scaling) {
        vec_ew_prod(work->scaling->E, work->data->u, work->data->u, work->data->m);
    }

    // Check if upper bound is greater than lower bound
    for (i = 0; i < work->data->m; i++) {
        if (work->data->u[i] < work->data->l[i]) {

            printf("lower bound must be lower than or equal to upper bound");

            return 1;
        }
    }

    // Reset solver information
    reset_info(work->info);


    // Update rho_vec and refactor if constraints type changes
    exitflag = update_rho_vec(work);



    work->info->update_time += qp_toc(work->timer);


    return exitflag;
}

int qp_warm_start(qpWorkspace *work, const float *x, const float *y) {

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);

    // Update warm_start setting to true
    if (!work->params->warm_start) work->params->warm_start = 1;

    // Copy primal and dual variables into the iterates
    prea_vec_copy(x, work->x, work->data->n);
    prea_vec_copy(y, work->y, work->data->m);

    // Scale iterates
    if (work->params->scaling) {
        vec_ew_prod(work->scaling->Dinv, work->x, work->x, work->data->n);
        vec_ew_prod(work->scaling->Einv, work->y, work->y, work->data->m);
        vec_mult_scalar(work->y, work->scaling->c, work->data->m);
    }

    // Compute Ax = z and store it in z
    mat_vec(work->data->A, work->x, work->z, 0);

    return 0;
}

int qp_warm_start_x(qpWorkspace *work, const float *x) {

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);

    // Update warm_start setting to true
    if (!work->params->warm_start) work->params->warm_start = 1;

    // Copy primal variable into the iterate x
    prea_vec_copy(x, work->x, work->data->n);

    // Scale iterate
    if (work->params->scaling) {
        vec_ew_prod(work->scaling->Dinv, work->x, work->x, work->data->n);
    }

    // Compute Ax = z and store it in z
    mat_vec(work->data->A, work->x, work->z, 0);

    return 0;
}

int qp_warm_start_y(qpWorkspace *work, const float *y) {

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);

    // Update warm_start setting to true
    if (!work->params->warm_start) work->params->warm_start = 1;

    // Copy primal variable into the iterate y
    prea_vec_copy(y, work->y, work->data->m);

    // Scale iterate
    if (work->params->scaling) {
        vec_ew_prod(work->scaling->Einv, work->y, work->y, work->data->m);
        vec_mult_scalar(work->y, work->scaling->c, work->data->m);
    }

    return 0;
}




int qp_update_P(qpWorkspace *work,
                const float *Px_new,
                const int   *Px_new_idx,
                int          P_new_n) {
    int i;        // For indexing
    int exitflag; // Exit flag
    int nnzP;     // Number of nonzeros in P

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);


    if (work->clear_update_time == 1) {
        work->clear_update_time = 0;
        work->info->update_time = 0.0;
    }
    qp_tic(work->timer); // Start timer


    nnzP = work->data->P->p[work->data->P->n];

    if (Px_new_idx) { // Passing the index of elements changed
        // Check if number of elements is less or equal than the total number of
        // nonzeros in P
        if (P_new_n > nnzP) {

            printf("new number of elements (%i) greater than elements in P (%i)",
                   (int)P_new_n,
                   (int)nnzP);

            return 1;
        }
    }

    if (work->params->scaling) {
        // Unscale data
        unscale_data(work);
    }

    // Update P elements
    if (Px_new_idx) { // Change only Px_new_idx
        for (i = 0; i < P_new_n; i++) {
            work->data->P->x[Px_new_idx[i]] = Px_new[i];
        }
    }
    else // Change whole P
    {
        for (i = 0; i < nnzP; i++) {
            work->data->P->x[i] = Px_new[i];
        }
    }

    if (work->params->scaling) {
        // Scale data
        scale_data(work);
    }

    // Update linear system structure with new data
    exitflag = work->linsys_solver->update_matrices(work->linsys_solver,
               work->data->P,
               work->data->A);

    // Reset solver information
    reset_info(work->info);



    if (exitflag < 0) {
        printf("new KKT matrix is not quasidefinite");
    }



    work->info->update_time += qp_toc(work->timer);


    return exitflag;
}


int qp_update_A(qpWorkspace *work,
                const float *Ax_new,
                const int   *Ax_new_idx,
                int          A_new_n) {
    int i;        // For indexing
    int exitflag; // Exit flag
    int nnzA;     // Number of nonzeros in A

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);


    if (work->clear_update_time == 1) {
        work->clear_update_time = 0;
        work->info->update_time = 0.0;
    }
    qp_tic(work->timer); // Start timer


    nnzA = work->data->A->p[work->data->A->n];

    if (Ax_new_idx) { // Passing the index of elements changed
        // Check if number of elements is less or equal than the total number of
        // nonzeros in A
        if (A_new_n > nnzA) {

            printf("new number of elements (%i) greater than elements in A (%i)",
                   (int)A_new_n,
                   (int)nnzA);

            return 1;
        }
    }

    if (work->params->scaling) {
        // Unscale data
        unscale_data(work);
    }

    // Update A elements
    if (Ax_new_idx) { // Change only Ax_new_idx
        for (i = 0; i < A_new_n; i++) {
            work->data->A->x[Ax_new_idx[i]] = Ax_new[i];
        }
    }
    else { // Change whole A
        for (i = 0; i < nnzA; i++) {
            work->data->A->x[i] = Ax_new[i];
        }
    }

    if (work->params->scaling) {
        // Scale data
        scale_data(work);
    }

    // Update linear system structure with new data
    exitflag = work->linsys_solver->update_matrices(work->linsys_solver,
               work->data->P,
               work->data->A);

    // Reset solver information
    reset_info(work->info);



    if (exitflag < 0) {
        printf("new KKT matrix is not quasidefinite");
    }



    work->info->update_time += qp_toc(work->timer);


    return exitflag;
}


int qp_update_P_A(qpWorkspace *work,
                  const float *Px_new,
                  const int   *Px_new_idx,
                  int          P_new_n,
                  const float *Ax_new,
                  const int   *Ax_new_idx,
                  int          A_new_n) {
    int i;          // For indexing
    int exitflag;   // Exit flag
    int nnzP, nnzA; // Number of nonzeros in P and A

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);


    if (work->clear_update_time == 1) {
        work->clear_update_time = 0;
        work->info->update_time = 0.0;
    }
    qp_tic(work->timer); // Start timer


    nnzP = work->data->P->p[work->data->P->n];
    nnzA = work->data->A->p[work->data->A->n];


    if (Px_new_idx) { // Passing the index of elements changed
        // Check if number of elements is less or equal than the total number of
        // nonzeros in P
        if (P_new_n > nnzP) {

            printf("new number of elements (%i) greater than elements in P (%i)",
                   (int)P_new_n,
                   (int)nnzP);

            return 1;
        }
    }


    if (Ax_new_idx) { // Passing the index of elements changed
        // Check if number of elements is less or equal than the total number of
        // nonzeros in A
        if (A_new_n > nnzA) {

            printf("new number of elements (%i) greater than elements in A (%i)",
                   (int)A_new_n,
                   (int)nnzA);

            return 2;
        }
    }

    if (work->params->scaling) {
        // Unscale data
        unscale_data(work);
    }

    // Update P elements
    if (Px_new_idx) { // Change only Px_new_idx
        for (i = 0; i < P_new_n; i++) {
            work->data->P->x[Px_new_idx[i]] = Px_new[i];
        }
    }
    else // Change whole P
    {
        for (i = 0; i < nnzP; i++) {
            work->data->P->x[i] = Px_new[i];
        }
    }

    // Update A elements
    if (Ax_new_idx) { // Change only Ax_new_idx
        for (i = 0; i < A_new_n; i++) {
            work->data->A->x[Ax_new_idx[i]] = Ax_new[i];
        }
    }
    else { // Change whole A
        for (i = 0; i < nnzA; i++) {
            work->data->A->x[i] = Ax_new[i];
        }
    }

    if (work->params->scaling) {
        // Scale data
        scale_data(work);
    }

    // Update linear system structure with new data
    exitflag = work->linsys_solver->update_matrices(work->linsys_solver,
               work->data->P,
               work->data->A);

    // Reset solver information
    reset_info(work->info);



    if (exitflag < 0) {
        printf("new KKT matrix is not quasidefinite");
    }



    work->info->update_time += qp_toc(work->timer);


    return exitflag;
}

int qp_update_rho(qpWorkspace *work, float rho_new) {
    int exitflag, i;

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);

    // Check value of rho
    if (rho_new <= 0) {

        printf("rho must be positive");

        return 1;
    }


    if (work->rho_update_from_solve == 0) {
        if (work->clear_update_time == 1) {
            work->clear_update_time = 0;
            work->info->update_time = 0.0;
        }
        qp_tic(work->timer); // Start timer
    }


    // Update rho in params
    work->params->rho = c_min(c_max(rho_new, RHO_MIN), RHO_MAX);

    // Update rho_vec and rho_inv_vec
    for (i = 0; i < work->data->m; i++) {
        if (work->constr_type[i] == 0) {
            // Inequalities
            work->rho_vec[i]     = work->params->rho;
            work->rho_inv_vec[i] = 1. / work->params->rho;
        }
        else if (work->constr_type[i] == 1) {
            // Equalities
            work->rho_vec[i]     = RHO_EQ_OVER_RHO_INEQ * work->params->rho;
            work->rho_inv_vec[i] = 1. / work->rho_vec[i];
        }
    }

    // Update rho_vec in KKT matrix
    exitflag = work->linsys_solver->update_rho_vec(work->linsys_solver,
               work->rho_vec);


    if (work->rho_update_from_solve == 0)
        work->info->update_time += qp_toc(work->timer);


    return exitflag;
}

/****************************
* Update problem params  *
****************************/
int qp_update_max_iter(qpWorkspace *work, int max_iter_new) {

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);

    // Check that max_iter is positive
    if (max_iter_new <= 0) {

        printf("max_iter must be positive");

        return 1;
    }

    // Update max_iter
    work->params->max_iter = max_iter_new;

    return 0;
}

int qp_update_eps_abs(qpWorkspace *work, float eps_abs_new) {

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);

    // Check that eps_abs is positive
    if (eps_abs_new < 0.) {

        printf("eps_abs must be nonnegative");

        return 1;
    }

    // Update eps_abs
    work->params->eps_abs = eps_abs_new;

    return 0;
}

int qp_update_eps_rel(qpWorkspace *work, float eps_rel_new) {

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);

    // Check that eps_rel is positive
    if (eps_rel_new < 0.) {

        printf("eps_rel must be nonnegative");

        return 1;
    }

    // Update eps_rel
    work->params->eps_rel = eps_rel_new;

    return 0;
}

int qp_update_eps_prim_inf(qpWorkspace *work, float eps_prim_inf_new) {

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);

    // Check that eps_prim_inf is positive
    if (eps_prim_inf_new < 0.) {

        printf("eps_prim_inf must be nonnegative");

        return 1;
    }

    // Update eps_prim_inf
    work->params->eps_prim_inf = eps_prim_inf_new;

    return 0;
}

int qp_update_eps_dual_inf(qpWorkspace *work, float eps_dual_inf_new) {

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);

    // Check that eps_dual_inf is positive
    if (eps_dual_inf_new < 0.) {

        printf("eps_dual_inf must be nonnegative");

        return 1;
    }

    // Update eps_dual_inf
    work->params->eps_dual_inf = eps_dual_inf_new;


    return 0;
}

int qp_update_alpha(qpWorkspace *work, float alpha_new) {

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);

    // Check that alpha is between 0 and 2
    if ((alpha_new <= 0.) || (alpha_new >= 2.)) {

        printf("alpha must be between 0 and 2");

        return 1;
    }

    // Update alpha
    work->params->alpha = alpha_new;

    return 0;
}

int qp_update_warm_start(qpWorkspace *work, int warm_start_new) {

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);

    // Check that warm_start is either 0 or 1
    if ((warm_start_new != 0) && (warm_start_new != 1)) {

        printf("warm_start should be either 0 or 1");

        return 1;
    }

    // Update warm_start
    work->params->warm_start = warm_start_new;

    return 0;
}

int qp_update_scaled_termination(qpWorkspace *work, int scaled_termination_new) {

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);

    // Check that scaled_termination is either 0 or 1
    if ((scaled_termination_new != 0) && (scaled_termination_new != 1)) {

        printf("scaled_termination should be either 0 or 1");

        return 1;
    }

    // Update scaled_termination
    work->params->scaled_termination = scaled_termination_new;

    return 0;
}

int qp_update_check_termination(qpWorkspace *work, int check_termination_new) {

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);

    // Check that check_termination is nonnegative
    if (check_termination_new < 0) {

        printf("check_termination should be nonnegative");

        return 1;
    }

    // Update check_termination
    work->params->check_termination = check_termination_new;

    return 0;
}



int qp_update_delta(qpWorkspace *work, float delta_new) {

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);

    // Check that delta is positive
    if (delta_new <= 0.) {

        printf("delta must be positive");

        return 1;
    }

    // Update delta
    work->params->delta = delta_new;

    return 0;
}

int qp_update_polish(qpWorkspace *work, int polish_new) {

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);

    // Check that polish is either 0 or 1
    if ((polish_new != 0) && (polish_new != 1)) {

        printf("polish should be either 0 or 1");

        return 1;
    }

    // Update polish
    work->params->polish = polish_new;



    // Reset polish time to zero
    work->info->polish_time = 0.0;


    return 0;
}

int qp_update_polish_refine_iter(qpWorkspace *work, int polish_refine_iter_new) {

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);

    // Check that polish_refine_iter is nonnegative
    if (polish_refine_iter_new < 0) {

        printf("polish_refine_iter must be nonnegative");

        return 1;
    }

    // Update polish_refine_iter
    work->params->polish_refine_iter = polish_refine_iter_new;

    return 0;
}

int qp_update_verbose(qpWorkspace *work, int verbose_new) {

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);

    // Check that verbose is either 0 or 1
    if ((verbose_new != 0) && (verbose_new != 1)) {

        printf("verbose should be either 0 or 1");

        return 1;
    }

    // Update verbose
    work->params->verbose = verbose_new;

    return 0;
}




int qp_update_time_limit(qpWorkspace *work, float time_limit_new) {

    // Check if workspace has been initialized
    if (!work) return qp_error(QP_WORKSPACE_NOT_INIT_ERROR);

    // Check that time_limit is nonnegative
    if (time_limit_new < 0.) {

        printf("time_limit must be nonnegative\n");

        return 1;
    }

    // Update time_limit
    work->params->time_limit = time_limit_new;

    return 0;
}

