#include "../include/qp.h" // For qp rho update
#include "../include/auxil.h"
// #include "../include/timer.h"
// #include "proj.h"
// #include "lin_alg.h"
// #include "constants.h"
// #include "scaling.h"
// #include "util.h"

/***********************************************************
* Auxiliary functions needed to compute ADMM iterations * *
***********************************************************/
float compute_rho_estimate(qpWorkspace *work) {
    int   n, m;                       // Dimensions
    float pri_res, dua_res;           // Primal and dual residuals
    float pri_res_norm, dua_res_norm; // Normalization for the residuals
    float temp_res_norm;              // Temporary residual norm
    float rho_estimate;               // Rho estimate value

    // Get problem dimensions
    n = work->data->n;
    m = work->data->m;

    // Get primal and dual residuals
    pri_res = vec_norm_inf(work->z_prev, m);
    dua_res = vec_norm_inf(work->x_prev, n);

    // Normalize primal residual
    pri_res_norm  = vec_norm_inf(work->z, m);           // ||z||
    temp_res_norm = vec_norm_inf(work->Ax, m);          // ||Ax||
    pri_res_norm  = c_max(pri_res_norm, temp_res_norm); // max (||z||,||Ax||)
    pri_res      /= (pri_res_norm + 1e-10);             // Normalize primal
    // residual (prevent 0
    // division)

    // Normalize dual residual
    dua_res_norm  = vec_norm_inf(work->data->q, n);     // ||q||
    temp_res_norm = vec_norm_inf(work->Aty, n);         // ||A' y||
    dua_res_norm  = c_max(dua_res_norm, temp_res_norm);
    temp_res_norm = vec_norm_inf(work->Px, n);          //  ||P x||
    dua_res_norm  = c_max(dua_res_norm, temp_res_norm); // max(||q||,||A' y||,||P
    // x||)
    dua_res      /= (dua_res_norm + 1e-10);             // Normalize dual residual
    // (prevent 0 division)


    // Return rho estimate
    rho_estimate = work->params->rho * c_sqrt(pri_res / (dua_res + 1e-10)); // (prevent
    // 0
    // division)
    rho_estimate = c_min(c_max(rho_estimate, RHO_MIN), RHO_MAX);              // Constrain
    // rho
    // values
    return rho_estimate;
}

int adapt_rho(qpWorkspace *work) {
    int   exitflag; // Exitflag
    float rho_new;  // New rho value

    exitflag = 0;     // Initialize exitflag to 0

    // Compute new rho
    rho_new = compute_rho_estimate(work);

    // Set rho estimate in info
    work->info->rho_estimate = rho_new;

    // Check if the new rho is large or small enough and update it in case
    if ((rho_new > work->params->rho * work->params->adaptive_rho_tolerance) ||
            (rho_new < work->params->rho /  work->params->adaptive_rho_tolerance)) {
        exitflag                 = qp_update_rho(work, rho_new);
        work->info->rho_updates += 1;
    }

    return exitflag;
}

void set_rho_vec(qpWorkspace *work) {
    int i;

    work->params->rho = c_min(c_max(work->params->rho, RHO_MIN), RHO_MAX);

    for (i = 0; i < work->data->m; i++) {
        if ((work->data->l[i] < -qpINFTY * MIN_SCALING) &&
                (work->data->u[i] > qpINFTY * MIN_SCALING)) {
            // Loose bounds
            work->constr_type[i] = -1;
            work->rho_vec[i]     = RHO_MIN;
        } else if (work->data->u[i] - work->data->l[i] < RHO_TOL) {
            // Equality constraints
            work->constr_type[i] = 1;
            work->rho_vec[i]     = RHO_EQ_OVER_RHO_INEQ * work->params->rho;
        } else {
            // Inequality constraints
            work->constr_type[i] = 0;
            work->rho_vec[i]     = work->params->rho;
        }
        work->rho_inv_vec[i] = 1. / work->rho_vec[i];
    }
}

int update_rho_vec(qpWorkspace *work) {
    int i, exitflag, constr_type_changed;

    exitflag            = 0;
    constr_type_changed = 0;

    for (i = 0; i < work->data->m; i++) {
        if ((work->data->l[i] < -qpINFTY * MIN_SCALING) &&
                (work->data->u[i] > qpINFTY * MIN_SCALING)) {
            // Loose bounds
            if (work->constr_type[i] != -1) {
                work->constr_type[i] = -1;
                work->rho_vec[i]     = RHO_MIN;
                work->rho_inv_vec[i] = 1. / RHO_MIN;
                constr_type_changed  = 1;
            }
        } else if (work->data->u[i] - work->data->l[i] < RHO_TOL) {
            // Equality constraints
            if (work->constr_type[i] != 1) {
                work->constr_type[i] = 1;
                work->rho_vec[i]     = RHO_EQ_OVER_RHO_INEQ * work->params->rho;
                work->rho_inv_vec[i] = 1. / work->rho_vec[i];
                constr_type_changed  = 1;
            }
        } else {
            // Inequality constraints
            if (work->constr_type[i] != 0) {
                work->constr_type[i] = 0;
                work->rho_vec[i]     = work->params->rho;
                work->rho_inv_vec[i] = 1. / work->params->rho;
                constr_type_changed  = 1;
            }
        }
    }

    // Update rho_vec in KKT matrix if constraints type has changed
    if (constr_type_changed == 1) {
        exitflag = work->linsys_solver->update_rho_vec(work->linsys_solver,
                   work->rho_vec);
    }

    return exitflag;
}

void swap_vectors(float **a, float **b) {
    float *temp;

    temp = *b;
    *b   = *a;
    *a   = temp;
}

void cold_start(qpWorkspace *work) {
    vec_set_scalar(work->x, 0., work->data->n);
    vec_set_scalar(work->z, 0., work->data->m);
    vec_set_scalar(work->y, 0., work->data->m);
}

static void compute_rhs(qpWorkspace *work) {
    int i; // Index

    for (i = 0; i < work->data->n; i++) {
        // Cycle over part related to x variables
        work->xz_tilde[i] = work->params->sigma * work->x_prev[i] -
                            work->data->q[i];
    }

    for (i = 0; i < work->data->m; i++) {
        // Cycle over dual variable in the first step (nu)
        work->xz_tilde[i + work->data->n] = work->z_prev[i] - work->rho_inv_vec[i] *
                                            work->y[i];
    }
}

void update_xz_tilde(qpWorkspace *work) {
    // Compute right-hand side
    compute_rhs(work);

    // Solve linear system
    work->linsys_solver->solve(work->linsys_solver, work->xz_tilde);
}

void update_x(qpWorkspace *work) {
    int i;

    // update x
    for (i = 0; i < work->data->n; i++) {
        work->x[i] = work->params->alpha * work->xz_tilde[i] +
                     ((float)1.0 - work->params->alpha) * work->x_prev[i];
    }

    // update delta_x
    for (i = 0; i < work->data->n; i++) {
        work->delta_x[i] = work->x[i] - work->x_prev[i];
    }
}

void update_z(qpWorkspace *work) {
    int i;

    // update z
    for (i = 0; i < work->data->m; i++) {
        work->z[i] = work->params->alpha * work->xz_tilde[i + work->data->n] +
                     ((float)1.0 - work->params->alpha) * work->z_prev[i] +
                     work->rho_inv_vec[i] * work->y[i];
    }

    // project z
    project(work, work->z);
}

void update_y(qpWorkspace *work) {
    int i; // Index

    for (i = 0; i < work->data->m; i++) {
        work->delta_y[i] = work->rho_vec[i] *
                           (work->params->alpha *
                            work->xz_tilde[i + work->data->n] +
                            ((float)1.0 - work->params->alpha) * work->z_prev[i] -
                            work->z[i]);
        work->y[i] += work->delta_y[i];
    }
}

float compute_obj_val(qpWorkspace *work, float *x) {
    float obj_val;

    obj_val = quad_form(work->data->P, x) +
              vec_prod(work->data->q, x, work->data->n);

    if (work->params->scaling) {
        obj_val *= work->scaling->cinv;
    }

    return obj_val;
}

float compute_pri_res(qpWorkspace *work, float *x, float *z) {
    // NB: Use z_prev as working vector
    // pr = Ax - z

    mat_vec(work->data->A, x, work->Ax, 0); // Ax
    vec_add_scaled(work->z_prev, work->Ax, z, work->data->m, -1);

    // If scaling active -> rescale residual
    if (work->params->scaling && !work->params->scaled_termination) {
        return vec_scaled_norm_inf(work->scaling->Einv, work->z_prev, work->data->m);
    }

    // Return norm of the residual
    return vec_norm_inf(work->z_prev, work->data->m);
}

float compute_pri_tol(qpWorkspace *work, float eps_abs, float eps_rel) {
    float max_rel_eps, temp_rel_eps;

    // max_rel_eps = max(||z||, ||A x||)
    if (work->params->scaling && !work->params->scaled_termination) {
        // ||Einv * z||
        max_rel_eps =
            vec_scaled_norm_inf(work->scaling->Einv, work->z, work->data->m);

        // ||Einv * A * x||
        temp_rel_eps = vec_scaled_norm_inf(work->scaling->Einv,
                                           work->Ax,
                                           work->data->m);

        // Choose maximum
        max_rel_eps = c_max(max_rel_eps, temp_rel_eps);
    } else { // No unscaling required
        // ||z||
        max_rel_eps = vec_norm_inf(work->z, work->data->m);

        // ||A * x||
        temp_rel_eps = vec_norm_inf(work->Ax, work->data->m);

        // Choose maximum
        max_rel_eps = c_max(max_rel_eps, temp_rel_eps);
    }

    // eps_prim
    return eps_abs + eps_rel * max_rel_eps;
}

float compute_dua_res(qpWorkspace *work, float *x, float *y) {
    // NB: Use x_prev as temporary vector
    // NB: Only upper triangular part of P is stored.
    // dr = q + A'*y + P*x

    // dr = q
    prea_vec_copy(work->data->q, work->x_prev, work->data->n);

    // P * x (upper triangular part)
    mat_vec(work->data->P, x, work->Px, 0);

    // P' * x (lower triangular part with no diagonal)
    mat_tpose_vec(work->data->P, x, work->Px, 1, 1);

    // dr += P * x (full P matrix)
    vec_add_scaled(work->x_prev, work->x_prev, work->Px, work->data->n, 1);

    // dr += A' * y
    if (work->data->m > 0) {
        mat_tpose_vec(work->data->A, y, work->Aty, 0, 0);
        vec_add_scaled(work->x_prev, work->x_prev, work->Aty, work->data->n, 1);
    }

    // If scaling active -> rescale residual
    if (work->params->scaling && !work->params->scaled_termination) {
        return work->scaling->cinv * vec_scaled_norm_inf(work->scaling->Dinv,
                work->x_prev,
                work->data->n);
    }

    return vec_norm_inf(work->x_prev, work->data->n);
}

float compute_dua_tol(qpWorkspace *work, float eps_abs, float eps_rel) {
    float max_rel_eps, temp_rel_eps;

    // max_rel_eps = max(||q||, ||A' y|, ||P x||)
    if (work->params->scaling && !work->params->scaled_termination) {
        // || Dinv q||
        max_rel_eps = vec_scaled_norm_inf(work->scaling->Dinv,
                                          work->data->q,
                                          work->data->n);

        // || Dinv A' y ||
        temp_rel_eps = vec_scaled_norm_inf(work->scaling->Dinv,
                                           work->Aty,
                                           work->data->n);
        max_rel_eps = c_max(max_rel_eps, temp_rel_eps);

        // || Dinv P x||
        temp_rel_eps = vec_scaled_norm_inf(work->scaling->Dinv,
                                           work->Px,
                                           work->data->n);
        max_rel_eps = c_max(max_rel_eps, temp_rel_eps);

        // Multiply by cinv
        max_rel_eps *= work->scaling->cinv;
    } else { // No scaling required
        // ||q||
        max_rel_eps = vec_norm_inf(work->data->q, work->data->n);

        // ||A'*y||
        temp_rel_eps = vec_norm_inf(work->Aty, work->data->n);
        max_rel_eps  = c_max(max_rel_eps, temp_rel_eps);

        // ||P*x||
        temp_rel_eps = vec_norm_inf(work->Px, work->data->n);
        max_rel_eps  = c_max(max_rel_eps, temp_rel_eps);
    }

    // eps_dual
    return eps_abs + eps_rel * max_rel_eps;
}

int is_primal_infeasible(qpWorkspace *work, float eps_prim_inf) {
    // This function checks for the primal infeasibility termination criteria.
    //
    // 1) A' * delta_y < eps * ||delta_y||
    //
    // 2) u'*max(delta_y, 0) + l'*min(delta_y, 0) < -eps * ||delta_y||
    //

    int i; // Index for loops
    float norm_delta_y;
    float ineq_lhs = 0.0;

    // Project delta_y onto the polar of the recession cone of [l,u]
    for (i = 0; i < work->data->m; i++) {
        if (work->data->u[i] > qpINFTY * MIN_SCALING) {          // Infinite upper bound
            if (work->data->l[i] < -qpINFTY * MIN_SCALING) {       // Infinite lower bound
                // Both bounds infinite
                work->delta_y[i] = 0.0;
            } else {
                // Only upper bound infinite
                work->delta_y[i] = c_min(work->delta_y[i], 0.0);
            }
        } else if (work->data->l[i] < -qpINFTY * MIN_SCALING) {  // Infinite lower bound
            // Only lower bound infinite
            work->delta_y[i] = c_max(work->delta_y[i], 0.0);
        }
    }

    // Compute infinity norm of delta_y (unscale if necessary)
    if (work->params->scaling && !work->params->scaled_termination) {
        // Use work->Adelta_x as temporary vector
        vec_ew_prod(work->scaling->E, work->delta_y, work->Adelta_x, work->data->m);
        norm_delta_y = vec_norm_inf(work->Adelta_x, work->data->m);
    } else {
        norm_delta_y = vec_norm_inf(work->delta_y, work->data->m);
    }

    if (norm_delta_y > eps_prim_inf) { // ||delta_y|| > 0

        for (i = 0; i < work->data->m; i++) {
            ineq_lhs += work->data->u[i] * c_max(work->delta_y[i], 0) + \
                        work->data->l[i] * c_min(work->delta_y[i], 0);
        }

        // Check if the condition is satisfied: ineq_lhs < -eps
        if (ineq_lhs < -eps_prim_inf * norm_delta_y) {
            // Compute and return ||A'delta_y|| < eps_prim_inf
            mat_tpose_vec(work->data->A, work->delta_y, work->Atdelta_y, 0, 0);

            // Unscale if necessary
            if (work->params->scaling && !work->params->scaled_termination) {
                vec_ew_prod(work->scaling->Dinv,
                            work->Atdelta_y,
                            work->Atdelta_y,
                            work->data->n);
            }

            return vec_norm_inf(work->Atdelta_y, work->data->n) < eps_prim_inf * norm_delta_y;
        }
    }

    // Conditions not satisfied -> not primal infeasible
    return 0;
}

int is_dual_infeasible(qpWorkspace *work, float eps_dual_inf) {
    // This function checks for the scaled dual infeasibility termination
    // criteria.
    //
    // 1) q * delta_x < - eps * || delta_x ||
    //
    // 2) ||P * delta_x || < eps * || delta_x ||
    //
    // 3) -> (A * delta_x)_i > -eps * || delta_x ||,    l_i != -inf
    //    -> (A * delta_x)_i <  eps * || delta_x ||,    u_i != inf
    //


    int   i; // Index for loops
    float norm_delta_x;
    float cost_scaling;

    // Compute norm of delta_x
    if (work->params->scaling && !work->params->scaled_termination) { // Unscale
        // if
        // necessary
        norm_delta_x = vec_scaled_norm_inf(work->scaling->D,
                                           work->delta_x,
                                           work->data->n);
        cost_scaling = work->scaling->c;
    } else {
        norm_delta_x = vec_norm_inf(work->delta_x, work->data->n);
        cost_scaling = 1.0;
    }

    // Prevent 0 division || delta_x || > 0
    if (norm_delta_x > eps_dual_inf) {
        // Normalize delta_x by its norm

        /* vec_mult_scalar(work->delta_x, 1./norm_delta_x, work->data->n); */

        // Check first if q'*delta_x < 0
        if (vec_prod(work->data->q, work->delta_x, work->data->n) <
                -cost_scaling * eps_dual_inf * norm_delta_x) {
            // Compute product P * delta_x (NB: P is store in upper triangular form)
            mat_vec(work->data->P, work->delta_x, work->Pdelta_x, 0);
            mat_tpose_vec(work->data->P, work->delta_x, work->Pdelta_x, 1, 1);

            // Scale if necessary
            if (work->params->scaling && !work->params->scaled_termination) {
                vec_ew_prod(work->scaling->Dinv,
                            work->Pdelta_x,
                            work->Pdelta_x,
                            work->data->n);
            }

            // Check if || P * delta_x || = 0
            if (vec_norm_inf(work->Pdelta_x, work->data->n) <
                    cost_scaling * eps_dual_inf * norm_delta_x) {
                // Compute A * delta_x
                mat_vec(work->data->A, work->delta_x, work->Adelta_x, 0);

                // Scale if necessary
                if (work->params->scaling && !work->params->scaled_termination) {
                    vec_ew_prod(work->scaling->Einv,
                                work->Adelta_x,
                                work->Adelta_x,
                                work->data->m);
                }

                // De Morgan Law Applied to dual infeasibility conditions for A * x
                // NB: Note that MIN_SCALING is used to adjust the infinity value
                //     in case the problem is scaled.
                for (i = 0; i < work->data->m; i++) {
                    if (((work->data->u[i] < qpINFTY * MIN_SCALING) &&
                            (work->Adelta_x[i] >  eps_dual_inf * norm_delta_x)) ||
                            ((work->data->l[i] > -qpINFTY * MIN_SCALING) &&
                             (work->Adelta_x[i] < -eps_dual_inf * norm_delta_x))) {
                        // At least one condition not satisfied -> not dual infeasible
                        return 0;
                    }
                }

                // All conditions passed -> dual infeasible
                return 1;
            }
        }
    }

    // Conditions not satisfied -> not dual infeasible
    return 0;
}

int has_solution(qpInfo * info) {

    return ((info->status_val != qp_PRIMAL_INFEASIBLE) &&
            (info->status_val != qp_PRIMAL_INFEASIBLE_INACCURATE) &&
            (info->status_val != qp_DUAL_INFEASIBLE) &&
            (info->status_val != qp_DUAL_INFEASIBLE_INACCURATE) &&
            (info->status_val != qp_NON_CVX));

}

void store_solution(qpWorkspace *work) {
    float norm_vec;

    if (has_solution(work->info)) {
        prea_vec_copy(work->x, work->solution->x, work->data->n); // primal
        prea_vec_copy(work->y, work->solution->y, work->data->m); // dual

        // Unscale solution if scaling has been performed
        if (work->params->scaling)
            unscale_solution(work);
    } else {
        // No solution present. Solution is NaN
        vec_set_scalar(work->solution->x, qpNAN, work->data->n);
        vec_set_scalar(work->solution->y, qpNAN, work->data->m);


        // Normalize infeasibility certificates if embedded is off
        // NB: It requires a division
        if ((work->info->status_val == qp_PRIMAL_INFEASIBLE) ||
                ((work->info->status_val == qp_PRIMAL_INFEASIBLE_INACCURATE))) {
            norm_vec = vec_norm_inf(work->delta_y, work->data->m);
            vec_mult_scalar(work->delta_y, 1. / norm_vec, work->data->m);
        }

        if ((work->info->status_val == qp_DUAL_INFEASIBLE) ||
                ((work->info->status_val == qp_DUAL_INFEASIBLE_INACCURATE))) {
            norm_vec = vec_norm_inf(work->delta_x, work->data->n);
            vec_mult_scalar(work->delta_x, 1. / norm_vec, work->data->n);
        }


        // Cold start iterates to 0 for next runs (they cannot start from NaN)
        cold_start(work);
    }
}

void update_info(qpWorkspace *work,
                 int          iter,
                 int          compute_objective,
                 int          polish) {
    float *x, *z, *y;                   // Allocate pointers to variables
    float *obj_val, *pri_res, *dua_res; // objective value, residuals

    float *run_time;                    // Execution time


    if (polish) {
        x       = work->pol->x;
        y       = work->pol->y;
        z       = work->pol->z;
        obj_val = &work->pol->obj_val;
        pri_res = &work->pol->pri_res;
        dua_res = &work->pol->dua_res;
        run_time = &work->info->polish_time;
    } else {
        x                = work->x;
        y                = work->y;
        z                = work->z;
        obj_val          = &work->info->obj_val;
        pri_res          = &work->info->pri_res;
        dua_res          = &work->info->dua_res;
        work->info->iter = iter; // Update iteration number
        run_time = &work->info->solve_time;

    }



    // Compute the objective if needed
    if (compute_objective) {
        *obj_val = compute_obj_val(work, x);
    }

    // Compute primal residual
    if (work->data->m == 0) {
        // No constraints -> Always primal feasible
        *pri_res = 0.;
    } else {
        *pri_res = compute_pri_res(work, x, z);
    }

    // Compute dual residual
    *dua_res = compute_dua_res(work, x, y);

    // Update timing
    *run_time = qp_toc(work->timer);

    work->summary_printed = 0; // The just updated info have not been printed
}


void reset_info(qpInfo *info) {

    // Initialize info values.
    info->solve_time = 0.0;  // Solve time to zero
    info->polish_time = 0.0; // Polish time to zero


    update_status(info, qp_UNSOLVED); // Problem is unsolved

    info->rho_updates = 0;              // Rho updates are now 0
}

void update_status(qpInfo *info, int status_val) {
    // Update status value
    info->status_val = status_val;

    // Update status string depending on status val
    if (status_val == qp_SOLVED) strcpy(info->status, "solved");

    if (status_val == qp_SOLVED_INACCURATE) strcpy(info->status,
                "solved inaccurate");
    else if (status_val == qp_PRIMAL_INFEASIBLE) strcpy(info->status,
                "primal infeasible");
    else if (status_val == qp_PRIMAL_INFEASIBLE_INACCURATE) strcpy(info->status,
                "primal infeasible inaccurate");
    else if (status_val == qp_UNSOLVED) strcpy(info->status, "unsolved");
    else if (status_val == qp_DUAL_INFEASIBLE) strcpy(info->status,
                "dual infeasible");
    else if (status_val == qp_DUAL_INFEASIBLE_INACCURATE) strcpy(info->status,
                "dual infeasible inaccurate");
    else if (status_val == qp_MAX_ITER_REACHED) strcpy(info->status,
                "maximum iterations reached");
    else if (status_val == qp_TIME_LIMIT_REACHED) strcpy(info->status,
                "run time limit reached");
    /* ifdef PROFILING */
    else if (status_val == qp_SIGINT) strcpy(info->status, "interrupted");

    else if (status_val == qp_NON_CVX) strcpy(info->status, "problem non convex");

}

int check_termination(qpWorkspace *work, int approximate) {
    float eps_prim, eps_dual, eps_prim_inf, eps_dual_inf;
    int   exitflag;
    int   prim_res_check, dual_res_check, prim_inf_check, dual_inf_check;
    float eps_abs, eps_rel;

    // Initialize variables to 0
    exitflag       = 0;
    prim_res_check = 0;
    dual_res_check = 0;
    prim_inf_check = 0;
    dual_inf_check = 0;

    // Initialize tolerances
    eps_abs      = work->params->eps_abs;
    eps_rel      = work->params->eps_rel;
    eps_prim_inf = work->params->eps_prim_inf;
    eps_dual_inf = work->params->eps_dual_inf;

    // If residuals are too large, the problem is probably non convex
    if ((work->info->pri_res > qpINFTY) ||
            (work->info->dua_res > qpINFTY)) {
        // Looks like residuals are diverging. Probably the problem is non convex!
        // Terminate and report it
        update_status(work->info, qp_NON_CVX);
        work->info->obj_val = qpNAN;
        return 1;
    }

    // If approximate solution required, increase tolerances by 10
    if (approximate) {
        eps_abs      *= 10;
        eps_rel      *= 10;
        eps_prim_inf *= 10;
        eps_dual_inf *= 10;
    }

    // Check residuals
    if (work->data->m == 0) {
        prim_res_check = 1; // No constraints -> Primal feasibility always satisfied
    }
    else {
        // Compute primal tolerance
        eps_prim = compute_pri_tol(work, eps_abs, eps_rel);

        // Primal feasibility check
        if (work->info->pri_res < eps_prim) {
            prim_res_check = 1;
        } else {
            // Primal infeasibility check
            prim_inf_check = is_primal_infeasible(work, eps_prim_inf);
        }
    } // End check if m == 0

    // Compute dual tolerance
    eps_dual = compute_dua_tol(work, eps_abs, eps_rel);

    // Dual feasibility check
    if (work->info->dua_res < eps_dual) {
        dual_res_check = 1;
    } else {
        // Check dual infeasibility
        dual_inf_check = is_dual_infeasible(work, eps_dual_inf);
    }

    // Compare checks to determine solver status
    if (prim_res_check && dual_res_check) {
        // Update final information
        if (approximate) {
            update_status(work->info, qp_SOLVED_INACCURATE);
        } else {
            update_status(work->info, qp_SOLVED);
        }
        exitflag = 1;
    }
    else if (prim_inf_check) {
        // Update final information
        if (approximate) {
            update_status(work->info, qp_PRIMAL_INFEASIBLE_INACCURATE);
        } else {
            update_status(work->info, qp_PRIMAL_INFEASIBLE);
        }

        if (work->params->scaling && !work->params->scaled_termination) {
            // Update infeasibility certificate
            vec_ew_prod(work->scaling->E, work->delta_y, work->delta_y, work->data->m);
        }
        work->info->obj_val = qpINFTY;
        exitflag            = 1;
    }
    else if (dual_inf_check) {
        // Update final information
        if (approximate) {
            update_status(work->info, qp_DUAL_INFEASIBLE_INACCURATE);
        } else {
            update_status(work->info, qp_DUAL_INFEASIBLE);
        }

        if (work->params->scaling && !work->params->scaled_termination) {
            // Update infeasibility certificate
            vec_ew_prod(work->scaling->D, work->delta_x, work->delta_x, work->data->n);
        }
        work->info->obj_val = -qpINFTY;
        exitflag            = 1;
    }

    return exitflag;
}


int validate_data(const qpData *data) {
    int j, ptr;

    if (!data) {
        printf("Missing data");
        return 1;
    }

    if (!(data->P)) {
        printf("Missing matrix P");
        return 1;
    }

    if (!(data->A)) {
        printf("Missing matrix A");
        return 1;
    }

    // General dimensions Tests
    if ((data->n <= 0) || (data->m < 0)) {
        printf("n must be positive and m nonnegative; n = %i, m = %i",
               (int)data->n, (int)data->m);
        return 1;
    }

    // Matrix P
    if (data->P->m != data->n) {
        printf("P does not have dimension n x n with n = %i", (int)data->n);
        return 1;
    }

    if (data->P->m != data->P->n) {
        printf("P is not square");
        return 1;
    }

    for (j = 0; j < data->n; j++) { // COLUMN
        for (ptr = data->P->p[j]; ptr < data->P->p[j + 1]; ptr++) {
            if (data->P->i[ptr] > j) {  // if ROW > COLUMN
                printf("P is not upper triangular");
                return 1;
            }
        }
    }

    // Matrix A
    if ((data->A->m != data->m) || (data->A->n != data->n)) {
        printf("A does not have dimension %i x %i", (int)data->m, (int)data->n);
        return 1;
    }

    // Lower and upper bounds
    for (j = 0; j < data->m; j++) {
        if (data->l[j] > data->u[j]) {
            printf("Lower bound at index %d is greater than upper bound: %.4e > %.4e",
                   (int)j, data->l[j], data->u[j]);
            return 1;
        }
    }

    // TODO: Complete with other checks

    return 0;
}

int validate_linsys_solver(int linsys_solver) {
    if ((linsys_solver != LDL_SOLVER) ) {
        return 1;
    }

    // TODO: Add more solvers in case

    // Valid solver
    return 0;
}

int validate_params(const qpParams *params) {
    if (!params) {
        printf("Missing params!");
        return 1;
    }

    if (params->scaling < 0) {
        printf("scaling must be nonnegative");
        return 1;
    }

    if ((params->adaptive_rho != 0) && (params->adaptive_rho != 1)) {
        printf("adaptive_rho must be either 0 or 1");
        return 1;
    }

    if (params->adaptive_rho_interval < 0) {
        printf("adaptive_rho_interval must be nonnegative");
        return 1;
    }

    if (params->adaptive_rho_fraction <= 0) {
        printf("adaptive_rho_fraction must be positive");
        return 1;
    }

    if (params->adaptive_rho_tolerance < 1.0) {

        printf("adaptive_rho_tolerance must be >= 1");
        /* ifdef PRINTING */
        return 1;
    }

    if (params->polish_refine_iter < 0) {

        printf("polish_refine_iter must be nonnegative");
        /* ifdef PRINTING */
        return 1;
    }

    if (params->rho <= 0.0) {

        printf("rho must be positive");
        /* ifdef PRINTING */
        return 1;
    }

    if (params->sigma <= 0.0) {

        printf("sigma must be positive");
        /* ifdef PRINTING */
        return 1;
    }

    if (params->delta <= 0.0) {

        printf("delta must be positive");
        /* ifdef PRINTING */
        return 1;
    }

    if (params->max_iter <= 0) {

        printf("max_iter must be positive");
        /* ifdef PRINTING */
        return 1;
    }

    if (params->eps_abs < 0.0) {

        printf("eps_abs must be nonnegative");
        /* ifdef PRINTING */
        return 1;
    }

    if (params->eps_rel < 0.0) {

        printf("eps_rel must be nonnegative");
        /* ifdef PRINTING */
        return 1;
    }

    if ((params->eps_rel == 0.0) &&
            (params->eps_abs == 0.0)) {

        printf("at least one of eps_abs and eps_rel must be positive");
        /* ifdef PRINTING */
        return 1;
    }

    if (params->eps_prim_inf <= 0.0) {

        printf("eps_prim_inf must be positive");
        /* ifdef PRINTING */
        return 1;
    }

    if (params->eps_dual_inf <= 0.0) {

        printf("eps_dual_inf must be positive");
        /* ifdef PRINTING */
        return 1;
    }

    if ((params->alpha <= 0.0) ||
            (params->alpha >= 2.0)) {

        printf("alpha must be strictly between 0 and 2");
        /* ifdef PRINTING */
        return 1;
    }

    if (validate_linsys_solver(params->linsys_solver)) {

        printf("linsys_solver not recognized");
        /* ifdef PRINTING */
        return 1;
    }

    if ((params->verbose != 0) &&
            (params->verbose != 1)) {

        printf("verbose must be either 0 or 1");
        /* ifdef PRINTING */
        return 1;
    }

    if ((params->scaled_termination != 0) &&
            (params->scaled_termination != 1)) {

        printf("scaled_termination must be either 0 or 1");
        /* ifdef PRINTING */
        return 1;
    }

    if (params->check_termination < 0) {

        printf("check_termination must be nonnegative");
        /* ifdef PRINTING */
        return 1;
    }

    if ((params->warm_start != 0) &&
            (params->warm_start != 1)) {

        printf("warm_start must be either 0 or 1");
        /* ifdef PRINTING */
        return 1;
    }


    if (params->time_limit < 0.0) {

        printf("time_limit must be nonnegative\n");
        /* ifdef PRINTING */
        return 1;
    }
    /* ifdef PROFILING */

    return 0;
}

// #ifndef EMBEDDED