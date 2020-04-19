/* qp.c */
#include "./include/qp.h"
#include "./include/ops.h"
#include "./include/linalg.h"
#include "./include/kkt.h"
#include "./include/csc.h"
#include "./include/timer.h"
#include "./include/params.h"


void project(qpInstance *problem, float *z)
{
    int ii, m;

    m = problem->data->m;

    for (ii = 0; ii < m; ii++) {
        z[ii] = c_min(c_max(z[ii], problem->data->l[ii]), problem->data->u[ii]);
    }
}

void project_normalcone(qpInstance *problem, float *z, float *y)
{
    int ii, m;

    m = problem->data->m;

    for (ii = 0; ii < m; ii++)
    {
        problem->z_prev[ii] = z[ii] + y[ii];

        z[ii] = c_min(c_max(problem->z_prev[ii], problem->data->l[ii]), problem->data->u[ii]);
        y[ii] = problem->z_prev[ii] - z[ii];
    }
}

void limit_scaling(float *D, int n)
{
    int i;

    for (i = 0; i < n; i++) {
        D[i] = D[i] < MIN_SCALING ? 1.0 : D[i];
        D[i] = D[i] > MAX_SCALING ? MAX_SCALING : D[i];
    }
}

void compute_inf_norm_cols_KKT(const csc *P, const csc *A,
                               float *D, float *D_temp_A,
                               float *E, int n) {
    // First half
    //  [ P ]
    //  [ A ]
    mat_inf_norm_cols_sym_triu(P, D);
    mat_inf_norm_cols(A, D_temp_A);
    vec_ew_max_vec(D, D_temp_A, D, n);

    // Second half
    //  [ A']
    //  [ 0 ]
    mat_inf_norm_rows(A, E);
}

int scale_data(qpInstance *problem) {
    // Scale KKT matrix
    //
    //    [ P   A']
    //    [ A   0 ]
    //
    // with diagonal matrix
    //
    //  S = [ D    ]
    //      [    E ]
    //

    int   i;          // Iterations index
    int   n, m;       // Number of constraints and variables
    float c_temp;     // Cost function scaling
    float inf_norm_q; // Infinity norm of q

    n = problem->data->n;
    m = problem->data->m;

    // Initialize scaling to 1
    problem->scaling->c = 1.0;
    vec_set_scalar(problem->scaling->D,    1., problem->data->n);
    vec_set_scalar(problem->scaling->Dinv, 1., problem->data->n);
    vec_set_scalar(problem->scaling->E,    1., problem->data->m);
    vec_set_scalar(problem->scaling->Einv, 1., problem->data->m);


    for (i = 0; i < problem->settings->scaling; i++) {
        //
        // First Ruiz step
        //

        // Compute norm of KKT columns
        compute_inf_norm_cols_KKT(problem->data->P, problem->data->A,
                                  problem->D_temp, problem->D_temp_A,
                                  problem->E_temp, n);

        // Set to 1 values with 0 norms (avoid crazy scaling)
        limit_scaling(problem->D_temp, n);
        limit_scaling(problem->E_temp, m);

        // Take square root of norms
        vec_ew_sqrt(problem->D_temp, n);
        vec_ew_sqrt(problem->E_temp, m);

        // Divide scalings D and E by themselves
        vec_ew_recipr(problem->D_temp, problem->D_temp, n);
        vec_ew_recipr(problem->E_temp, problem->E_temp, m);

        // Equilibrate matrices P and A and vector q
        // P <- DPD
        mat_premult_diag(problem->data->P, problem->D_temp);
        mat_postmult_diag(problem->data->P, problem->D_temp);

        // A <- EAD
        mat_premult_diag(problem->data->A, problem->E_temp);
        mat_postmult_diag(problem->data->A, problem->D_temp);

        // q <- Dq
        vec_ew_prod(problem->D_temp,     problem->data->q, problem->data->q,    n);

        // Update equilibration matrices D and E
        vec_ew_prod(problem->scaling->D, problem->D_temp,  problem->scaling->D, n);
        vec_ew_prod(problem->scaling->E, problem->E_temp,  problem->scaling->E, m);

        //
        // Cost normalization step
        //

        // Compute avg norm of cols of P
        mat_inf_norm_cols_sym_triu(problem->data->P, problem->D_temp);
        c_temp = vec_mean(problem->D_temp, n);

        // Compute inf norm of q
        inf_norm_q = vec_norm_inf(problem->data->q, n);

        // If norm_q == 0, set it to 1 (ignore it in the scaling)
        // NB: Using the same function as with vectors here
        limit_scaling(&inf_norm_q, 1);

        // Compute max between avg norm of cols of P and inf norm of q
        c_temp = c_max(c_temp, inf_norm_q);

        // Limit scaling (use same function as with vectors)
        limit_scaling(&c_temp, 1);

        // Invert scaling c = 1 / cost_measure
        c_temp = 1. / c_temp;

        // Scale P
        mat_mult_scalar(problem->data->P, c_temp);

        // Scale q
        vec_mult_scalar(problem->data->q, c_temp, n);

        // Update cost scaling
        problem->scaling->c *= c_temp;
    }


    // Store cinv, Dinv, Einv
    problem->scaling->cinv = 1. / problem->scaling->c;
    vec_ew_recipr(problem->scaling->D, problem->scaling->Dinv, problem->data->n);
    vec_ew_recipr(problem->scaling->E, problem->scaling->Einv, problem->data->m);


    // Scale problem vectors l, u
    vec_ew_prod(problem->scaling->E, problem->data->l, problem->data->l, problem->data->m);
    vec_ew_prod(problem->scaling->E, problem->data->u, problem->data->u, problem->data->m);

    return 0;
}

int unscale_data(qpInstance *problem) {
    // Unscale cost
    mat_mult_scalar(problem->data->P, problem->scaling->cinv);
    mat_premult_diag(problem->data->P, problem->scaling->Dinv);
    mat_postmult_diag(problem->data->P, problem->scaling->Dinv);
    vec_mult_scalar(problem->data->q, problem->scaling->cinv, problem->data->n);
    vec_ew_prod(problem->scaling->Dinv, problem->data->q, problem->data->q, problem->data->n);

    // Unscale constraints
    mat_premult_diag(problem->data->A, problem->scaling->Einv);
    mat_postmult_diag(problem->data->A, problem->scaling->Dinv);
    vec_ew_prod(problem->scaling->Einv, problem->data->l, problem->data->l, problem->data->m);
    vec_ew_prod(problem->scaling->Einv, problem->data->u, problem->data->u, problem->data->m);

    return 0;
}

int unscale_solution(qpInstance *problem) {
    // primal
    vec_ew_prod(problem->scaling->D,
                problem->solution->x,
                problem->solution->x,
                problem->data->n);

    // dual
    vec_ew_prod(problem->scaling->E,
                problem->solution->y,
                problem->solution->y,
                problem->data->m);
    vec_mult_scalar(problem->solution->y, problem->scaling->cinv, problem->data->m);

    return 0;
}

static int form_Ared(qpInstance problem) {
    int j, ptr;
    int Ared_nnz = 0;

    // Initialize counters for active constraints
    problem->pol->n_low = 0;
    problem->pol->n_upp = 0;

    /* Guess which linear constraints are lower-active, upper-active and free
     *    A_to_Alow[j] = -1    (if j-th row of A is not inserted in Alow)
     *    A_to_Alow[j] =  i    (if j-th row of A is inserted at i-th row of Alow)
     * Aupp is formed in the equivalent way.
     * Ared is formed by stacking vertically Alow and Aupp.
     */
    for (j = 0; j < problem->data->m; j++) {
        if (problem->z[j] - problem->data->l[j] < -problem->y[j]) { // lower-active
            problem->pol->Alow_to_A[problem->pol->n_low] = j;
            problem->pol->A_to_Alow[j]                = problem->pol->n_low++;
        } else {
            problem->pol->A_to_Alow[j] = -1;
        }
    }

    for (j = 0; j < problem->data->m; j++) {
        if (problem->data->u[j] - problem->z[j] < problem->y[j]) { // upper-active
            problem->pol->Aupp_to_A[problem->pol->n_upp] = j;
            problem->pol->A_to_Aupp[j]                = problem->pol->n_upp++;
        } else {
            problem->pol->A_to_Aupp[j] = -1;
        }
    }

    // Check if there are no active constraints
    if (problem->pol->n_low + problem->pol->n_upp == 0) {
        // Form empty Ared
        problem->pol->Ared = csc_spalloc(0, problem->data->n, 0, 1, 0);
        if (!(problem->pol->Ared)) return -1;
        int_vec_set_scalar(problem->pol->Ared->p, 0, problem->data->n + 1);
        return 0; // mred = 0
    }

    // Count number of elements in Ared
    for (j = 0; j < problem->data->A->p[problem->data->A->n]; j++) {
        if ((problem->pol->A_to_Alow[problem->data->A->i[j]] != -1) ||
                (problem->pol->A_to_Aupp[problem->data->A->i[j]] != -1)) Ared_nnz++;
    }

    // Form Ared
    // Ared = vstack[Alow, Aupp]
    problem->pol->Ared = csc_spalloc(problem->pol->n_low + problem->pol->n_upp,
                                     problem->data->n, Ared_nnz, 1, 0);
    if (!(problem->pol->Ared)) return -1;
    Ared_nnz = 0; // counter

    for (j = 0; j < problem->data->n; j++) { // Cycle over columns of A
        problem->pol->Ared->p[j] = Ared_nnz;

        for (ptr = problem->data->A->p[j]; ptr < problem->data->A->p[j + 1]; ptr++) {
            // Cycle over elements in j-th column
            if (problem->pol->A_to_Alow[problem->data->A->i[ptr]] != -1) {
                // Lower-active rows of A
                problem->pol->Ared->i[Ared_nnz] =
                    problem->pol->A_to_Alow[problem->data->A->i[ptr]];
                problem->pol->Ared->x[Ared_nnz++] = problem->data->A->x[ptr];
            } else if (problem->pol->A_to_Aupp[problem->data->A->i[ptr]] != -1) {
                // Upper-active rows of A
                problem->pol->Ared->i[Ared_nnz] = problem->pol->A_to_Aupp[problem->data->A->i[ptr]] \
                                                  + problem->pol->n_low;
                problem->pol->Ared->x[Ared_nnz++] = problem->data->A->x[ptr];
            }
        }
    }

    // Update the last element in Ared->p
    problem->pol->Ared->p[problem->data->n] = Ared_nnz;

    // Return number of rows in Ared
    return problem->pol->n_low + problem->pol->n_upp;
}

/**
 * Form reduced right-hand side rhs_red = vstack[-q, l_low, u_upp]
 * @param  problem problemspace
 * @param  rhs  right-hand-side
 * @return      reduced rhs
 */
static void form_rhs_red(qpInstance *problem, float *rhs) {
    int j;

    // Form the rhs of the reduced KKT linear system
    for (j = 0; j < problem->data->n; j++) { // -q
        rhs[j] = -problem->data->q[j];
    }

    for (j = 0; j < problem->pol->n_low; j++) { // l_low
        rhs[problem->data->n + j] = problem->data->l[problem->pol->Alow_to_A[j]];
    }

    for (j = 0; j < problem->pol->n_upp; j++) { // u_upp
        rhs[problem->data->n + problem->pol->n_low + j] =
            problem->data->u[problem->pol->Aupp_to_A[j]];
    }
}

static int iterative_refinement(qpInstance *problem,
                                LinSysSolver  *p,
                                float       *z,
                                float       *b) {
    int i, j, n;
    float *rhs;

    if (problem->settings->polish_refine_iter > 0) {

        // Assign dimension n
        n = problem->data->n + problem->pol->Ared->m;

        // Allocate rhs vector
        rhs = (float *)c_malloc(sizeof(float) * n);

        if (!rhs) {
            //   return osqp_error(OSQP_MEM_ALLOC_ERROR);
            return -5;
        } else {
            for (i = 0; i < problem->settings->polish_refine_iter; i++) {
                // Form the RHS for the iterative refinement:  b - K*z
                prea_vec_copy(b, rhs, n);

                // Upper Part: R^{n}
                // -= Px (upper triang)
                mat_vec(problem->data->P, z, rhs, -1);

                // -= Px (lower triang)
                mat_tpose_vec(problem->data->P, z, rhs, -1, 1);

                // -= Ared'*y_red
                mat_tpose_vec(problem->pol->Ared, z + problem->data->n, rhs, -1, 0);

                // Lower Part: R^{m}
                mat_vec(problem->pol->Ared, z, rhs + problem->data->n, -1);

                // Solve linear system. Store solution in rhs
                p->solve(p, rhs);

                // Update solution
                for (j = 0; j < n; j++) {
                    z[j] += rhs[j];
                }
            }
        }
        if (rhs) c_free(rhs);
    }
    return 0;
}
/**
 * Compute dual variable y from yred
 * @param problem problemspace
 * @param yred Dual variables associated to active constraints
 */
static void get_ypol_from_yred(qpInstance *problem, float *yred) {
    int j;

    // If there are no active constraints
    if (problem->pol->n_low + problem->pol->n_upp == 0) {
        vec_set_scalar(problem->pol->y, 0., problem->data->m);
        return;
    }

    // NB: yred = vstack[ylow, yupp]
    for (j = 0; j < problem->data->m; j++) {
        if (problem->pol->A_to_Alow[j] != -1) {
            // lower-active
            problem->pol->y[j] = yred[problem->pol->A_to_Alow[j]];
        } else if (problem->pol->A_to_Aupp[j] != -1) {
            // upper-active
            problem->pol->y[j] = yred[problem->pol->A_to_Aupp[j] + problem->pol->n_low];
        } else {
            // inactive
            problem->pol->y[j] = 0.0;
        }
    }
}

int polish(qpInstance *problem) {
    int mred, polish_successful, exitflag;
    float *rhs_red;
    LinSysSolver *plsh;
    float *pol_sol; // Polished solution

    qp_tic(problem->timer); // Start timer

    // Form Ared by assuming the active constraints and store in problem->pol->Ared
    mred = form_Ared(problem);
    if (mred < 0) { // problem->pol->red = OSQP_NULL
        // Polishing failed
        problem->info->status_polish = -1;

        return -1;
    }

    // Form and factorize reduced KKT -> FIX THIS
    exitflag = init_linsys_solver(&plsh, problem->data->P, problem->pol->Ared,
                                  problem->settings->delta, 0,
                                  problem->settings->linsys_solver, 1);

    if (exitflag) {
        // Polishing failed
        problem->info->status_polish = -1;

        // Memory clean-up
        if (problem->pol->Ared) csc_spfree(problem->pol->Ared);

        return 1;
    }

    // Form reduced right-hand side rhs_red
    rhs_red = c_malloc(sizeof(float) * (problem->data->n + mred));
    if (!rhs_red) {
        // Polishing failed
        problem->info->status_polish = -1;

        // Memory clean-up
        csc_spfree(problem->pol->Ared);

        return -1;
    }
    form_rhs_red(problem, rhs_red);

    pol_sol = vec_copy(rhs_red, problem->data->n + mred);
    if (!pol_sol) {
        // Polishing failed
        problem->info->status_polish = -1;

        // Memory clean-up
        csc_spfree(problem->pol->Ared);
        c_free(rhs_red);

        return -1;
    }

    // Solve the reduced KKT system
    plsh->solve(plsh, pol_sol);

    // Perform iterative refinement to compensate for the regularization error
    exitflag = iterative_refinement(problem, plsh, pol_sol, rhs_red);

    if (exitflag) {
        // Polishing failed
        problem->info->status_polish = -1;

        // Memory clean-up
        csc_spfree(problem->pol->Ared);
        c_free(rhs_red);
        c_free(pol_sol);

        return -1;
    }

    // Store the polished solution (x,z,y)
    prea_vec_copy(pol_sol, problem->pol->x, problem->data->n);   // pol->x
    mat_vec(problem->data->A, problem->pol->x, problem->pol->z, 0); // pol->z
    get_ypol_from_yred(problem, pol_sol + problem->data->n);     // pol->y

    // Ensure (z,y) satisfies normal cone constraint
    project_normalcone(problem, problem->pol->z, problem->pol->y);

    // Compute primal and dual residuals at the polished solution
    update_info(problem, 0, 1, 1);

    // Check if polish was successful
    polish_successful = (problem->pol->pri_res < problem->info->pri_res &&
                         problem->pol->dua_res < problem->info->dua_res) || // Residuals
                        // are
                        // reduced
                        (problem->pol->pri_res < problem->info->pri_res &&
                         problem->info->dua_res < 1e-10) ||              // Dual
                        // residual
                        // already
                        // tiny
                        (problem->pol->dua_res < problem->info->dua_res &&
                         problem->info->pri_res < 1e-10);                // Primal
    // residual
    // already
    // tiny

    if (polish_successful) {
        // Update solver information
        problem->info->obj_val       = problem->pol->obj_val;
        problem->info->pri_res       = problem->pol->pri_res;
        problem->info->dua_res       = problem->pol->dua_res;
        problem->info->status_polish = 1;

        // Update (x, z, y) in ADMM iterations
        // NB: z needed for warm starting
        prea_vec_copy(problem->pol->x, problem->x, problem->data->n);
        prea_vec_copy(problem->pol->z, problem->z, problem->data->m);
        prea_vec_copy(problem->pol->y, problem->y, problem->data->m);


    } else { // Polishing failed
        problem->info->status_polish = -1;

        // TODO: Try to find a better solution on the line connecting ADMM
        //       and polished solution
    }

    // Memory clean-up
    plsh->free(plsh);

    // Checks that they are not NULL are already performed earlier
    csc_spfree(problem->pol->Ared);
    c_free(rhs_red);
    c_free(pol_sol);

    return 0;
}















/*
***********************************************************
* Auxiliary functions needed to compute ADMM iterations * *
***********************************************************/
float compute_rho_estimate(qpInstance *problem)
{
    int   n, m;                       // Dimensions
    float pri_res, dua_res;           // Primal and dual residuals
    float pri_res_norm, dua_res_norm; // Normalization for the residuals
    float temp_res_norm;              // Temporary residual norm
    float rho_estimate;               // Rho estimate value

    // Get problem dimensions
    n = problem->data->n;
    m = problem->data->m;

    // Get primal and dual residuals
    pri_res = vec_norm_inf(problem->z_prev, m);
    dua_res = vec_norm_inf(problem->x_prev, n);

    // Normalize primal residual
    pri_res_norm  = vec_norm_inf(problem->z, m);           // ||z||
    temp_res_norm = vec_norm_inf(problem->Ax, m);          // ||Ax||
    pri_res_norm  = c_max(pri_res_norm, temp_res_norm); // max (||z||,||Ax||)
    pri_res      /= (pri_res_norm + 1e-10);             // Normalize primal
    // residual (prevent 0
    // division)

    // Normalize dual residual
    dua_res_norm  = vec_norm_inf(problem->data->q, n);     // ||q||
    temp_res_norm = vec_norm_inf(problem->Aty, n);         // ||A' y||
    dua_res_norm  = c_max(dua_res_norm, temp_res_norm);
    temp_res_norm = vec_norm_inf(problem->Px, n);          //  ||P x||
    dua_res_norm  = c_max(dua_res_norm, temp_res_norm); // max(||q||,||A' y||,||P
    // x||)
    dua_res      /= (dua_res_norm + 1e-10);             // Normalize dual residual
    // (prevent 0 division)


    // Return rho estimate
    rho_estimate = problem->settings->rho * c_sqrt(pri_res / (dua_res + 1e-10)); // (prevent
    // 0
    // division)
    rho_estimate = c_min(c_max(rho_estimate, RHO_MIN), RHO_MAX);              // Constrain
    // rho
    // values
    return rho_estimate;
}

int adapt_rho(qpInstance *problem) {
    int   exitflag; // Exitflag
    float rho_new;  // New rho value

    exitflag = 0;     // Initialize exitflag to 0

    // Compute new rho
    rho_new = compute_rho_estimate(problem);

    // Set rho estimate in info
    problem->info->rho_estimate = rho_new;

    // Check if the new rho is large or small enough and update it in case
    if ((rho_new > problem->settings->rho * problem->settings->adaptive_rho_tolerance) ||
            (rho_new < problem->settings->rho /  problem->settings->adaptive_rho_tolerance)) {
        exitflag                 = osqp_update_rho(problem, rho_new);
        problem->info->rho_updates += 1;
    }

    return exitflag;
}

void set_rho_vec(qpInstance *problem) {
    int i;

    problem->settings->rho = c_min(c_max(problem->settings->rho, RHO_MIN), RHO_MAX);

    for (i = 0; i < problem->data->m; i++) {
        if ((problem->data->l[i] < -INFTY * MIN_SCALING) &&
                (problem->data->u[i] > INFTY * MIN_SCALING)) {
            // Loose bounds
            problem->constr_type[i] = -1;
            problem->rho_vec[i]     = RHO_MIN;
        } else if (problem->data->u[i] - problem->data->l[i] < RHO_TOL) {
            // Equality constraints
            problem->constr_type[i] = 1;
            problem->rho_vec[i]     = RHO_EQ_OVER_RHO_INEQ * problem->settings->rho;
        } else {
            // Inequality constraints
            problem->constr_type[i] = 0;
            problem->rho_vec[i]     = problem->settings->rho;
        }
        problem->rho_inv_vec[i] = 1. / problem->rho_vec[i];
    }
}

int update_rho_vec(qpInstance *problem) {
    int i, exitflag, constr_type_changed;

    exitflag            = 0;
    constr_type_changed = 0;

    for (i = 0; i < problem->data->m; i++) {
        if ((problem->data->l[i] < -INFTY * MIN_SCALING) &&
                (problem->data->u[i] > INFTY * MIN_SCALING)) {
            // Loose bounds
            if (problem->constr_type[i] != -1) {
                problem->constr_type[i] = -1;
                problem->rho_vec[i]     = RHO_MIN;
                problem->rho_inv_vec[i] = 1. / RHO_MIN;
                constr_type_changed  = 1;
            }
        } else if (problem->data->u[i] - problem->data->l[i] < RHO_TOL) {
            // Equality constraints
            if (problem->constr_type[i] != 1) {
                problem->constr_type[i] = 1;
                problem->rho_vec[i]     = RHO_EQ_OVER_RHO_INEQ * problem->settings->rho;
                problem->rho_inv_vec[i] = 1. / problem->rho_vec[i];
                constr_type_changed  = 1;
            }
        } else {
            // Inequality constraints
            if (problem->constr_type[i] != 0) {
                problem->constr_type[i] = 0;
                problem->rho_vec[i]     = problem->settings->rho;
                problem->rho_inv_vec[i] = 1. / problem->settings->rho;
                constr_type_changed  = 1;
            }
        }
    }

    // Update rho_vec in KKT matrix if constraints type has changed
    if (constr_type_changed == 1) {
        exitflag = problem->linsys_solver->update_rho_vec(problem->linsys_solver,
                   problem->rho_vec);
    }

    return exitflag;
}

void swap_vectors(float **a, float **b) {
    float *temp;

    temp = *b;
    *b   = *a;
    *a   = temp;
}

void cold_start(qpInstance *problem) {
    vec_set_scalar(problem->x, 0., problem->data->n);
    vec_set_scalar(problem->z, 0., problem->data->m);
    vec_set_scalar(problem->y, 0., problem->data->m);
}

static void compute_rhs(qpInstance *problem) {
    int i; // Index

    for (i = 0; i < problem->data->n; i++) {
        // Cycle over part related to x variables
        problem->xz_tilde[i] = problem->settings->sigma * problem->x_prev[i] -
                               problem->data->q[i];
    }

    for (i = 0; i < problem->data->m; i++) {
        // Cycle over dual variable in the first step (nu)
        problem->xz_tilde[i + problem->data->n] = problem->z_prev[i] - problem->rho_inv_vec[i] *
                problem->y[i];
    }
}

void update_xz_tilde(qpInstance *problem) {
    // Compute right-hand side
    compute_rhs(problem);

    // Solve linear system
    problem->linsys_solver->solve(problem->linsys_solver, problem->xz_tilde);
}

void update_x(qpInstance *problem) {
    int i;

    // update x
    for (i = 0; i < problem->data->n; i++) {
        problem->x[i] = problem->settings->alpha * problem->xz_tilde[i] +
                        ((float)1.0 - problem->settings->alpha) * problem->x_prev[i];
    }

    // update delta_x
    for (i = 0; i < problem->data->n; i++) {
        problem->delta_x[i] = problem->x[i] - problem->x_prev[i];
    }
}

void update_z(qpInstance *problem) {
    int i;

    // update z
    for (i = 0; i < problem->data->m; i++) {
        problem->z[i] = problem->settings->alpha * problem->xz_tilde[i + problem->data->n] +
                        ((float)1.0 - problem->settings->alpha) * problem->z_prev[i] +
                        problem->rho_inv_vec[i] * problem->y[i];
    }

    // project z
    project(problem, problem->z);
}

void update_y(qpInstance *problem) {
    int i; // Index

    for (i = 0; i < problem->data->m; i++) {
        problem->delta_y[i] = problem->rho_vec[i] *
                              (problem->settings->alpha *
                               problem->xz_tilde[i + problem->data->n] +
                               ((float)1.0 - problem->settings->alpha) * problem->z_prev[i] -
                               problem->z[i]);
        problem->y[i] += problem->delta_y[i];
    }
}

float compute_obj_val(qpInstance *problem, float *x) {
    float obj_val;

    obj_val = quad_form(problem->data->P, x) +
              vec_prod(problem->data->q, x, problem->data->n);

    if (problem->settings->scaling) {
        obj_val *= problem->scaling->cinv;
    }

    return obj_val;
}

float compute_pri_res(qpInstance *problem, float *x, float *z) {
    // NB: Use z_prev as probleming vector
    // pr = Ax - z

    mat_vec(problem->data->A, x, problem->Ax, 0); // Ax
    vec_add_scaled(problem->z_prev, problem->Ax, z, problem->data->m, -1);

    // If scaling active -> rescale residual
    if (problem->settings->scaling && !problem->settings->scaled_termination) {
        return vec_scaled_norm_inf(problem->scaling->Einv, problem->z_prev, problem->data->m);
    }

    // Return norm of the residual
    return vec_norm_inf(problem->z_prev, problem->data->m);
}

float compute_pri_tol(qpInstance *problem, float eps_abs, float eps_rel) {
    float max_rel_eps, temp_rel_eps;

    // max_rel_eps = max(||z||, ||A x||)
    if (problem->settings->scaling && !problem->settings->scaled_termination) {
        // ||Einv * z||
        max_rel_eps =
            vec_scaled_norm_inf(problem->scaling->Einv, problem->z, problem->data->m);

        // ||Einv * A * x||
        temp_rel_eps = vec_scaled_norm_inf(problem->scaling->Einv,
                                           problem->Ax,
                                           problem->data->m);

        // Choose maximum
        max_rel_eps = c_max(max_rel_eps, temp_rel_eps);
    } else { // No unscaling required
        // ||z||
        max_rel_eps = vec_norm_inf(problem->z, problem->data->m);

        // ||A * x||
        temp_rel_eps = vec_norm_inf(problem->Ax, problem->data->m);

        // Choose maximum
        max_rel_eps = c_max(max_rel_eps, temp_rel_eps);
    }

    // eps_prim
    return eps_abs + eps_rel * max_rel_eps;
}

float compute_dua_res(qpInstance *problem, float *x, float *y) {
    // NB: Use x_prev as temporary vector
    // NB: Only upper triangular part of P is stored.
    // dr = q + A'*y + P*x

    // dr = q
    prea_vec_copy(problem->data->q, problem->x_prev, problem->data->n);

    // P * x (upper triangular part)
    mat_vec(problem->data->P, x, problem->Px, 0);

    // P' * x (lower triangular part with no diagonal)
    mat_tpose_vec(problem->data->P, x, problem->Px, 1, 1);

    // dr += P * x (full P matrix)
    vec_add_scaled(problem->x_prev, problem->x_prev, problem->Px, problem->data->n, 1);

    // dr += A' * y
    if (problem->data->m > 0) {
        mat_tpose_vec(problem->data->A, y, problem->Aty, 0, 0);
        vec_add_scaled(problem->x_prev, problem->x_prev, problem->Aty, problem->data->n, 1);
    }

    // If scaling active -> rescale residual
    if (problem->settings->scaling && !problem->settings->scaled_termination) {
        return problem->scaling->cinv * vec_scaled_norm_inf(problem->scaling->Dinv,
                problem->x_prev,
                problem->data->n);
    }

    return vec_norm_inf(problem->x_prev, problem->data->n);
}

float compute_dua_tol(qpInstance *problem, float eps_abs, float eps_rel) {
    float max_rel_eps, temp_rel_eps;

    // max_rel_eps = max(||q||, ||A' y|, ||P x||)
    if (problem->settings->scaling && !problem->settings->scaled_termination) {
        // || Dinv q||
        max_rel_eps = vec_scaled_norm_inf(problem->scaling->Dinv,
                                          problem->data->q,
                                          problem->data->n);

        // || Dinv A' y ||
        temp_rel_eps = vec_scaled_norm_inf(problem->scaling->Dinv,
                                           problem->Aty,
                                           problem->data->n);
        max_rel_eps = c_max(max_rel_eps, temp_rel_eps);

        // || Dinv P x||
        temp_rel_eps = vec_scaled_norm_inf(problem->scaling->Dinv,
                                           problem->Px,
                                           problem->data->n);
        max_rel_eps = c_max(max_rel_eps, temp_rel_eps);

        // Multiply by cinv
        max_rel_eps *= problem->scaling->cinv;
    } else { // No scaling required
        // ||q||
        max_rel_eps = vec_norm_inf(problem->data->q, problem->data->n);

        // ||A'*y||
        temp_rel_eps = vec_norm_inf(problem->Aty, problem->data->n);
        max_rel_eps  = c_max(max_rel_eps, temp_rel_eps);

        // ||P*x||
        temp_rel_eps = vec_norm_inf(problem->Px, problem->data->n);
        max_rel_eps  = c_max(max_rel_eps, temp_rel_eps);
    }

    // eps_dual
    return eps_abs + eps_rel * max_rel_eps;
}

int is_primal_infeasible(qpInstance *problem, float eps_prim_inf) {
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
    for (i = 0; i < problem->data->m; i++) {
        if (problem->data->u[i] > INFTY * MIN_SCALING) {          // Infinite upper bound
            if (problem->data->l[i] < -INFTY * MIN_SCALING) {       // Infinite lower bound
                // Both bounds infinite
                problem->delta_y[i] = 0.0;
            } else {
                // Only upper bound infinite
                problem->delta_y[i] = c_min(problem->delta_y[i], 0.0);
            }
        } else if (problem->data->l[i] < -INFTY * MIN_SCALING) {  // Infinite lower bound
            // Only lower bound infinite
            problem->delta_y[i] = c_max(problem->delta_y[i], 0.0);
        }
    }

    // Compute infinity norm of delta_y (unscale if necessary)
    if (problem->settings->scaling && !problem->settings->scaled_termination) {
        // Use problem->Adelta_x as temporary vector
        vec_ew_prod(problem->scaling->E, problem->delta_y, problem->Adelta_x, problem->data->m);
        norm_delta_y = vec_norm_inf(problem->Adelta_x, problem->data->m);
    } else {
        norm_delta_y = vec_norm_inf(problem->delta_y, problem->data->m);
    }

    if (norm_delta_y > eps_prim_inf) { // ||delta_y|| > 0

        for (i = 0; i < problem->data->m; i++) {
            ineq_lhs += problem->data->u[i] * c_max(problem->delta_y[i], 0) + \
                        problem->data->l[i] * c_min(problem->delta_y[i], 0);
        }

        // Check if the condition is satisfied: ineq_lhs < -eps
        if (ineq_lhs < -eps_prim_inf * norm_delta_y) {
            // Compute and return ||A'delta_y|| < eps_prim_inf
            mat_tpose_vec(problem->data->A, problem->delta_y, problem->Atdelta_y, 0, 0);

            // Unscale if necessary
            if (problem->settings->scaling && !problem->settings->scaled_termination) {
                vec_ew_prod(problem->scaling->Dinv,
                            problem->Atdelta_y,
                            problem->Atdelta_y,
                            problem->data->n);
            }

            return vec_norm_inf(problem->Atdelta_y, problem->data->n) < eps_prim_inf * norm_delta_y;
        }
    }

    // Conditions not satisfied -> not primal infeasible
    return 0;
}

int is_dual_infeasible(qpInstance *problem, float eps_dual_inf) {
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
    if (problem->settings->scaling && !problem->settings->scaled_termination) { // Unscale
        // if
        // necessary
        norm_delta_x = vec_scaled_norm_inf(problem->scaling->D,
                                           problem->delta_x,
                                           problem->data->n);
        cost_scaling = problem->scaling->c;
    } else {
        norm_delta_x = vec_norm_inf(problem->delta_x, problem->data->n);
        cost_scaling = 1.0;
    }

    // Prevent 0 division || delta_x || > 0
    if (norm_delta_x > eps_dual_inf) {
        // Normalize delta_x by its norm

        /* vec_mult_scalar(problem->delta_x, 1./norm_delta_x, problem->data->n); */

        // Check first if q'*delta_x < 0
        if (vec_prod(problem->data->q, problem->delta_x, problem->data->n) <
                -cost_scaling * eps_dual_inf * norm_delta_x) {
            // Compute product P * delta_x (NB: P is store in upper triangular form)
            mat_vec(problem->data->P, problem->delta_x, problem->Pdelta_x, 0);
            mat_tpose_vec(problem->data->P, problem->delta_x, problem->Pdelta_x, 1, 1);

            // Scale if necessary
            if (problem->settings->scaling && !problem->settings->scaled_termination) {
                vec_ew_prod(problem->scaling->Dinv,
                            problem->Pdelta_x,
                            problem->Pdelta_x,
                            problem->data->n);
            }

            // Check if || P * delta_x || = 0
            if (vec_norm_inf(problem->Pdelta_x, problem->data->n) <
                    cost_scaling * eps_dual_inf * norm_delta_x) {
                // Compute A * delta_x
                mat_vec(problem->data->A, problem->delta_x, problem->Adelta_x, 0);

                // Scale if necessary
                if (problem->settings->scaling && !problem->settings->scaled_termination) {
                    vec_ew_prod(problem->scaling->Einv,
                                problem->Adelta_x,
                                problem->Adelta_x,
                                problem->data->m);
                }

                // De Morgan Law Applied to dual infeasibility conditions for A * x
                // NB: Note that MIN_SCALING is used to adjust the infinity value
                //     in case the problem is scaled.
                for (i = 0; i < problem->data->m; i++) {
                    if (((problem->data->u[i] < INFTY * MIN_SCALING) &&
                            (problem->Adelta_x[i] >  eps_dual_inf * norm_delta_x)) ||
                            ((problem->data->l[i] > -INFTY * MIN_SCALING) &&
                             (problem->Adelta_x[i] < -eps_dual_inf * norm_delta_x))) {
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

    return ((info->status_val != OSQP_PRIMAL_INFEASIBLE) &&
            (info->status_val != OSQP_PRIMAL_INFEASIBLE_INACCURATE) &&
            (info->status_val != OSQP_DUAL_INFEASIBLE) &&
            (info->status_val != OSQP_DUAL_INFEASIBLE_INACCURATE) &&
            (info->status_val != OSQP_NON_CVX));

}

void store_solution(qpInstance *problem) {
    float norm_vec;

    if (has_solution(problem->info)) {
        prea_vec_copy(problem->x, problem->solution->x, problem->data->n); // primal
        prea_vec_copy(problem->y, problem->solution->y, problem->data->m); // dual

        // Unscale solution if scaling has been performed
        if (problem->settings->scaling)
            unscale_solution(problem);
    } else {
        // No solution present. Solution is NaN
        vec_set_scalar(problem->solution->x, qpNAN, problem->data->n);
        vec_set_scalar(problem->solution->y, qpNAN, problem->data->m);


        // Normalize infeasibility certificates if embedded is off
        // NB: It requires a division
        if ((problem->info->status_val == OSQP_PRIMAL_INFEASIBLE) ||
                ((problem->info->status_val == OSQP_PRIMAL_INFEASIBLE_INACCURATE))) {
            norm_vec = vec_norm_inf(problem->delta_y, problem->data->m);
            vec_mult_scalar(problem->delta_y, 1. / norm_vec, problem->data->m);
        }

        if ((problem->info->status_val == OSQP_DUAL_INFEASIBLE) ||
                ((problem->info->status_val == OSQP_DUAL_INFEASIBLE_INACCURATE))) {
            norm_vec = vec_norm_inf(problem->delta_x, problem->data->n);
            vec_mult_scalar(problem->delta_x, 1. / norm_vec, problem->data->n);
        }


        // Cold start iterates to 0 for next runs (they cannot start from NaN)
        cold_start(problem);
    }
}

void update_info(qpInstance *problem,
                 int          iter,
                 int          compute_objective,
                 int          polish) {
    float *x, *z, *y;                   // Allocate pointers to variables
    float *obj_val, *pri_res, *dua_res; // objective value, residuals

    float *run_time;                    // Execution time


    if (polish) {
        x       = problem->pol->x;
        y       = problem->pol->y;
        z       = problem->pol->z;
        obj_val = &problem->pol->obj_val;
        pri_res = &problem->pol->pri_res;
        dua_res = &problem->pol->dua_res;
        run_time = &problem->info->polish_time;
    } else {
        x                = problem->x;
        y                = problem->y;
        z                = problem->z;
        obj_val          = &problem->info->obj_val;
        pri_res          = &problem->info->pri_res;
        dua_res          = &problem->info->dua_res;
        problem->info->iter = iter; // Update iteration number
        run_time = &problem->info->solve_time;
    }



    // Compute the objective if needed
    if (compute_objective) {
        *obj_val = compute_obj_val(problem, x);
    }

    // Compute primal residual
    if (problem->data->m == 0) {
        // No constraints -> Always primal feasible
        *pri_res = 0.;
    } else {
        *pri_res = compute_pri_res(problem, x, z);
    }

    // Compute dual residual
    *dua_res = compute_dua_res(problem, x, y);

    // Update timing
    *run_time = qp_toc(problem->timer);

}


void reset_info(qpInfo *info) {

    // Initialize info values.
    info->solve_time = 0.0;  // Solve time to zero
    info->polish_time = 0.0; // Polish time to zero



    update_status(info, OSQP_UNSOLVED); // Problem is unsolved

    info->rho_updates = 0;              // Rho updates are now 0
}

void update_status(qpInfo *info, int status_val) {
    // Update status value
    info->status_val = status_val;

    // Update status string depending on status val
    if (status_val == OSQP_SOLVED) c_strcpy(info->status, "solved");

    if (status_val == OSQP_SOLVED_INACCURATE) c_strcpy(info->status,
                "solved inaccurate");
    else if (status_val == OSQP_PRIMAL_INFEASIBLE) c_strcpy(info->status,
                "primal infeasible");
    else if (status_val == OSQP_PRIMAL_INFEASIBLE_INACCURATE) c_strcpy(info->status,
                "primal infeasible inaccurate");
    else if (status_val == OSQP_UNSOLVED) c_strcpy(info->status, "unsolved");
    else if (status_val == OSQP_DUAL_INFEASIBLE) c_strcpy(info->status,
                "dual infeasible");
    else if (status_val == OSQP_DUAL_INFEASIBLE_INACCURATE) c_strcpy(info->status,
                "dual infeasible inaccurate");
    else if (status_val == OSQP_MAX_ITER_REACHED) c_strcpy(info->status,
                "maximum iterations reached");
    else if (status_val == OSQP_TIME_LIMIT_REACHED) c_strcpy(info->status,
                "run time limit reached");
    else if (status_val == OSQP_SIGINT) c_strcpy(info->status, "interrupted");

    else if (status_val == OSQP_NON_CVX) c_strcpy(info->status, "problem non convex");

}

int check_termination(qpInstance *problem, int approximate) {
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
    eps_abs      = problem->settings->eps_abs;
    eps_rel      = problem->settings->eps_rel;
    eps_prim_inf = problem->settings->eps_prim_inf;
    eps_dual_inf = problem->settings->eps_dual_inf;

    // If residuals are too large, the problem is probably non convex
    if ((problem->info->pri_res > INFTY) ||
            (problem->info->dua_res > INFTY)) {
        // Looks like residuals are diverging. Probably the problem is non convex!
        // Terminate and report it
        update_status(problem->info, OSQP_NON_CVX);
        problem->info->obj_val = qpNAN;
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
    if (problem->data->m == 0) {
        prim_res_check = 1; // No constraints -> Primal feasibility always satisfied
    }
    else {
        // Compute primal tolerance
        eps_prim = compute_pri_tol(problem, eps_abs, eps_rel);

        // Primal feasibility check
        if (problem->info->pri_res < eps_prim) {
            prim_res_check = 1;
        } else {
            // Primal infeasibility check
            prim_inf_check = is_primal_infeasible(problem, eps_prim_inf);
        }
    } // End check if m == 0

    // Compute dual tolerance
    eps_dual = compute_dua_tol(problem, eps_abs, eps_rel);

    // Dual feasibility check
    if (problem->info->dua_res < eps_dual) {
        dual_res_check = 1;
    } else {
        // Check dual infeasibility
        dual_inf_check = is_dual_infeasible(problem, eps_dual_inf);
    }

    // Compare checks to determine solver status
    if (prim_res_check && dual_res_check) {
        // Update final information
        if (approximate) {
            update_status(problem->info, OSQP_SOLVED_INACCURATE);
        } else {
            update_status(problem->info, OSQP_SOLVED);
        }
        exitflag = 1;
    }
    else if (prim_inf_check) {
        // Update final information
        if (approximate) {
            update_status(problem->info, OSQP_PRIMAL_INFEASIBLE_INACCURATE);
        } else {
            update_status(problem->info, OSQP_PRIMAL_INFEASIBLE);
        }

        if (problem->settings->scaling && !problem->settings->scaled_termination) {
            // Update infeasibility certificate
            vec_ew_prod(problem->scaling->E, problem->delta_y, problem->delta_y, problem->data->m);
        }
        problem->info->obj_val = INFTY;
        exitflag            = 1;
    }
    else if (dual_inf_check) {
        // Update final information
        if (approximate) {
            update_status(problem->info, OSQP_DUAL_INFEASIBLE_INACCURATE);
        } else {
            update_status(problem->info, OSQP_DUAL_INFEASIBLE);
        }

        if (problem->settings->scaling && !problem->settings->scaled_termination) {
            // Update infeasibility certificate
            vec_ew_prod(problem->scaling->D, problem->delta_x, problem->delta_x, problem->data->n);
        }
        problem->info->obj_val = -INFTY;
        exitflag            = 1;
    }

    return exitflag;
}


