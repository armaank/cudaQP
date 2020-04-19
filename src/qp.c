/* qp.c */
#include "./include/qp.h"
#include "./include/ops.h"
#include "./include/linalg.h"
#include "./include/kkt.h"
#include "./include/csc.h"
#include "./include/timer.h"
#include "./include/params.h"
#include "./include/error.h"



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


    for (i = 0; i < problem->hparams->scaling; i++) {
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

    if (problem->hparams->polish_refine_iter > 0) {

        // Assign dimension n
        n = problem->data->n + problem->pol->Ared->m;

        // Allocate rhs vector
        rhs = (float *)c_malloc(sizeof(float) * n);

        if (!rhs) {
            //   return qp_error(QP_MEM_ALLOC_ERROR);
            return -5;
        } else {
            for (i = 0; i < problem->hparams->polish_refine_iter; i++) {
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
    if (mred < 0) { // problem->pol->red = QP_NULL
        // Polishing failed
        problem->info->status_polish = -1;

        return -1;
    }

    // Form and factorize reduced KKT -> FIX THIS
    exitflag = init_linsys_solver(&plsh, problem->data->P, problem->pol->Ared,
                                  problem->hparams->delta, 0,
                                  problem->hparams->linsys_solver, 1);

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
    rho_estimate = problem->hparams->rho * c_sqrt(pri_res / (dua_res + 1e-10)); // (prevent
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
    if ((rho_new > problem->hparams->rho * problem->hparams->adaptive_rho_tolerance) ||
            (rho_new < problem->hparams->rho /  problem->hparams->adaptive_rho_tolerance)) {
        exitflag                 = qp_update_rho(problem, rho_new);
        problem->info->rho_updates += 1;
    }

    return exitflag;
}

void set_rho_vec(qpInstance *problem) {
    int i;

    problem->hparams->rho = c_min(c_max(problem->hparams->rho, RHO_MIN), RHO_MAX);

    for (i = 0; i < problem->data->m; i++) {
        if ((problem->data->l[i] < -INFTY * MIN_SCALING) &&
                (problem->data->u[i] > INFTY * MIN_SCALING)) {
            // Loose bounds
            problem->constr_type[i] = -1;
            problem->rho_vec[i]     = RHO_MIN;
        } else if (problem->data->u[i] - problem->data->l[i] < RHO_TOL) {
            // Equality constraints
            problem->constr_type[i] = 1;
            problem->rho_vec[i]     = RHO_EQ_OVER_RHO_INEQ * problem->hparams->rho;
        } else {
            // Inequality constraints
            problem->constr_type[i] = 0;
            problem->rho_vec[i]     = problem->hparams->rho;
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
                problem->rho_vec[i]     = RHO_EQ_OVER_RHO_INEQ * problem->hparams->rho;
                problem->rho_inv_vec[i] = 1. / problem->rho_vec[i];
                constr_type_changed  = 1;
            }
        } else {
            // Inequality constraints
            if (problem->constr_type[i] != 0) {
                problem->constr_type[i] = 0;
                problem->rho_vec[i]     = problem->hparams->rho;
                problem->rho_inv_vec[i] = 1. / problem->hparams->rho;
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
        problem->xz_tilde[i] = problem->hparams->sigma * problem->x_prev[i] -
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
        problem->x[i] = problem->hparams->alpha * problem->xz_tilde[i] +
                        ((float)1.0 - problem->hparams->alpha) * problem->x_prev[i];
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
        problem->z[i] = problem->hparams->alpha * problem->xz_tilde[i + problem->data->n] +
                        ((float)1.0 - problem->hparams->alpha) * problem->z_prev[i] +
                        problem->rho_inv_vec[i] * problem->y[i];
    }

    // project z
    project(problem, problem->z);
}

void update_y(qpInstance *problem) {
    int i; // Index

    for (i = 0; i < problem->data->m; i++) {
        problem->delta_y[i] = problem->rho_vec[i] *
                              (problem->hparams->alpha *
                               problem->xz_tilde[i + problem->data->n] +
                               ((float)1.0 - problem->hparams->alpha) * problem->z_prev[i] -
                               problem->z[i]);
        problem->y[i] += problem->delta_y[i];
    }
}

float compute_obj_val(qpInstance *problem, float *x) {
    float obj_val;

    obj_val = quad_form(problem->data->P, x) +
              vec_prod(problem->data->q, x, problem->data->n);

    if (problem->hparams->scaling) {
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
    if (problem->hparams->scaling && !problem->hparams->scaled_termination) {
        return vec_scaled_norm_inf(problem->scaling->Einv, problem->z_prev, problem->data->m);
    }

    // Return norm of the residual
    return vec_norm_inf(problem->z_prev, problem->data->m);
}

float compute_pri_tol(qpInstance *problem, float eps_abs, float eps_rel) {
    float max_rel_eps, temp_rel_eps;

    // max_rel_eps = max(||z||, ||A x||)
    if (problem->hparams->scaling && !problem->hparams->scaled_termination) {
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
    if (problem->hparams->scaling && !problem->hparams->scaled_termination) {
        return problem->scaling->cinv * vec_scaled_norm_inf(problem->scaling->Dinv,
                problem->x_prev,
                problem->data->n);
    }

    return vec_norm_inf(problem->x_prev, problem->data->n);
}

float compute_dua_tol(qpInstance *problem, float eps_abs, float eps_rel) {
    float max_rel_eps, temp_rel_eps;

    // max_rel_eps = max(||q||, ||A' y|, ||P x||)
    if (problem->hparams->scaling && !problem->hparams->scaled_termination) {
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
    if (problem->hparams->scaling && !problem->hparams->scaled_termination) {
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
            if (problem->hparams->scaling && !problem->hparams->scaled_termination) {
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
    if (problem->hparams->scaling && !problem->hparams->scaled_termination) { // Unscale
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
            if (problem->hparams->scaling && !problem->hparams->scaled_termination) {
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
                if (problem->hparams->scaling && !problem->hparams->scaled_termination) {
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

    return ((info->status_val != QP_PRIMAL_INFEASIBLE) &&
            (info->status_val != QP_PRIMAL_INFEASIBLE_INACCURATE) &&
            (info->status_val != QP_DUAL_INFEASIBLE) &&
            (info->status_val != QP_DUAL_INFEASIBLE_INACCURATE) &&
            (info->status_val != QP_NON_CVX));

}

void store_solution(qpInstance *problem) {
    float norm_vec;

    if (has_solution(problem->info)) {
        prea_vec_copy(problem->x, problem->solution->x, problem->data->n); // primal
        prea_vec_copy(problem->y, problem->solution->y, problem->data->m); // dual

        // Unscale solution if scaling has been performed
        if (problem->hparams->scaling)
            unscale_solution(problem);
    } else {
        // No solution present. Solution is NaN
        vec_set_scalar(problem->solution->x, qpNAN, problem->data->n);
        vec_set_scalar(problem->solution->y, qpNAN, problem->data->m);


        // Normalize infeasibility certificates if embedded is off
        // NB: It requires a division
        if ((problem->info->status_val == QP_PRIMAL_INFEASIBLE) ||
                ((problem->info->status_val == QP_PRIMAL_INFEASIBLE_INACCURATE))) {
            norm_vec = vec_norm_inf(problem->delta_y, problem->data->m);
            vec_mult_scalar(problem->delta_y, 1. / norm_vec, problem->data->m);
        }

        if ((problem->info->status_val == QP_DUAL_INFEASIBLE) ||
                ((problem->info->status_val == QP_DUAL_INFEASIBLE_INACCURATE))) {
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



    update_status(info, QP_UNSOLVED); // Problem is unsolved

    info->rho_updates = 0;              // Rho updates are now 0
}

void update_status(qpInfo *info, int status_val) {
    // Update status value
    info->status_val = status_val;

    // Update status string depending on status val
    if (status_val == QP_SOLVED) c_strcpy(info->status, "solved");

    if (status_val == QP_SOLVED_INACCURATE) c_strcpy(info->status,
                "solved inaccurate");
    else if (status_val == QP_PRIMAL_INFEASIBLE) c_strcpy(info->status,
                "primal infeasible");
    else if (status_val == QP_PRIMAL_INFEASIBLE_INACCURATE) c_strcpy(info->status,
                "primal infeasible inaccurate");
    else if (status_val == QP_UNSOLVED) c_strcpy(info->status, "unsolved");
    else if (status_val == QP_DUAL_INFEASIBLE) c_strcpy(info->status,
                "dual infeasible");
    else if (status_val == QP_DUAL_INFEASIBLE_INACCURATE) c_strcpy(info->status,
                "dual infeasible inaccurate");
    else if (status_val == QP_MAX_ITER_REACHED) c_strcpy(info->status,
                "maximum iterations reached");
    else if (status_val == QP_TIME_LIMIT_REACHED) c_strcpy(info->status,
                "run time limit reached");
    else if (status_val == QP_SIGINT) c_strcpy(info->status, "interrupted");

    else if (status_val == QP_NON_CVX) c_strcpy(info->status, "problem non convex");

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
    eps_abs      = problem->hparams->eps_abs;
    eps_rel      = problem->hparams->eps_rel;
    eps_prim_inf = problem->hparams->eps_prim_inf;
    eps_dual_inf = problem->hparams->eps_dual_inf;

    // If residuals are too large, the problem is probably non convex
    if ((problem->info->pri_res > INFTY) ||
            (problem->info->dua_res > INFTY)) {
        // Looks like residuals are diverging. Probably the problem is non convex!
        // Terminate and report it
        update_status(problem->info, QP_NON_CVX);
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
            update_status(problem->info, QP_SOLVED_INACCURATE);
        } else {
            update_status(problem->info, QP_SOLVED);
        }
        exitflag = 1;
    }
    else if (prim_inf_check) {
        // Update final information
        if (approximate) {
            update_status(problem->info, QP_PRIMAL_INFEASIBLE_INACCURATE);
        } else {
            update_status(problem->info, QP_PRIMAL_INFEASIBLE);
        }

        if (problem->hparams->scaling && !problem->hparams->scaled_termination) {
            // Update infeasibility certificate
            vec_ew_prod(problem->scaling->E, problem->delta_y, problem->delta_y, problem->data->m);
        }
        problem->info->obj_val = INFTY;
        exitflag            = 1;
    }
    else if (dual_inf_check) {
        // Update final information
        if (approximate) {
            update_status(problem->info, QP_DUAL_INFEASIBLE_INACCURATE);
        } else {
            update_status(problem->info, QP_DUAL_INFEASIBLE);
        }

        if (problem->hparams->scaling && !problem->hparams->scaled_termination) {
            // Update infeasibility certificate
            vec_ew_prod(problem->scaling->D, problem->delta_x, problem->delta_x, problem->data->n);
        }
        problem->info->obj_val = -INFTY;
        exitflag            = 1;
    }

    return exitflag;
}
























/**********************
* Main API Functions *
**********************/
void qp_set_default_hparams(qpHyperparams *hparams) {

  hparams->rho           = (float)RHO;            /* ADMM step */
  hparams->sigma         = (float)SIGMA;          /* ADMM step */
  hparams->scaling = SCALING;                       /* heuristic problem scaling */
#if EMBEDDED != 1
  hparams->adaptive_rho           = ADAPTIVE_RHO;
  hparams->adaptive_rho_interval  = ADAPTIVE_RHO_INTERVAL;
  hparams->adaptive_rho_tolerance = (float)ADAPTIVE_RHO_TOLERANCE;

# ifdef PROFILING
  hparams->adaptive_rho_fraction = (float)ADAPTIVE_RHO_FRACTION;
# endif /* ifdef PROFILING */
#endif  /* if EMBEDDED != 1 */

  hparams->max_iter      = MAX_ITER;                /* maximum iterations to
                                                        take */
  hparams->eps_abs       = (float)EPS_ABS;        /* absolute convergence
                                                        tolerance */
  hparams->eps_rel       = (float)EPS_REL;        /* relative convergence
                                                        tolerance */
  hparams->eps_prim_inf  = (float)EPS_PRIM_INF;   /* primal infeasibility
                                                        tolerance */
  hparams->eps_dual_inf  = (float)EPS_DUAL_INF;   /* dual infeasibility
                                                        tolerance */
  hparams->alpha         = (float)ALPHA;          /* relaxation parameter */
  hparams->linsys_solver = LINSYS_SOLVER;           /* relaxation parameter */

#ifndef EMBEDDED
  hparams->delta              = DELTA;              /* regularization parameter
                                                        for polish */
  hparams->polish             = POLISH;             /* ADMM solution polish: 1
                                                      */
  hparams->polish_refine_iter = POLISH_REFINE_ITER; /* iterative refinement
                                                        steps in polish */
  hparams->verbose            = VERBOSE;            /* print output */
#endif /* ifndef EMBEDDED */

  hparams->scaled_termination = SCALED_TERMINATION; /* Evaluate scaled
                                                        termination criteria*/
  hparams->check_termination  = CHECK_TERMINATION;  /* Interval for evaluating
                                                        termination criteria */
  hparams->warm_start         = WARM_START;         /* warm starting */

#ifdef PROFILING
  hparams->time_limit = TIME_LIMIT;
#endif /* ifdef PROFILING */
}

#ifndef EMBEDDED


int qp_setup(qpInstance** problem, const qpData *data, const qpHyperparams *hparams) {
  int exitflag;

  qpInstance * problem;

  // Validate data
  if (validate_data(data)) return qp_error(QP_DATA_VALIDATION_ERROR);

  // Validate hparams
  if (validate_hparams(hparams)) return qp_error(QP_SETTINGS_VALIDATION_ERROR);

  // Allocate empty problemspace
  problem = c_calloc(1, sizeof(qpInstance));
  if (!(problem)) return qp_error(QP_MEM_ALLOC_ERROR);
  *problemp = problem;

  // Start and allocate directly timer
# ifdef PROFILING
  problem->timer = c_malloc(sizeof(QPTimer));
  if (!(problem->timer)) return QPqp_error(QPqpQP_MEM_ALLOC_ERROR);
  qp_tic(problem->timer);
# endif /* ifdef PROFILING */

  // Copy problem data into problemspace
  problem->data = c_malloc(sizeof(qpData));
  if (!(problem->data)) return qp_error(qp_MEM_ALLOC_ERROR);
  problem->data->n = data->n;
  problem->data->m = data->m;

  // Cost function
  problem->data->P = copy_csc_mat(data->P);
  problem->data->q = vec_copy(data->q, data->n);
  if (!(problem->data->P) || !(problem->data->q)) return qp_error(qp_MEM_ALLOC_ERROR);

  // Constraints
  problem->data->A = copy_csc_mat(data->A);
  if (!(problem->data->A)) return qp_error(QP_MEM_ALLOC_ERROR);
  problem->data->l = vec_copy(data->l, data->m);
  problem->data->u = vec_copy(data->u, data->m);
  if ( data->m && (!(problem->data->l) || !(problem->data->u)) )
    return qp_error(QP_MEM_ALLOC_ERROR);

  // Vectorized rho parameter
  problem->rho_vec     = c_malloc(data->m * sizeof(float));
  problem->rho_inv_vec = c_malloc(data->m * sizeof(float));
  if ( data->m && (!(problem->rho_vec) || !(problem->rho_inv_vec)) )
    return qp_error(QP_MEM_ALLOC_ERROR);

  // Type of constraints
  problem->constr_type = c_calloc(data->m, sizeof(int));
  if (data->m && !(problem->constr_type)) return qp_error(QP_MEM_ALLOC_ERROR);

  // Allocate internal solver variables (ADMM steps)
  problem->x        = c_calloc(data->n, sizeof(float));
  problem->z        = c_calloc(data->m, sizeof(float));
  problem->xz_tilde = c_calloc(data->n + data->m, sizeof(float));
  problem->x_prev   = c_calloc(data->n, sizeof(float));
  problem->z_prev   = c_calloc(data->m, sizeof(float));
  problem->y        = c_calloc(data->m, sizeof(float));
  if (!(problem->x) || !(problem->xz_tilde) || !(problem->x_prev))
    return qp_error(QP_MEM_ALLOC_ERROR);
  if ( data->m && (!(problem->z) || !(problem->z_prev) || !(problem->y)) )
    return qp_error(QP_MEM_ALLOC_ERROR);

  // Initialize variables x, y, z to 0
  cold_start(problem);

  // Primal and dual residuals variables
  problem->Ax  = c_calloc(data->m, sizeof(float));
  problem->Px  = c_calloc(data->n, sizeof(float));
  problem->Aty = c_calloc(data->n, sizeof(float));

  // Primal infeasibility variables
  problem->delta_y   = c_calloc(data->m, sizeof(float));
  problem->Atdelta_y = c_calloc(data->n, sizeof(float));

  // Dual infeasibility variables
  problem->delta_x  = c_calloc(data->n, sizeof(float));
  problem->Pdelta_x = c_calloc(data->n, sizeof(float));
  problem->Adelta_x = c_calloc(data->m, sizeof(float));

  if (!(problem->Px) || !(problem->Aty) || !(problem->Atdelta_y) ||
      !(problem->delta_x) || !(problem->Pdelta_x))
    return qp_error(QP_MEM_ALLOC_ERROR);
  if ( data->m && (!(problem->Ax) || !(problem->delta_y) || !(problem->Adelta_x)) )
    return qp_error(QP_MEM_ALLOC_ERROR);

  // Copy hparams
  problem->hparams = copy_hparams(hparams);
  if (!(problem->hparams)) return qp_error(QP_MEM_ALLOC_ERROR);

  // Perform scaling
  if (hparams->scaling) {
    // Allocate scaling structure
    problem->scaling = c_malloc(sizeof(QPScaling));
    if (!(problem->scaling)) return qp_error(QP_MEM_ALLOC_ERROR);
    problem->scaling->D    = c_malloc(data->n * sizeof(float));
    problem->scaling->Dinv = c_malloc(data->n * sizeof(float));
    problem->scaling->E    = c_malloc(data->m * sizeof(float));
    problem->scaling->Einv = c_malloc(data->m * sizeof(float));
    if (!(problem->scaling->D) || !(problem->scaling->Dinv))
      return qp_error(QP_MEM_ALLOC_ERROR);
    if ( data->m && (!(problem->scaling->E) || !(problem->scaling->Einv)) )
      return qp_error(QP_MEM_ALLOC_ERROR);


    // Allocate problemspace variables used in scaling
    problem->D_temp   = c_malloc(data->n * sizeof(float));
    problem->D_temp_A = c_malloc(data->n * sizeof(float));
    problem->E_temp   = c_malloc(data->m * sizeof(float));
    // if (!(problem->D_temp) || !(problem->D_temp_A) || !(problem->E_temp))
    //   return qp_error(QP_MEM_ALLOC_ERROR);
    if (!(problem->D_temp) || !(problem->D_temp_A)) return qp_error(QP_MEM_ALLOC_ERROR);
    if (data->m && !(problem->E_temp))           return qp_error(QP_MEM_ALLOC_ERROR);

    // Scale data
    scale_data(problem);
  } else {
    problem->scaling  = QP_NULL;
    problem->D_temp   = QP_NULL;
    problem->D_temp_A = QP_NULL;
    problem->E_temp   = QP_NULL;
  }

  // Set type of constraints
  set_rho_vec(problem);

  // Load linear system solver
  if (load_linsys_solver(problem->hparams->linsys_solver)) return qp_error(QP_LINSYS_SOLVER_LOAD_ERROR);

  // Initialize linear system solver structure
  exitflag = init_linsys_solver(&(problem->linsys_solver), problem->data->P, problem->data->A,
                                problem->hparams->sigma, problem->rho_vec,
                                problem->hparams->linsys_solver, 0);

  if (exitflag) {
    return qp_error(exitflag);
  }

  // Initialize active constraints structure
  problem->pol = c_malloc(sizeof(QPPolish));
  if (!(problem->pol)) return qp_error(QP_MEM_ALLOC_ERROR);
  problem->pol->Alow_to_A = c_malloc(data->m * sizeof(int));
  problem->pol->Aupp_to_A = c_malloc(data->m * sizeof(int));
  problem->pol->A_to_Alow = c_malloc(data->m * sizeof(int));
  problem->pol->A_to_Aupp = c_malloc(data->m * sizeof(int));
  problem->pol->x         = c_malloc(data->n * sizeof(float));
  problem->pol->z         = c_malloc(data->m * sizeof(float));
  problem->pol->y         = c_malloc(data->m * sizeof(float));
  if (!(problem->pol->x)) return qp_error(QP_MEM_ALLOC_ERROR);
  if ( data->m && (!(problem->pol->Alow_to_A) || !(problem->pol->Aupp_to_A) ||
      !(problem->pol->A_to_Alow) || !(problem->pol->A_to_Aupp) ||
      !(problem->pol->z) || !(problem->pol->y)) )
    return qp_error(QP_MEM_ALLOC_ERROR);

  // Allocate solution
  problem->solution = c_calloc(1, sizeof(QPSolution));
  if (!(problem->solution)) return qp_error(QP_MEM_ALLOC_ERROR);
  problem->solution->x = c_calloc(1, data->n * sizeof(float));
  problem->solution->y = c_calloc(1, data->m * sizeof(float));
  if (!(problem->solution->x))            return qp_error(QP_MEM_ALLOC_ERROR);
  if (data->m && !(problem->solution->y)) return qp_error(QP_MEM_ALLOC_ERROR);

  // Allocate and initialize information
  problem->info = c_calloc(1, sizeof(QPInfo));
  if (!(problem->info)) return qp_error(QP_MEM_ALLOC_ERROR);
  problem->info->status_polish = 0;              // Polishing not performed
  update_status(problem->info, QP_UNSOLVED);
# ifdef PROFILING
  problem->info->solve_time  = 0.0;                   // Solve time to zero
  problem->info->update_time = 0.0;                   // Update time to zero
  problem->info->polish_time = 0.0;                   // Polish time to zero
  problem->info->run_time    = 0.0;                   // Total run time to zero
  problem->info->setup_time  = qp_toc(problem->timer); // Update timer information

  problem->first_run         = 1;
  problem->clear_update_time = 0;
  problem->rho_update_from_solve = 0;
# endif /* ifdef PROFILING */
  problem->info->rho_updates  = 0;                    // Rho updates set to 0
  problem->info->rho_estimate = problem->hparams->rho;  // Best rho estimate

  // Print header
# ifdef PRINTING
  if (problem->hparams->verbose) print_setup_header(problem);
  problem->summary_printed = 0; // Initialize last summary  to not printed
# endif /* ifdef PRINTING */


  // If adaptive rho and automatic interval, but profiling disabled, we need to
  // set the interval to a default value
# ifndef PROFILING
  if (problem->hparams->adaptive_rho && !problem->hparams->adaptive_rho_interval) {
    if (problem->hparams->check_termination) {
      // If check_termination is enabled, we set it to a multiple of the check
      // termination interval
      problem->hparams->adaptive_rho_interval = ADAPTIVE_RHO_MULTIPLE_TERMINATION *
                                              problem->hparams->check_termination;
    } else {
      // If check_termination is disabled we set it to a predefined fix number
      problem->hparams->adaptive_rho_interval = ADAPTIVE_RHO_FIXED;
    }
  }
# endif /* ifndef PROFILING */

  // Return exit flag
  return 0;
}

#endif // #ifndef EMBEDDED


int qp_solve(qpInstance *problem) {

  int exitflag;
  int iter;
  int compute_cost_function; // Boolean: compute the cost function in the loop or not
  int can_check_termination; // Boolean: check termination or not

#ifdef PROFILING
  float temp_run_time;       // Temporary variable to store current run time
#endif /* ifdef PROFILING */

#ifdef PRINTING
  int can_print;             // Boolean whether you can print
#endif /* ifdef PRINTING */

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

#ifdef PROFILING
  if (problem->clear_update_time == 1)
    problem->info->update_time = 0.0;
  problem->rho_update_from_solve = 1;
#endif /* ifdef PROFILING */

  // Initialize variables
  exitflag              = 0;
  can_check_termination = 0;
#ifdef PRINTING
  can_print = problem->hparams->verbose;
#endif /* ifdef PRINTING */
#ifdef PRINTING
  compute_cost_function = problem->hparams->verbose; // Compute cost function only
                                                   // if verbose is on
#else /* ifdef PRINTING */
  compute_cost_function = 0;                       // Never compute cost
                                                   // function during the
                                                   // iterations if no printing
                                                   // enabled
#endif /* ifdef PRINTING */



#ifdef PROFILING
  qp_tic(problem->timer); // Start timer
#endif /* ifdef PROFILING */


#ifdef PRINTING

  if (problem->hparams->verbose) {
    // Print Header for every column
    print_header();
  }
#endif /* ifdef PRINTING */

#ifdef CTRLC

  // initialize Ctrl-C support
  qp_start_interrupt_listener();
#endif /* ifdef CTRLC */

  // Initialize variables (cold start or warm start depending on hparams)
  if (!problem->hparams->warm_start) cold_start(problem);  // If not warm start ->
                                                      // set x, z, y to zero

  // Main ADMM algorithm
  for (iter = 1; iter <= problem->hparams->max_iter; iter++) {
    // Update x_prev, z_prev (preallocated, no malloc)
    swap_vectors(&(problem->x), &(problem->x_prev));
    swap_vectors(&(problem->z), &(problem->z_prev));

    /* ADMM STEPS */
    /* Compute \tilde{x}^{k+1}, \tilde{z}^{k+1} */
    update_xz_tilde(problem);

    /* Compute x^{k+1} */
    update_x(problem);

    /* Compute z^{k+1} */
    update_z(problem);

    /* Compute y^{k+1} */
    update_y(problem);

    /* End of ADMM Steps */

#ifdef CTRLC

    // Check the interrupt signal
    if (qp_is_interrupted()) {
      update_status(problem->info, QP_SIGINT);
# ifdef PRINTING
      c_print("Solver interrupted\n");
# endif /* ifdef PRINTING */
      exitflag = 1;
      goto exit;
    }
#endif /* ifdef CTRLC */

#ifdef PROFILING

    // Check if solver time_limit is enabled. In case, check if the current
    // run time is more than the time_limit option.
    if (problem->first_run) {
      temp_run_time = problem->info->setup_time + qp_toc(problem->timer);
    }
    else {
      temp_run_time = problem->info->update_time + qp_toc(problem->timer);
    }

    if (problem->hparams->time_limit &&
        (temp_run_time >= problem->hparams->time_limit)) {
      update_status(problem->info, QP_TIME_LIMIT_REACHED);
# ifdef PRINTING
      if (problem->hparams->verbose) c_print("run time limit reached\n");
      can_print = 0;  // Not printing at this iteration
# endif /* ifdef PRINTING */
      break;
    }
#endif /* ifdef PROFILING */


    // Can we check for termination ?
    can_check_termination = problem->hparams->check_termination &&
                            (iter % problem->hparams->check_termination == 0);

#ifdef PRINTING

    // Can we print ?
    can_print = problem->hparams->verbose &&
                ((iter % PRINT_INTERVAL == 0) || (iter == 1));

    if (can_check_termination || can_print) { // Update status in either of
                                              // these cases
      // Update information
      update_info(problem, iter, compute_cost_function, 0);

      if (can_print) {
        // Print summary
        print_summary(problem);
      }

      if (can_check_termination) {
        // Check algorithm termination
        if (check_termination(problem, 0)) {
          // Terminate algorithm
          break;
        }
      }
    }
#else /* ifdef PRINTING */

    if (can_check_termination) {
      // Update information and compute also objective value
      update_info(problem, iter, compute_cost_function, 0);

      // Check algorithm termination
      if (check_termination(problem, 0)) {
        // Terminate algorithm
        break;
      }
    }
#endif /* ifdef PRINTING */


#if EMBEDDED != 1
# ifdef PROFILING

    // If adaptive rho with automatic interval, check if the solve time is a
    // certain fraction
    // of the setup time.
    if (problem->hparams->adaptive_rho && !problem->hparams->adaptive_rho_interval) {
      // Check time
      if (qp_toc(problem->timer) >
          problem->hparams->adaptive_rho_fraction * problem->info->setup_time) {
        // Enough time has passed. We now get the number of iterations between
        // the updates.
        if (problem->hparams->check_termination) {
          // If check_termination is enabled, we round the number of iterations
          // between
          // rho updates to the closest multiple of check_termination
          problem->hparams->adaptive_rho_interval = (int)c_roundmultiple(iter,
                                                                         problem->hparams->check_termination);
        } else {
          // If check_termination is disabled, we round the number of iterations
          // between
          // updates to the closest multiple of the default check_termination
          // interval.
          problem->hparams->adaptive_rho_interval = (int)c_roundmultiple(iter,
                                                                         CHECK_TERMINATION);
        }

        // Make sure the interval is not 0 and at least check_termination times
        problem->hparams->adaptive_rho_interval = c_max(
          problem->hparams->adaptive_rho_interval,
          problem->hparams->check_termination);
      } // If time condition is met
    }   // If adaptive rho enabled and interval set to auto
# endif // PROFILING

    // Adapt rho
    if (problem->hparams->adaptive_rho &&
        problem->hparams->adaptive_rho_interval &&
        (iter % problem->hparams->adaptive_rho_interval == 0)) {
      // Update info with the residuals if it hasn't been done before
# ifdef PRINTING

      if (!can_check_termination && !can_print) {
        // Information has not been computed neither for termination or printing
        // reasons
        update_info(problem, iter, compute_cost_function, 0);
      }
# else /* ifdef PRINTING */

      if (!can_check_termination) {
        // Information has not been computed before for termination check
        update_info(problem, iter, compute_cost_function, 0);
      }
# endif /* ifdef PRINTING */

      // Actually update rho
      if (adapt_rho(problem)) {
# ifdef PRINTING
        c_eprint("Failed rho update");
# endif // PRINTING
        exitflag = 1;
        goto exit;
      }
    }
#endif // EMBEDDED != 1

  }        // End of ADMM for loop


  // Update information and check termination condition if it hasn't been done
  // during last iteration (max_iter reached or check_termination disabled)
  if (!can_check_termination) {
    /* Update information */
#ifdef PRINTING

    if (!can_print) {
      // Update info only if it hasn't been updated before for printing
      // reasons
      update_info(problem, iter - 1, compute_cost_function, 0);
    }
#else /* ifdef PRINTING */

    // If no printing is enabled, update info directly
    update_info(problem, iter - 1, compute_cost_function, 0);
#endif /* ifdef PRINTING */

#ifdef PRINTING

    /* Print summary */
    if (problem->hparams->verbose && !problem->summary_printed) print_summary(problem);
#endif /* ifdef PRINTING */

    /* Check whether a termination criterion is triggered */
    check_termination(problem, 0);
  }

  // Compute objective value in case it was not
  // computed during the iterations
  if (!compute_cost_function && has_solution(problem->info)){
    problem->info->obj_val = compute_obj_val(problem, problem->x);
  }


#ifdef PRINTING
  /* Print summary for last iteration */
  if (problem->hparams->verbose && !problem->summary_printed) {
    print_summary(problem);
  }
#endif /* ifdef PRINTING */

  /* if max iterations reached, change status accordingly */
  if (problem->info->status_val == QP_UNSOLVED) {
    if (!check_termination(problem, 1)) { // Try to check for approximate
      update_status(problem->info, QP_MAX_ITER_REACHED);
    }
  }

#ifdef PROFILING
  /* if time-limit reached check termination and update status accordingly */
 if (problem->info->status_val == QP_TIME_LIMIT_REACHED) {
    if (!check_termination(problem, 1)) { // Try for approximate solutions
      update_status(problem->info, QP_TIME_LIMIT_REACHED); /* Change update status back to QP_TIME_LIMIT_REACHED */
    }
  }
#endif /* ifdef PROFILING */


#if EMBEDDED != 1
  /* Update rho estimate */
  problem->info->rho_estimate = compute_rho_estimate(problem);
#endif /* if EMBEDDED != 1 */

  /* Update solve time */
#ifdef PROFILING
  problem->info->solve_time = qp_toc(problem->timer);
#endif /* ifdef PROFILING */


#ifndef EMBEDDED
  // Polish the obtained solution
  if (problem->hparams->polish && (problem->info->status_val == QP_SOLVED))
    polish(problem);
#endif /* ifndef EMBEDDED */

#ifdef PROFILING
  /* Update total time */
  if (problem->first_run) {
    // total time: setup + solve + polish
    problem->info->run_time = problem->info->setup_time +
                           problem->info->solve_time +
                           problem->info->polish_time;
  } else {
    // total time: update + solve + polish
    problem->info->run_time = problem->info->update_time +
                           problem->info->solve_time +
                           problem->info->polish_time;
  }

  // Indicate that the solve function has already been executed
  if (problem->first_run) problem->first_run = 0;

  // Indicate that the update_time should be set to zero
  problem->clear_update_time = 1;

  // Indicate that qp_update_rho is not called from qp_solve
  problem->rho_update_from_solve = 0;
#endif /* ifdef PROFILING */

#ifdef PRINTING
  /* Print final footer */
  if (problem->hparams->verbose) print_footer(problem->info, problem->hparams->polish);
#endif /* ifdef PRINTING */

  // Store solution
  store_solution(problem);


// Define exit flag for quitting function
#if defined(PROFILING) || defined(CTRLC) || EMBEDDED != 1
exit:
#endif /* if defined(PROFILING) || defined(CTRLC) || EMBEDDED != 1 */

#ifdef CTRLC
  // Restore previous signal handler
  qp_end_interrupt_listener();
#endif /* ifdef CTRLC */

  return exitflag;
}


#ifndef EMBEDDED

int qp_cleanup(qpInstance *problem) {
  int exitflag = 0;

  if (problem) { // If problemspace has been allocated
    // Free Data
    if (problem->data) {
      if (problem->data->P) csc_spfree(problem->data->P);
      if (problem->data->A) csc_spfree(problem->data->A);
      if (problem->data->q) c_free(problem->data->q);
      if (problem->data->l) c_free(problem->data->l);
      if (problem->data->u) c_free(problem->data->u);
      c_free(problem->data);
    }

    // Free scaling variables
    if (problem->scaling){
      if (problem->scaling->D)    c_free(problem->scaling->D);
      if (problem->scaling->Dinv) c_free(problem->scaling->Dinv);
      if (problem->scaling->E)    c_free(problem->scaling->E);
      if (problem->scaling->Einv) c_free(problem->scaling->Einv);
      c_free(problem->scaling);
    }

    // Free temp problemspace variables for scaling
    if (problem->D_temp)   c_free(problem->D_temp);
    if (problem->D_temp_A) c_free(problem->D_temp_A);
    if (problem->E_temp)   c_free(problem->E_temp);

    // Free linear system solver structure
    if (problem->linsys_solver) {
      if (problem->linsys_solver->free) {
        problem->linsys_solver->free(problem->linsys_solver);
      }
    }

    // Unload linear system solver after free
    if (problem->hparams) {
      exitflag = unload_linsys_solver(problem->hparams->linsys_solver);
    }

#ifndef EMBEDDED
    // Free active constraints structure
    if (problem->pol) {
      if (problem->pol->Alow_to_A) c_free(problem->pol->Alow_to_A);
      if (problem->pol->Aupp_to_A) c_free(problem->pol->Aupp_to_A);
      if (problem->pol->A_to_Alow) c_free(problem->pol->A_to_Alow);
      if (problem->pol->A_to_Aupp) c_free(problem->pol->A_to_Aupp);
      if (problem->pol->x)         c_free(problem->pol->x);
      if (problem->pol->z)         c_free(problem->pol->z);
      if (problem->pol->y)         c_free(problem->pol->y);
      c_free(problem->pol);
    }
#endif /* ifndef EMBEDDED */

    // Free other Variables
    if (problem->rho_vec)     c_free(problem->rho_vec);
    if (problem->rho_inv_vec) c_free(problem->rho_inv_vec);
#if EMBEDDED != 1
    if (problem->constr_type) c_free(problem->constr_type);
#endif
    if (problem->x)           c_free(problem->x);
    if (problem->z)           c_free(problem->z);
    if (problem->xz_tilde)    c_free(problem->xz_tilde);
    if (problem->x_prev)      c_free(problem->x_prev);
    if (problem->z_prev)      c_free(problem->z_prev);
    if (problem->y)           c_free(problem->y);
    if (problem->Ax)          c_free(problem->Ax);
    if (problem->Px)          c_free(problem->Px);
    if (problem->Aty)         c_free(problem->Aty);
    if (problem->delta_y)     c_free(problem->delta_y);
    if (problem->Atdelta_y)   c_free(problem->Atdelta_y);
    if (problem->delta_x)     c_free(problem->delta_x);
    if (problem->Pdelta_x)    c_free(problem->Pdelta_x);
    if (problem->Adelta_x)    c_free(problem->Adelta_x);

    // Free Settings
    if (problem->hparams) c_free(problem->hparams);

    // Free solution
    if (problem->solution) {
      if (problem->solution->x) c_free(problem->solution->x);
      if (problem->solution->y) c_free(problem->solution->y);
      c_free(problem->solution);
    }

    // Free information
    if (problem->info) c_free(problem->info);

# ifdef PROFILING
    // Free timer
    if (problem->timer) c_free(problem->timer);
# endif /* ifdef PROFILING */

    // Free problem
    c_free(problem);
  }

  return exitflag;
}

#endif // #ifndef EMBEDDED


/************************
* Update problem data  *
************************/
int qp_update_lin_cost(qpInstance *problem, const float *q_new) {

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

#ifdef PROFILING
  if (problem->clear_update_time == 1) {
    problem->clear_update_time = 0;
    problem->info->update_time = 0.0;
  }
  qp_tic(problem->timer); // Start timer
#endif /* ifdef PROFILING */

  // Replace q by the new vector
  prea_vec_copy(q_new, problem->data->q, problem->data->n);

  // Scaling
  if (problem->hparams->scaling) {
    vec_ew_prod(problem->scaling->D, problem->data->q, problem->data->q, problem->data->n);
    vec_mult_scalar(problem->data->q, problem->scaling->c, problem->data->n);
  }

  // Reset solver information
  reset_info(problem->info);

#ifdef PROFILING
  problem->info->update_time += qp_toc(problem->timer);
#endif /* ifdef PROFILING */

  return 0;
}

int qp_update_bounds(qpInstance *problem,
                         const float *l_new,
                         const float *u_new) {
  int i, exitflag = 0;

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

#ifdef PROFILING
  if (problem->clear_update_time == 1) {
    problem->clear_update_time = 0;
    problem->info->update_time = 0.0;
  }
  qp_tic(problem->timer); // Start timer
#endif /* ifdef PROFILING */

  // Check if lower bound is smaller than upper bound
  for (i = 0; i < problem->data->m; i++) {
    if (l_new[i] > u_new[i]) {
#ifdef PRINTING
      c_eprint("lower bound must be lower than or equal to upper bound");
#endif /* ifdef PRINTING */
      return 1;
    }
  }

  // Replace l and u by the new vectors
  prea_vec_copy(l_new, problem->data->l, problem->data->m);
  prea_vec_copy(u_new, problem->data->u, problem->data->m);

  // Scaling
  if (problem->hparams->scaling) {
    vec_ew_prod(problem->scaling->E, problem->data->l, problem->data->l, problem->data->m);
    vec_ew_prod(problem->scaling->E, problem->data->u, problem->data->u, problem->data->m);
  }

  // Reset solver information
  reset_info(problem->info);

#if EMBEDDED != 1
  // Update rho_vec and refactor if constraints type changes
  exitflag = update_rho_vec(problem);
#endif // EMBEDDED != 1

#ifdef PROFILING
  problem->info->update_time += qp_toc(problem->timer);
#endif /* ifdef PROFILING */

  return exitflag;
}

int qp_update_lower_bound(qpInstance *problem, const float *l_new) {
  int i, exitflag = 0;

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

#ifdef PROFILING
  if (problem->clear_update_time == 1) {
    problem->clear_update_time = 0;
    problem->info->update_time = 0.0;
  }
  qp_tic(problem->timer); // Start timer
#endif /* ifdef PROFILING */

  // Replace l by the new vector
  prea_vec_copy(l_new, problem->data->l, problem->data->m);

  // Scaling
  if (problem->hparams->scaling) {
    vec_ew_prod(problem->scaling->E, problem->data->l, problem->data->l, problem->data->m);
  }

  // Check if lower bound is smaller than upper bound
  for (i = 0; i < problem->data->m; i++) {
    if (problem->data->l[i] > problem->data->u[i]) {
#ifdef PRINTING
      c_eprint("upper bound must be greater than or equal to lower bound");
#endif /* ifdef PRINTING */
      return 1;
    }
  }

  // Reset solver information
  reset_info(problem->info);

#if EMBEDDED != 1
  // Update rho_vec and refactor if constraints type changes
  exitflag = update_rho_vec(problem);
#endif // EMBEDDED ! =1

#ifdef PROFILING
  problem->info->update_time += qp_toc(problem->timer);
#endif /* ifdef PROFILING */

  return exitflag;
}

int qp_update_upper_bound(qpInstance *problem, const float *u_new) {
  int i, exitflag = 0;

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

#ifdef PROFILING
  if (problem->clear_update_time == 1) {
    problem->clear_update_time = 0;
    problem->info->update_time = 0.0;
  }
  qp_tic(problem->timer); // Start timer
#endif /* ifdef PROFILING */

  // Replace u by the new vector
  prea_vec_copy(u_new, problem->data->u, problem->data->m);

  // Scaling
  if (problem->hparams->scaling) {
    vec_ew_prod(problem->scaling->E, problem->data->u, problem->data->u, problem->data->m);
  }

  // Check if upper bound is greater than lower bound
  for (i = 0; i < problem->data->m; i++) {
    if (problem->data->u[i] < problem->data->l[i]) {
#ifdef PRINTING
      c_eprint("lower bound must be lower than or equal to upper bound");
#endif /* ifdef PRINTING */
      return 1;
    }
  }

  // Reset solver information
  reset_info(problem->info);

#if EMBEDDED != 1
  // Update rho_vec and refactor if constraints type changes
  exitflag = update_rho_vec(problem);
#endif // EMBEDDED != 1

#ifdef PROFILING
  problem->info->update_time += qp_toc(problem->timer);
#endif /* ifdef PROFILING */

  return exitflag;
}

int qp_warm_start(qpInstance *problem, const float *x, const float *y) {

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

  // Update warm_start setting to true
  if (!problem->hparams->warm_start) problem->hparams->warm_start = 1;

  // Copy primal and dual variables into the iterates
  prea_vec_copy(x, problem->x, problem->data->n);
  prea_vec_copy(y, problem->y, problem->data->m);

  // Scale iterates
  if (problem->hparams->scaling) {
    vec_ew_prod(problem->scaling->Dinv, problem->x, problem->x, problem->data->n);
    vec_ew_prod(problem->scaling->Einv, problem->y, problem->y, problem->data->m);
    vec_mult_scalar(problem->y, problem->scaling->c, problem->data->m);
  }

  // Compute Ax = z and store it in z
  mat_vec(problem->data->A, problem->x, problem->z, 0);

  return 0;
}

int qp_warm_start_x(qpInstance *problem, const float *x) {

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

  // Update warm_start setting to true
  if (!problem->hparams->warm_start) problem->hparams->warm_start = 1;

  // Copy primal variable into the iterate x
  prea_vec_copy(x, problem->x, problem->data->n);

  // Scale iterate
  if (problem->hparams->scaling) {
    vec_ew_prod(problem->scaling->Dinv, problem->x, problem->x, problem->data->n);
  }

  // Compute Ax = z and store it in z
  mat_vec(problem->data->A, problem->x, problem->z, 0);

  return 0;
}

int qp_warm_start_y(qpInstance *problem, const float *y) {

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

  // Update warm_start setting to true
  if (!problem->hparams->warm_start) problem->hparams->warm_start = 1;

  // Copy primal variable into the iterate y
  prea_vec_copy(y, problem->y, problem->data->m);

  // Scale iterate
  if (problem->hparams->scaling) {
    vec_ew_prod(problem->scaling->Einv, problem->y, problem->y, problem->data->m);
    vec_mult_scalar(problem->y, problem->scaling->c, problem->data->m);
  }

  return 0;
}


#if EMBEDDED != 1

int qp_update_P(qpInstance *problem,
                    const float *Px_new,
                    const int   *Px_new_idx,
                    int          P_new_n) {
  int i;        // For indexing
  int exitflag; // Exit flag
  int nnzP;     // Number of nonzeros in P

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

#ifdef PROFILING
  if (problem->clear_update_time == 1) {
    problem->clear_update_time = 0;
    problem->info->update_time = 0.0;
  }
  qp_tic(problem->timer); // Start timer
#endif /* ifdef PROFILING */

  nnzP = problem->data->P->p[problem->data->P->n];

  if (Px_new_idx) { // Passing the index of elements changed
    // Check if number of elements is less or equal than the total number of
    // nonzeros in P
    if (P_new_n > nnzP) {
# ifdef PRINTING
      c_eprint("new number of elements (%i) greater than elements in P (%i)",
               (int)P_new_n,
               (int)nnzP);
# endif /* ifdef PRINTING */
      return 1;
    }
  }

  if (problem->hparams->scaling) {
    // Unscale data
    unscale_data(problem);
  }

  // Update P elements
  if (Px_new_idx) { // Change only Px_new_idx
    for (i = 0; i < P_new_n; i++) {
      problem->data->P->x[Px_new_idx[i]] = Px_new[i];
    }
  }
  else // Change whole P
  {
    for (i = 0; i < nnzP; i++) {
      problem->data->P->x[i] = Px_new[i];
    }
  }

  if (problem->hparams->scaling) {
    // Scale data
    scale_data(problem);
  }

  // Update linear system structure with new data
  exitflag = problem->linsys_solver->update_matrices(problem->linsys_solver,
                                                  problem->data->P,
                                                  problem->data->A);

  // Reset solver information
  reset_info(problem->info);

# ifdef PRINTING

  if (exitflag < 0) {
    c_eprint("new KKT matrix is not quasidefinite");
  }
# endif /* ifdef PRINTING */

#ifdef PROFILING
  problem->info->update_time += qp_toc(problem->timer);
#endif /* ifdef PROFILING */

  return exitflag;
}


int qp_update_A(qpInstance *problem,
                    const float *Ax_new,
                    const int   *Ax_new_idx,
                    int          A_new_n) {
  int i;        // For indexing
  int exitflag; // Exit flag
  int nnzA;     // Number of nonzeros in A

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

#ifdef PROFILING
  if (problem->clear_update_time == 1) {
    problem->clear_update_time = 0;
    problem->info->update_time = 0.0;
  }
  qp_tic(problem->timer); // Start timer
#endif /* ifdef PROFILING */

  nnzA = problem->data->A->p[problem->data->A->n];

  if (Ax_new_idx) { // Passing the index of elements changed
    // Check if number of elements is less or equal than the total number of
    // nonzeros in A
    if (A_new_n > nnzA) {
# ifdef PRINTING
      c_eprint("new number of elements (%i) greater than elements in A (%i)",
               (int)A_new_n,
               (int)nnzA);
# endif /* ifdef PRINTING */
      return 1;
    }
  }

  if (problem->hparams->scaling) {
    // Unscale data
    unscale_data(problem);
  }

  // Update A elements
  if (Ax_new_idx) { // Change only Ax_new_idx
    for (i = 0; i < A_new_n; i++) {
      problem->data->A->x[Ax_new_idx[i]] = Ax_new[i];
    }
  }
  else { // Change whole A
    for (i = 0; i < nnzA; i++) {
      problem->data->A->x[i] = Ax_new[i];
    }
  }

  if (problem->hparams->scaling) {
    // Scale data
    scale_data(problem);
  }

  // Update linear system structure with new data
  exitflag = problem->linsys_solver->update_matrices(problem->linsys_solver,
                                                  problem->data->P,
                                                  problem->data->A);

  // Reset solver information
  reset_info(problem->info);

# ifdef PRINTING

  if (exitflag < 0) {
    c_eprint("new KKT matrix is not quasidefinite");
  }
# endif /* ifdef PRINTING */

#ifdef PROFILING
  problem->info->update_time += qp_toc(problem->timer);
#endif /* ifdef PROFILING */

  return exitflag;
}


int qp_update_P_A(qpInstance *problem,
                      const float *Px_new,
                      const int   *Px_new_idx,
                      int          P_new_n,
                      const float *Ax_new,
                      const int   *Ax_new_idx,
                      int          A_new_n) {
  int i;          // For indexing
  int exitflag;   // Exit flag
  int nnzP, nnzA; // Number of nonzeros in P and A

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

#ifdef PROFILING
  if (problem->clear_update_time == 1) {
    problem->clear_update_time = 0;
    problem->info->update_time = 0.0;
  }
  qp_tic(problem->timer); // Start timer
#endif /* ifdef PROFILING */

  nnzP = problem->data->P->p[problem->data->P->n];
  nnzA = problem->data->A->p[problem->data->A->n];


  if (Px_new_idx) { // Passing the index of elements changed
    // Check if number of elements is less or equal than the total number of
    // nonzeros in P
    if (P_new_n > nnzP) {
# ifdef PRINTING
      c_eprint("new number of elements (%i) greater than elements in P (%i)",
               (int)P_new_n,
               (int)nnzP);
# endif /* ifdef PRINTING */
      return 1;
    }
  }


  if (Ax_new_idx) { // Passing the index of elements changed
    // Check if number of elements is less or equal than the total number of
    // nonzeros in A
    if (A_new_n > nnzA) {
# ifdef PRINTING
      c_eprint("new number of elements (%i) greater than elements in A (%i)",
               (int)A_new_n,
               (int)nnzA);
# endif /* ifdef PRINTING */
      return 2;
    }
  }

  if (problem->hparams->scaling) {
    // Unscale data
    unscale_data(problem);
  }

  // Update P elements
  if (Px_new_idx) { // Change only Px_new_idx
    for (i = 0; i < P_new_n; i++) {
      problem->data->P->x[Px_new_idx[i]] = Px_new[i];
    }
  }
  else // Change whole P
  {
    for (i = 0; i < nnzP; i++) {
      problem->data->P->x[i] = Px_new[i];
    }
  }

  // Update A elements
  if (Ax_new_idx) { // Change only Ax_new_idx
    for (i = 0; i < A_new_n; i++) {
      problem->data->A->x[Ax_new_idx[i]] = Ax_new[i];
    }
  }
  else { // Change whole A
    for (i = 0; i < nnzA; i++) {
      problem->data->A->x[i] = Ax_new[i];
    }
  }

  if (problem->hparams->scaling) {
    // Scale data
    scale_data(problem);
  }

  // Update linear system structure with new data
  exitflag = problem->linsys_solver->update_matrices(problem->linsys_solver,
                                                  problem->data->P,
                                                  problem->data->A);

  // Reset solver information
  reset_info(problem->info);

# ifdef PRINTING

  if (exitflag < 0) {
    c_eprint("new KKT matrix is not quasidefinite");
  }
# endif /* ifdef PRINTING */

#ifdef PROFILING
  problem->info->update_time += qp_toc(problem->timer);
#endif /* ifdef PROFILING */

  return exitflag;
}

int qp_update_rho(qpInstance *problem, float rho_new) {
  int exitflag, i;

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

  // Check value of rho
  if (rho_new <= 0) {
# ifdef PRINTING
    c_eprint("rho must be positive");
# endif /* ifdef PRINTING */
    return 1;
  }

#ifdef PROFILING
  if (problem->rho_update_from_solve == 0) {
    if (problem->clear_update_time == 1) {
      problem->clear_update_time = 0;
      problem->info->update_time = 0.0;
    }
    qp_tic(problem->timer); // Start timer
  }
#endif /* ifdef PROFILING */

  // Update rho in hparams
  problem->hparams->rho = c_min(c_max(rho_new, RHO_MIN), RHO_MAX);

  // Update rho_vec and rho_inv_vec
  for (i = 0; i < problem->data->m; i++) {
    if (problem->constr_type[i] == 0) {
      // Inequalities
      problem->rho_vec[i]     = problem->hparams->rho;
      problem->rho_inv_vec[i] = 1. / problem->hparams->rho;
    }
    else if (problem->constr_type[i] == 1) {
      // Equalities
      problem->rho_vec[i]     = RHO_EQ_OVER_RHO_INEQ * problem->hparams->rho;
      problem->rho_inv_vec[i] = 1. / problem->rho_vec[i];
    }
  }

  // Update rho_vec in KKT matrix
  exitflag = problem->linsys_solver->update_rho_vec(problem->linsys_solver,
                                                 problem->rho_vec);

#ifdef PROFILING
  if (problem->rho_update_from_solve == 0)
    problem->info->update_time += qp_toc(problem->timer);
#endif /* ifdef PROFILING */

  return exitflag;
}

#endif // EMBEDDED != 1

/****************************
* Update problem hparams  *
****************************/
int qp_update_max_iter(qpInstance *problem, int max_iter_new) {

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

  // Check that max_iter is positive
  if (max_iter_new <= 0) {
#ifdef PRINTING
    c_eprint("max_iter must be positive");
#endif /* ifdef PRINTING */
    return 1;
  }

  // Update max_iter
  problem->hparams->max_iter = max_iter_new;

  return 0;
}

int qp_update_eps_abs(qpInstance *problem, float eps_abs_new) {

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

  // Check that eps_abs is positive
  if (eps_abs_new < 0.) {
#ifdef PRINTING
    c_eprint("eps_abs must be nonnegative");
#endif /* ifdef PRINTING */
    return 1;
  }

  // Update eps_abs
  problem->hparams->eps_abs = eps_abs_new;

  return 0;
}

int qp_update_eps_rel(qpInstance *problem, float eps_rel_new) {

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

  // Check that eps_rel is positive
  if (eps_rel_new < 0.) {
#ifdef PRINTING
    c_eprint("eps_rel must be nonnegative");
#endif /* ifdef PRINTING */
    return 1;
  }

  // Update eps_rel
  problem->hparams->eps_rel = eps_rel_new;

  return 0;
}

int qp_update_eps_prim_inf(qpInstance *problem, float eps_prim_inf_new) {

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

  // Check that eps_prim_inf is positive
  if (eps_prim_inf_new < 0.) {
#ifdef PRINTING
    c_eprint("eps_prim_inf must be nonnegative");
#endif /* ifdef PRINTING */
    return 1;
  }

  // Update eps_prim_inf
  problem->hparams->eps_prim_inf = eps_prim_inf_new;

  return 0;
}

int qp_update_eps_dual_inf(qpInstance *problem, float eps_dual_inf_new) {

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

  // Check that eps_dual_inf is positive
  if (eps_dual_inf_new < 0.) {
#ifdef PRINTING
    c_eprint("eps_dual_inf must be nonnegative");
#endif /* ifdef PRINTING */
    return 1;
  }

  // Update eps_dual_inf
  problem->hparams->eps_dual_inf = eps_dual_inf_new;


  return 0;
}

int qp_update_alpha(qpInstance *problem, float alpha_new) {

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

  // Check that alpha is between 0 and 2
  if ((alpha_new <= 0.) || (alpha_new >= 2.)) {
#ifdef PRINTING
    c_eprint("alpha must be between 0 and 2");
#endif /* ifdef PRINTING */
    return 1;
  }

  // Update alpha
  problem->hparams->alpha = alpha_new;

  return 0;
}

int qp_update_warm_start(qpInstance *problem, int warm_start_new) {

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

  // Check that warm_start is either 0 or 1
  if ((warm_start_new != 0) && (warm_start_new != 1)) {
#ifdef PRINTING
    c_eprint("warm_start should be either 0 or 1");
#endif /* ifdef PRINTING */
    return 1;
  }

  // Update warm_start
  problem->hparams->warm_start = warm_start_new;

  return 0;
}

int qp_update_scaled_termination(qpInstance *problem, int scaled_termination_new) {

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

  // Check that scaled_termination is either 0 or 1
  if ((scaled_termination_new != 0) && (scaled_termination_new != 1)) {
#ifdef PRINTING
    c_eprint("scaled_termination should be either 0 or 1");
#endif /* ifdef PRINTING */
    return 1;
  }

  // Update scaled_termination
  problem->hparams->scaled_termination = scaled_termination_new;

  return 0;
}

int qp_update_check_termination(qpInstance *problem, int check_termination_new) {

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

  // Check that check_termination is nonnegative
  if (check_termination_new < 0) {
#ifdef PRINTING
    c_eprint("check_termination should be nonnegative");
#endif /* ifdef PRINTING */
    return 1;
  }

  // Update check_termination
  problem->hparams->check_termination = check_termination_new;

  return 0;
}

#ifndef EMBEDDED

int qp_update_delta(qpInstance *problem, float delta_new) {

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

  // Check that delta is positive
  if (delta_new <= 0.) {
# ifdef PRINTING
    c_eprint("delta must be positive");
# endif /* ifdef PRINTING */
    return 1;
  }

  // Update delta
  problem->hparams->delta = delta_new;

  return 0;
}

int qp_update_polish(qpInstance *problem, int polish_new) {

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

  // Check that polish is either 0 or 1
  if ((polish_new != 0) && (polish_new != 1)) {
# ifdef PRINTING
    c_eprint("polish should be either 0 or 1");
# endif /* ifdef PRINTING */
    return 1;
  }

  // Update polish
  problem->hparams->polish = polish_new;

# ifdef PROFILING

  // Reset polish time to zero
  problem->info->polish_time = 0.0;
# endif /* ifdef PROFILING */

  return 0;
}

int qp_update_polish_refine_iter(qpInstance *problem, int polish_refine_iter_new) {

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

  // Check that polish_refine_iter is nonnegative
  if (polish_refine_iter_new < 0) {
# ifdef PRINTING
    c_eprint("polish_refine_iter must be nonnegative");
# endif /* ifdef PRINTING */
    return 1;
  }

  // Update polish_refine_iter
  problem->hparams->polish_refine_iter = polish_refine_iter_new;

  return 0;
}

int qp_update_verbose(qpInstance *problem, int verbose_new) {

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

  // Check that verbose is either 0 or 1
  if ((verbose_new != 0) && (verbose_new != 1)) {
# ifdef PRINTING
    c_eprint("verbose should be either 0 or 1");
# endif /* ifdef PRINTING */
    return 1;
  }

  // Update verbose
  problem->hparams->verbose = verbose_new;

  return 0;
}

#endif // EMBEDDED

#ifdef PROFILING

int qp_update_time_limit(qpInstance *problem, float time_limit_new) {

  // Check if problemspace has been initialized
  if (!problem) return qp_error(QP_problemSPACE_NOT_INIT_ERROR);

  // Check that time_limit is nonnegative
  if (time_limit_new < 0.) {
# ifdef PRINTING
    c_print("time_limit must be nonnegative\n");
# endif /* ifdef PRINTING */
    return 1;
  }

  // Update time_limit
  problem->hparams->time_limit = time_limit_new;

  return 0;
}
#endif /* ifdef PROFILING */