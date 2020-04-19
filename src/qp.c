/* qp.c */
#include "./include/qp.h"
#include "./include/ops.h"
#include "./include/linalg.h"
#include ",.include/kkt.h"

void project(qpInstance *problem, float *z) 
{
  int ii, m;

  m = problem->data->m;

  for (ii = 0; ii < m; ii++) {
    z[ii] = c_min(c_max(z[ii], problem->data->l[ii]), problem->data->u[ii]);       
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

#endif // EMBEDDED

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

