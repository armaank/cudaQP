/* main.c 
* testing qp solver
*/
#include "./include/qp.h"


int main(int argc, char **argv) {
  // Load problem data
  float P_x[3] = { 4.0, 1.0, 2.0, };
  int   P_nnz  = 3;
  int   P_i[3] = { 0, 0, 1, };
  int   P_p[3] = { 0, 1, 3, };
  float q[2]   = { 1.0, 1.0, };
  float A_x[4] = { 1.0, 1.0, 1.0, 1.0, };
  int   A_nnz  = 4;
  int   A_i[4] = { 0, 1, 0, 2, };
  int   A_p[3] = { 0, 2, 4, };
  float l[3]   = { 1.0, 0.0, 0.0, };
  float u[3]   = { 1.0, 0.7, 0.7, };
  int n = 2;
  int m = 3;

  // Exitflag
  int exitflag = 0;

  // Workspace structures
  qpInstance *problem;
  qpHyperparams  *params = (qpHyperparams *)c_malloc(sizeof(qpHyperparams));
  qpData      *data     = (qpData *)c_malloc(sizeof(qpData));

  // Populate data
  if (data) {
    data->n = n;
    data->m = m;
    data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
    data->q = q;
    data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
    data->l = l;
    data->u = u;
  }

  // Define solver params as default
  if (params) qp_set_default_params(params);

  // Setup problemspace
  exitflag = qp_setup(&problem, data, params);

  // Solve Problem
  qp_solve(problem);

  // Clean problemspace
  qp_cleanup(problem);
  if (data) {
    if (data->A) c_free(data->A);
    if (data->P) c_free(data->P);
    c_free(data);
  }
  if (params)  c_free(params);

  return exitflag;
}


