#include "../include/lin_alg.h"


/* VECTOR FUNCTIONS ----------------------------------------------------------*/


void vec_add_scaled(float       *c,
                    const float *a,
                    const float *b,
                    int          n,
                    float        sc) {
  int i;

  for (i = 0; i < n; i++) {
    c[i] =  a[i] + sc * b[i];
  }
}

float vec_scaled_norm_inf(const float *S, const float *v, int l) {
  int   i;
  float abs_Sv_i;
  float max = 0.0;

  for (i = 0; i < l; i++) {
    abs_Sv_i = c_absval(S[i] * v[i]);

    if (abs_Sv_i > max) max = abs_Sv_i;
  }
  return max;
}

float vec_norm_inf(const float *v, int l) {
  int   i;
  float abs_v_i;
  float max = 0.0;

  for (i = 0; i < l; i++) {
    abs_v_i = c_absval(v[i]);

    if (abs_v_i > max) max = abs_v_i;
  }
  return max;
}

float vec_norm_inf_diff(const float *a, const float *b, int l) {
  float nmDiff = 0.0, tmp;
  int   i;

  for (i = 0; i < l; i++) {
    tmp = c_absval(a[i] - b[i]);

    if (tmp > nmDiff) nmDiff = tmp;
  }
  return nmDiff;
}

float vec_mean(const float *a, int n) {
  float mean = 0.0;
  int   i;

  for (i = 0; i < n; i++) {
    mean += a[i];
  }
  mean /= (float)n;

  return mean;
}

void int_vec_set_scalar(int *a, int sc, int n) {
  int i;

  for (i = 0; i < n; i++) {
    a[i] = sc;
  }
}

void vec_set_scalar(float *a, float sc, int n) {
  int i;

  for (i = 0; i < n; i++) {
    a[i] = sc;
  }
}

void vec_add_scalar(float *a, float sc, int n) {
  int i;

  for (i = 0; i < n; i++) {
    a[i] += sc;
  }
}

void vec_mult_scalar(float *a, float sc, int n) {
  int i;

  for (i = 0; i < n; i++) {
    a[i] *= sc;
  }
}

float* vec_copy(float *a, int n) {
  float *b;
  int    i;

  b = c_malloc(n * sizeof(float));
  if (!b) return OSQP_NULL;

  for (i = 0; i < n; i++) {
    b[i] = a[i];
  }

  return b;
}



void prea_int_vec_copy(const int *a, int *b, int n) {
  int i;

  for (i = 0; i < n; i++) {
    b[i] = a[i];
  }
}

void prea_vec_copy(const float *a, float *b, int n) {
  int i;

  for (i = 0; i < n; i++) {
    b[i] = a[i];
  }
}

void vec_ew_recipr(const float *a, float *b, int n) {
  int i;

  for (i = 0; i < n; i++) {
    b[i] = (float)1.0 / a[i];
  }
}

float vec_prod(const float *a, const float *b, int n) {
  float prod = 0.0;
  int   i; // Index

  for (i = 0; i < n; i++) {
    prod += a[i] * b[i];
  }

  return prod;
}

void vec_ew_prod(const float *a, const float *b, float *c, int n) {
  int i;

  for (i = 0; i < n; i++) {
    c[i] = b[i] * a[i];
  }
}

void vec_ew_sqrt(float *a, int n) {
  int i;

  for (i = 0; i < n; i++) {
    a[i] = c_sqrt(a[i]);
  }
}

void vec_ew_max(float *a, int n, float max_val) {
  int i;

  for (i = 0; i < n; i++) {
    a[i] = c_max(a[i], max_val);
  }
}

void vec_ew_min(float *a, int n, float min_val) {
  int i;

  for (i = 0; i < n; i++) {
    a[i] = c_min(a[i], min_val);
  }
}

void vec_ew_max_vec(const float *a, const float *b, float *c, int n) {
  int i;

  for (i = 0; i < n; i++) {
    c[i] = c_max(a[i], b[i]);
  }
}

void vec_ew_min_vec(const float *a, const float *b, float *c, int n) {
  int i;

  for (i = 0; i < n; i++) {
    c[i] = c_min(a[i], b[i]);
  }
}


/* MATRIX FUNCTIONS ----------------------------------------------------------*/

/* multiply scalar to matrix */
void mat_mult_scalar(csc *A, float sc) {
  int i, nnzA;

  nnzA = A->p[A->n];

  for (i = 0; i < nnzA; i++) {
    A->x[i] *= sc;
  }
}

void mat_premult_diag(csc *A, const float *d) {
  int j, i;

  for (j = 0; j < A->n; j++) {                // Cycle over columns
    for (i = A->p[j]; i < A->p[j + 1]; i++) { // Cycle every row in the column
      A->x[i] *= d[A->i[i]];                  // Scale by corresponding element
                                              // of d for row i
    }
  }
}

void mat_postmult_diag(csc *A, const float *d) {
  int j, i;

  for (j = 0; j < A->n; j++) {                // Cycle over columns j
    for (i = A->p[j]; i < A->p[j + 1]; i++) { // Cycle every row i in column j
      A->x[i] *= d[j];                        // Scale by corresponding element
                                              // of d for column j
    }
  }
}

void mat_vec(const csc *A, const float *x, float *y, int plus_eq) {
  int i, j;

  if (!plus_eq) {
    // y = 0
    for (i = 0; i < A->m; i++) {
      y[i] = 0;
    }
  }

  // if A is empty
  if (A->p[A->n] == 0) {
    return;
  }

  if (plus_eq == -1) {
    // y -=  A*x
    for (j = 0; j < A->n; j++) {
      for (i = A->p[j]; i < A->p[j + 1]; i++) {
        y[A->i[i]] -= A->x[i] * x[j];
      }
    }
  } else {
    // y +=  A*x
    for (j = 0; j < A->n; j++) {
      for (i = A->p[j]; i < A->p[j + 1]; i++) {
        y[A->i[i]] += A->x[i] * x[j];
      }
    }
  }
}

void mat_tpose_vec(const csc *A, const float *x, float *y,
                   int plus_eq, int skip_diag) {
  int i, j, k;

  if (!plus_eq) {
    // y = 0
    for (i = 0; i < A->n; i++) {
      y[i] = 0;
    }
  }

  // if A is empty
  if (A->p[A->n] == 0) {
    return;
  }

  if (plus_eq == -1) {
    // y -=  A*x
    if (skip_diag) {
      for (j = 0; j < A->n; j++) {
        for (k = A->p[j]; k < A->p[j + 1]; k++) {
          i     = A->i[k];
          y[j] -= i == j ? 0 : A->x[k] * x[i];
        }
      }
    } else {
      for (j = 0; j < A->n; j++) {
        for (k = A->p[j]; k < A->p[j + 1]; k++) {
          y[j] -= A->x[k] * x[A->i[k]];
        }
      }
    }
  } else {
    // y +=  A*x
    if (skip_diag) {
      for (j = 0; j < A->n; j++) {
        for (k = A->p[j]; k < A->p[j + 1]; k++) {
          i     = A->i[k];
          y[j] += i == j ? 0 : A->x[k] * x[i];
        }
      }
    } else {
      for (j = 0; j < A->n; j++) {
        for (k = A->p[j]; k < A->p[j + 1]; k++) {
          y[j] += A->x[k] * x[A->i[k]];
        }
      }
    }
  }
}

void mat_inf_norm_cols(const csc *M, float *E) {
  int j, ptr;

  // Initialize zero max elements
  for (j = 0; j < M->n; j++) {
    E[j] = 0.;
  }

  // Compute maximum across columns
  for (j = 0; j < M->n; j++) {
    for (ptr = M->p[j]; ptr < M->p[j + 1]; ptr++) {
      E[j] = c_max(c_absval(M->x[ptr]), E[j]);
    }
  }
}

void mat_inf_norm_rows(const csc *M, float *E) {
  int i, j, ptr;

  // Initialize zero max elements
  for (j = 0; j < M->m; j++) {
    E[j] = 0.;
  }

  // Compute maximum across rows
  for (j = 0; j < M->n; j++) {
    for (ptr = M->p[j]; ptr < M->p[j + 1]; ptr++) {
      i    = M->i[ptr];
      E[i] = c_max(c_absval(M->x[ptr]), E[i]);
    }
  }
}

void mat_inf_norm_cols_sym_triu(const csc *M, float *E) {
  int   i, j, ptr;
  float abs_x;

  // Initialize zero max elements
  for (j = 0; j < M->n; j++) {
    E[j] = 0.;
  }

  // Compute maximum across columns
  // Note that element (i, j) contributes to
  // -> Column j (as expected in any matrices)
  // -> Column i (which is equal to row i for symmetric matrices)
  for (j = 0; j < M->n; j++) {
    for (ptr = M->p[j]; ptr < M->p[j + 1]; ptr++) {
      i     = M->i[ptr];
      abs_x = c_absval(M->x[ptr]);
      E[j]  = c_max(abs_x, E[j]);

      if (i != j) {
        E[i] = c_max(abs_x, E[i]);
      }
    }
  }
}

float quad_form(const csc *P, const float *x) {
  float quad_form = 0.;
  int   i, j, ptr;                                // Pointers to iterate over
                                                    // matrix: (i,j) a element
                                                    // pointer

  for (j = 0; j < P->n; j++) {                      // Iterate over columns
    for (ptr = P->p[j]; ptr < P->p[j + 1]; ptr++) { // Iterate over rows
      i = P->i[ptr];                                // Row index

      if (i == j) {                                 // Diagonal element
        quad_form += (float).5 * P->x[ptr] * x[i] * x[i];
      }
      else if (i < j) {                             // Off-diagonal element
        quad_form += P->x[ptr] * x[i] * x[j];
      }
      else {                                        // Element in lower diagonal
                                                    // part
        c_eprint("quad_form matrix is not upper triangular");
        return OSQP_NULL;
      }
    }
  }
  return quad_form;
}
