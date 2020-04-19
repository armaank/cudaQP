#include "../include/cs.h"


static void* csc_malloc(int n, int size) {
  return c_malloc(n * size);
}

static void* csc_calloc(int n, int size) {
  return c_calloc(n, size);
}

csc* csc_matrix(int m, int n, int nzmax, float *x, int *i, int *p)
{
  csc *M = (csc *)c_malloc(sizeof(csc));

  if (!M) return OSQP_NULL;

  M->m     = m;
  M->n     = n;
  M->nz    = -1;
  M->nzmax = nzmax;
  M->x     = x;
  M->i     = i;
  M->p     = p;
  return M;
}

csc* csc_spalloc(int m, int n, int nzmax, int values, int triplet) {
  csc *A = csc_calloc(1, sizeof(csc)); /* allocate the csc struct */

  if (!A) return OSQP_NULL;            /* out of memory */

  A->m     = m;                        /* define dimensions and nzmax */
  A->n     = n;
  A->nzmax = nzmax = c_max(nzmax, 1);
  A->nz    = triplet ? 0 : -1;         /* allocate triplet or comp.col */
  A->p     = csc_malloc(triplet ? nzmax : n + 1, sizeof(int));
  A->i     = csc_malloc(nzmax,  sizeof(int));
  A->x     = values ? csc_malloc(nzmax,  sizeof(float)) : OSQP_NULL;
  if (!A->p || !A->i || (values && !A->x)){
    csc_spfree(A);
    return OSQP_NULL;
  } else return A;
}

void csc_spfree(csc *A) {
  if (A){
    if (A->p) c_free(A->p);
    if (A->i) c_free(A->i);
    if (A->x) c_free(A->x);
    c_free(A);
  }
}

csc* triplet_to_csc(const csc *T, int *TtoC) {
  int m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj;
  float *Cx, *Tx;
  csc     *C;

  m  = T->m;
  n  = T->n;
  Ti = T->i;
  Tj = T->p;
  Tx = T->x;
  nz = T->nz;
  C  = csc_spalloc(m, n, nz, Tx != OSQP_NULL, 0);     /* allocate result */
  w  = csc_calloc(n, sizeof(int));                  /* get workspace */

  if (!C || !w) return csc_done(C, w, OSQP_NULL, 0);  /* out of memory */

  Cp = C->p;
  Ci = C->i;
  Cx = C->x;

  for (k = 0; k < nz; k++) w[Tj[k]]++;  /* column counts */
  csc_cumsum(Cp, w, n);                 /* column pointers */

  for (k = 0; k < nz; k++) {
    Ci[p = w[Tj[k]]++] = Ti[k];         /* A(i,j) is the pth entry in C */

    if (Cx) {
      Cx[p] = Tx[k];

      if (TtoC != OSQP_NULL) TtoC[k] = p;  // Assign vector of indices
    }
  }
  return csc_done(C, w, OSQP_NULL, 1);     /* success; free w and return C */
}

csc* triplet_to_csr(const csc *T, int *TtoC) {
  int m, n, nz, p, k, *Cp, *Cj, *w, *Ti, *Tj;
  float *Cx, *Tx;
  csc     *C;

  m  = T->m;
  n  = T->n;
  Ti = T->i;
  Tj = T->p;
  Tx = T->x;
  nz = T->nz;
  C  = csc_spalloc(m, n, nz, Tx != OSQP_NULL, 0);     /* allocate result */
  w  = csc_calloc(m, sizeof(int));                  /* get workspace */

  if (!C || !w) return csc_done(C, w, OSQP_NULL, 0);  /* out of memory */

  Cp = C->p;
  Cj = C->i;
  Cx = C->x;

  for (k = 0; k < nz; k++) w[Ti[k]]++;  /* row counts */
  csc_cumsum(Cp, w, m);                 /* row pointers */

  for (k = 0; k < nz; k++) {
    Cj[p = w[Ti[k]]++] = Tj[k];         /* A(i,j) is the pth entry in C */

    if (Cx) {
      Cx[p] = Tx[k];

      if (TtoC != OSQP_NULL) TtoC[k] = p;  // Assign vector of indices
    }
  }
  return csc_done(C, w, OSQP_NULL, 1);     /* success; free w and return C */
}

int csc_cumsum(int *p, int *c, int n) {
  int i, nz = 0;

  if (!p || !c) return -1;  /* check inputs */

  for (i = 0; i < n; i++)
  {
    p[i] = nz;
    nz  += c[i];
    c[i] = p[i];
  }
  p[n] = nz;
  return nz; /* return sum (c [0..n-1]) */
}

int* csc_pinv(int const *p, int n) {
  int k, *pinv;

  if (!p) return OSQP_NULL;                /* p = OSQP_NULL denotes identity */

  pinv = csc_malloc(n, sizeof(int));     /* allocate result */

  if (!pinv) return OSQP_NULL;             /* out of memory */

  for (k = 0; k < n; k++) pinv[p[k]] = k;  /* invert the permutation */
  return pinv;                             /* return result */
}

csc* csc_symperm(const csc *A, const int *pinv, int *AtoC, int values) {
  int i, j, p, q, i2, j2, n, *Ap, *Ai, *Cp, *Ci, *w;
  float *Cx, *Ax;
  csc     *C;

  n  = A->n;
  Ap = A->p;
  Ai = A->i;
  Ax = A->x;
  C  = csc_spalloc(n, n, Ap[n], values && (Ax != OSQP_NULL),
                   0);                                /* alloc result*/
  w = csc_calloc(n, sizeof(int));                   /* get workspace */

  if (!C || !w) return csc_done(C, w, OSQP_NULL, 0);  /* out of memory */

  Cp = C->p;
  Ci = C->i;
  Cx = C->x;

  for (j = 0; j < n; j++)    /* count entries in each column of C */
  {
    j2 = pinv ? pinv[j] : j; /* column j of A is column j2 of C */

    for (p = Ap[j]; p < Ap[j + 1]; p++) {
      i = Ai[p];

      if (i > j) continue;     /* skip lower triangular part of A */
      i2 = pinv ? pinv[i] : i; /* row i of A is row i2 of C */
      w[c_max(i2, j2)]++;      /* column count of C */
    }
  }
  csc_cumsum(Cp, w, n);        /* compute column pointers of C */

  for (j = 0; j < n; j++) {
    j2 = pinv ? pinv[j] : j;   /* column j of A is column j2 of C */

    for (p = Ap[j]; p < Ap[j + 1]; p++) {
      i = Ai[p];

      if (i > j) continue;                             /* skip lower triangular
                                                          part of A*/
      i2                         = pinv ? pinv[i] : i; /* row i of A is row i2
                                                          of C */
      Ci[q = w[c_max(i2, j2)]++] = c_min(i2, j2);

      if (Cx) Cx[q] = Ax[p];

      if (AtoC) { // If vector AtoC passed, store values of the mappings
        AtoC[p] = q;
      }
    }
  }
  return csc_done(C, w, OSQP_NULL, 1); /* success; free workspace, return C */
}

csc* copy_csc_mat(const csc *A) {
  csc *B = csc_spalloc(A->m, A->n, A->p[A->n], 1, 0);

  if (!B) return OSQP_NULL;

  prea_int_vec_copy(A->p, B->p, A->n + 1);
  prea_int_vec_copy(A->i, B->i, A->p[A->n]);
  prea_vec_copy(A->x, B->x, A->p[A->n]);

  return B;
}

void prea_copy_csc_mat(const csc *A, csc *B) {
  prea_int_vec_copy(A->p, B->p, A->n + 1);
  prea_int_vec_copy(A->i, B->i, A->p[A->n]);
  prea_vec_copy(A->x, B->x, A->p[A->n]);

  B->nzmax = A->nzmax;
}

csc* csc_done(csc *C, void *w, void *x, int ok) {
  c_free(w);                   /* free workspace */
  c_free(x);
  if (ok) return C;
  else {
    csc_spfree(C);
    return OSQP_NULL;
  }
}

csc* csc_to_triu(csc *M) {
  csc  *M_trip;    // Matrix in triplet format
  csc  *M_triu;    // Resulting upper triangular matrix
  int nnzorigM;  // Number of nonzeros from original matrix M
  int nnzmaxM;   // Estimated maximum number of elements of upper triangular M
  int n;         // Dimension of M
  int ptr, i, j; // Counters for (i,j) and index in M
  int z_M = 0;   // Counter for elements in M_trip


  // Check if matrix is square
  if (M->m != M->n) {
    c_eprint("Matrix M not square");
    return OSQP_NULL;
  }
  n = M->n;

  // Get number of nonzeros full M
  nnzorigM = M->p[n];

  // Estimate nnzmaxM
  // Number of nonzero elements in original M + diagonal part.
  // -> Full matrix M as input: estimate is half the number of total elements +
  // diagonal = .5 * (nnzorigM + n)
  // -> Upper triangular matrix M as input: estimate is the number of total
  // elements + diagonal = nnzorigM + n
  // The maximum between the two is nnzorigM + n
  nnzmaxM = nnzorigM + n;

  // OLD
  // nnzmaxM = n*(n+1)/2;  // Full upper triangular matrix (This version
  // allocates too much memory!)
  // nnzmaxM = .5 * (nnzorigM + n);  // half of the total elements + diagonal

  // Allocate M_trip
  M_trip = csc_spalloc(n, n, nnzmaxM, 1, 1); // Triplet format

  if (!M_trip) {
    c_eprint("Upper triangular matrix extraction failed (out of memory)");
    return OSQP_NULL;
  }

  // Fill M_trip with only elements in M which are in the upper triangular
  for (j = 0; j < n; j++) { // Cycle over columns
    for (ptr = M->p[j]; ptr < M->p[j + 1]; ptr++) {
      // Get row index
      i = M->i[ptr];

      // Assign element only if in the upper triangular
      if (i <= j) {
        // c_print("\nM(%i, %i) = %.4f", M->i[ptr], j, M->x[ptr]);

        M_trip->i[z_M] = i;
        M_trip->p[z_M] = j;
        M_trip->x[z_M] = M->x[ptr];

        // Increase counter for the number of elements
        z_M++;
      }
    }
  }

  // Set number of nonzeros
  M_trip->nz = z_M;

  // Convert triplet matrix to csc format
  M_triu = triplet_to_csc(M_trip, OSQP_NULL);

  // Assign number of nonzeros of full matrix to triu M
  M_triu->nzmax = nnzmaxM;

  // Cleanup and return result
  csc_spfree(M_trip);

  // Return matrix in triplet form
  return M_triu;
}
