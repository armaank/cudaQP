#ifndef CSC_H
#define CSC_H

#include <stdlib.h>
#include <stdio.h>
#include "linalg.h"
#include "ops.h"
/* create Compressed-Column-Sparse matrix from existing arrays */
csc* csc_matrix(int m, int n, int nzmax, float *x, int *i, int *p);

/* create uninitialized CSC matrix  */
csc* csc_spalloc(int m, int n, int nzmax, int values, int triplet);

/* free sparse matrix */
void csc_spfree(csc *A);

/*  Copy sparse CSC matrix A to output */
csc* copy_csc_mat(const csc *A);

/* copy sparse CSC matrix A to B (B is preallocated, NO MALOC) */
void prea_copy_csc_mat(const csc *A, csc *B);

/* C = compressed-column CSC from matrix T in triplet form */
csc* triplet_to_csc(const csc *T, int *TtoC);

/* C = compressed-row CSR from matrix T in triplet form */
csc* triplet_to_csr(const csc *T, int *TtoC);

/* convert sparse to dense */
float* csc_to_dns(csc *M);

/* convert square CSC matrix into upper triangular one */
csc* csc_to_triu(csc *M);

/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
int csc_cumsum(int *p, int *c, int n);

/* compute inverse of permutation matrix stored in the vector p */
int* csc_pinv(int const *p, int n);

/* C = A(p,p)= PAP' where A and C are symmetric the upper part stored */
csc* csc_symperm(const csc *A, const int *pinv, int *AtoC, int values);

csc* csc_done(csc  *C,
              void *w,
              void *x,
              int ok);

#endif //CSC_H
