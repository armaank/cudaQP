#ifndef CSC_H
# define CSC_H

# include "types.h"   // CSC matrix type
# include "linalg.h" // Vector copy operations

/* create Compressed-Column-Sparse matrix from existing arrays */
csc* csc_matrix(c_int m, c_int n, c_int nzmax, c_float *x, c_int *i, c_int *p);


/* create uninitialized CSC matrix  */
csc* csc_spalloc(c_int m, c_int n, c_int nzmax, c_int values, c_int triplet);

/* free sparse matrix */
void csc_spfree(csc *A);

/*  Copy sparse CSC matrix A to output */
csc* copy_csc_mat(const csc *A);

/* copy sparse CSC matrix A to B (B is preallocated, NO MALOC) */
void prea_copy_csc_mat(const csc *A, csc *B);

/* C = compressed-column CSC from matrix T in triplet form */
csc* triplet_to_csc(const csc *T, c_int *TtoC);

/* C = compressed-row CSR from matrix T in triplet form */
csc* triplet_to_csr(const csc *T, c_int *TtoC);

/* convert sparse to dense */
c_float* csc_to_dns(csc *M);

/* convert square CSC matrix into upper triangular one */
csc* csc_to_triu(csc *M);

/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
c_int csc_cumsum(c_int *p, c_int *c, c_int n);

/* compute inverse of permutation matrix stored in the vector p */
c_int* csc_pinv(c_int const *p, c_int n);

/* C = A(p,p)= PAP' where A and C are symmetric the upper part stored */
csc* csc_symperm(const csc *A, const c_int *pinv, c_int *AtoC, c_int values);

#endif // ifndef CSC_H
