/* csc.c */
#include "../include/csc.h"
#include <stdlib.h>
#include <stdio.h>
#include "../include/ops.h"
// #include "../include/linalg.h"

static void* csc_malloc(int n, int size)
{
    return malloc(n * size);
}

static void* csc_calloc(int n, int size)
{
    return calloc(n, size);
}

csc* csc_spalloc(int m, int n, int nzmax, int values, int triplet)
{
    /* allocate the csc struct */
    csc *A = csc_calloc(1, sizeof(csc));

    /* out of memory */
    if (!A) return 0;

    /* define dimensions and nzmax */
    A->m = m;
    A->n = n;
    A->nzmax = nzmax = c_max(nzmax, 1);
    A->nz = triplet ? 0 : -1; /* allocate triplet or comp.col */
    A->p = csc_malloc(triplet ? nzmax : n + 1, sizeof(int));
    A->i = csc_malloc(nzmax, sizeof(int));
    A->x = values ? csc_malloc(nzmax, sizeof(float)) : 0;

    if (!A->p || !A->i || (values && !A->x))
    {
        csc_spfree(A);
        return 0;
    }
    else
        return A;
}

void csc_spfree(csc *A)
{
    if (A)
    {
        if (A->p) free(A->p);
        if (A->i) free(A->i);
        if (A->x) free(A->x);
        free(A);
    }
}

int csc_cumsum(int *p, int *c, int n)
{
    int ii, nz = 0;

    /* check inputs */
    if (!p || !c)
        return -1;

    for (ii = 0; ii < n; ii++)
    {
        p[ii] = nz;
        nz  += c[ii];
        c[ii] = p[ii];
    }
    p[n] = nz;

    /* return sum (c [0..n-1]) */
    return nz;
}

int* csc_pinv(int const *p, int n)
{
    int kk, *pinv;

    /*0 denotes identity - might need to return null here */
    if (!p) return 0;

    /* allocate result */
    pinv = csc_malloc(n, sizeof(int));

    /* out of memory */
    if (!pinv) return 0;

    for (kk = 0; kk < n; kk++)
        pinv[p[kk]] = kk; // invert the permutation

    /* return result */
    return pinv;
}

csc* copy_csc_mat(const csc *A)
{
    csc *B = csc_spalloc(A->m, A->n, A->p[A->n], 1, 0);

    if (!B)
        return 0;

    prea_int_vec_copy(A->p, B->p, A->n + 1);
    prea_int_vec_copy(A->i, B->i, A->p[A->n]);
    prea_vec_copy(A->x, B->x, A->p[A->n]);

    return B;
}

void prea_copy_csc_mat(const csc *A, csc *B)
{
    prea_int_vec_copy(A->p, B->p, A->n + 1);
    prea_int_vec_copy(A->i, B->i, A->p[A->n]);
    prea_vec_copy(A->x, B->x, A->p[A->n]);

    B->nzmax = A->nzmax;
}
csc* csc_done(csc *C, void *w, void *x, int ok) {
    free(w);                   /* free workspace */
    free(x);
    if (ok) return C;
    else {
        csc_spfree(C);
        return 0;
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
    C  = csc_spalloc(m, n, nz, Tx != 0, 0);     /* allocate result */
    w  = csc_calloc(n, sizeof(int));                  /* get workspace */

    if (!C || !w) return csc_done(C, w, 0, 0);  /* out of memory */

    Cp = C->p;
    Ci = C->i;
    Cx = C->x;

    for (k = 0; k < nz; k++) w[Tj[k]]++;  /* column counts */
    csc_cumsum(Cp, w, n);                 /* column pointers */

    for (k = 0; k < nz; k++) {
        Ci[p = w[Tj[k]]++] = Ti[k];         /* A(i,j) is the pth entry in C */

        if (Cx) {
            Cx[p] = Tx[k];

            if (TtoC != 0) TtoC[k] = p;  // Assign vector of indices
        }
    }
    return csc_done(C, w, 0, 1);     /* success; free w and return C */
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
        printf("Matrix M not square");
        return 0;
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
        printf("Upper triangular matrix extraction failed (out of memory)");
        return 0;
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
    M_triu = triplet_to_csc(M_trip, 0);

    // Assign number of nonzeros of full matrix to triu M
    M_triu->nzmax = nnzmaxM;

    // Cleanup and return result
    csc_spfree(M_trip);

    // Return matrix in triplet form
    return M_triu;
}


