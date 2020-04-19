/* csc.c */

#include "../include/csc.h"

static void* csc_malloc(int n, int size)
{
    return c_malloc(n * size);
}

static void* csc_calloc(int n, int size)
{
    return c_calloc(n, size);
}

csc* csc_matrix(int m, int n, int nzmax, float *x, int *i, int *p)
{
    csc *M = (csc *)c_malloc(sizeof(csc));

    if (!M) return 0;

    M->m = m;
    M->n = n;
    M->nz = -1;
    M->nzmax = nzmax;
    M->x = x;
    M->i = i;
    M->p = p;
    return M;
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
        if (A->p) c_free(A->p);
        if (A->i) c_free(A->i);
        if (A->x) c_free(A->x);
        c_free(A);
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





