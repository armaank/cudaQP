/* csc.c */

#include "../include/csc.h"



static void* csc_malloc(c_int n, c_int size)
{
    return c_malloc(n * size);
}

static void* csc_calloc(c_int n, c_int size)
{
    return c_calloc(n, size);
}

csc* csc_matrix(c_int m, c_int n, c_int nzmax, c_float *x, c_int *i, c_int *p)
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

csc* csc_spalloc(c_int m, c_int n, c_int nzmax, c_int values, c_int triplet)
{
    csc *A = csc_calloc(1, sizeof(csc)); /* allocate the csc struct */

    if (!A) return NULL; /* out of memory */

    A->m = m; /* define dimensions and nzmax */
    A->n = n;
    A->nzmax = nzmax = c_max(nzmax, 1);
    A->nz = triplet ? 0 : -1; /* allocate triplet or comp.col */
    A->p = csc_malloc(triplet ? nzmax : n + 1, sizeof(c_int));
    A->i = csc_malloc(nzmax, sizeof(c_int));
    A->x = values ? csc_malloc(nzmax, sizeof(c_float)) : NULL;

    if (!A->p || !A->i || (values && !A->x))
    {
        csc_spfree(A);
        return NULL;
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

c_int csc_cumsum(c_int *p, c_int *c, c_int n)
{
    c_int ii, nz = 0;

    if (!p || !c) return -1;  /* check inputs */

    for (ii = 0; ii < n; ii++)
    {
        p[ii] = nz;
        nz  += c[ii];
        c[ii] = p[ii];
    }
    p[n] = nz;
    return nz; /* return sum (c [0..n-1]) */
}

c_int* csc_pinv(c_int const *p, c_int n)
{
    c_int kk, *pinv;

    if (!p) return NULL; /*_NULL denotes identity */

    pinv = csc_malloc(n, sizeof(c_int)); /* allocate result */

    if (!pinv) return NULL; /* out of memory */

    for (kk = 0; kk < n; kk++)
        pinv[p[kk]] = kk; /* invert the permutation */
    return pinv; /* return result */
}

csc* copy_csc_mat(const csc *A)
{
    csc *B = csc_spalloc(A->m, A->n, A->p[A->n], 1, 0);

    if (!B) return NULL;

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





