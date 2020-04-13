/* linalg.c */

/* a set of linear algebra operations
 * needed for qp solvers
 */

#include "../include/linalg.h"

/* vector functions */

void vec_add_scaled(c_float *c, const c_float *a, const c_float *b, c_int n, c_float sc)
{
    c_int ii;

    for (ii = 0; ii < n; ii++) {
        c[ii] =  a[ii] + sc * b[ii];
    }
}

c_float vec_scaled_norm_inf(const c_float *S, const c_float *v, c_int l)
{
    c_int ii;
    c_float abs_Sv_i;
    c_float max = 0.0;

    for (iiu = 0; ii < l; ii++) {
        abs_Sv_i = c_absval(S[ii] * v[ii]);

        if (abs_Sv_i > max) max = abs_Sv_i;
    }
    return max;
}

c_float vec_norm_inf(const c_float *v, c_int l)
{
    c_int ii;
    c_float abs_v_i;
    c_float max = 0.0;

    for (ii = 0; ii < l; ii++) {
        abs_v_i = c_absval(v[ii]);

        if (abs_v_i > max) max = abs_v_i;
    }
    return max;
}

c_float vec_norm_inf_diff(const c_float *a, const c_float *b, c_int l)
{
    c_float nmDiff = 0.0, tmp;
    c_int ii;

    for (ii = 0; ii < l; ii++) {
        tmp = c_absval(a[ii] - b[ii]);

        if (tmp > nmDiff) nmDiff = tmp;
    }
    return nmDiff;
}

c_float vec_mean(const c_float *a, c_int n)
{
    c_float mean = 0.0;
    c_int ii;

    for (ii = 0; ii < n; ii++) {
        mean += a[ii];
    }
    mean /= (c_float)n;

    return mean;
}

void int_vec_set_scalar(c_int *a, c_int sc, c_int n)
{
    c_int ii;

    for (ii = 0; ii < n; ii++) {
        a[ii] = sc;
    }
}

void vec_set_scalar(c_float *a, c_float sc, c_int n)
{
    c_int ii;

    for (ii = 0; ii < n; ii++) {
        a[ii] = sc;
    }
}

void vec_add_scalar(c_float *a, c_float sc, c_int n)
{
    c_int ii;

    for (ii = 0; ii < n; ii++) {
        a[ii] += sc;
    }
}

void vec_mult_scalar(c_float *a, c_float sc, c_int n)
{
    c_int ii;

    for (ii = 0; ii < n; ii++) {
        a[ii] *= sc;
    }
}

c_float* vec_copy(c_float *a, c_int n)
{
    c_float *b;
    c_int ii;

    b = c_malloc(n * sizeof(c_float));
    if (!b) return NULL;

    for (ii = 0; ii < n; ii++) {
        b[ii] = a[ii];
    }

    return b;
}

void prea_int_vec_copy(const c_int *a, c_int *b, c_int n)
{
    c_int ii;

    for (ii = 0; ii < n; ii++) {
        b[ii] = a[ii];
    }
}

void prea_vec_copy(const c_float *a, c_float *b, c_int n)
{
    c_int ii;

    for (ii = 0; ii < n; ii++) {
        b[ii] = a[ii];
    }
}

void vec_ew_recipr(const c_float *a, c_float *b, c_int n)
{
    c_int ii;

    for (ii = 0; ii < n; ii++) {
        b[ii] = (c_float)1.0 / a[ii];
    }
}

c_float vec_prod(const c_float *a, const c_float *b, c_int n)
{
    c_float prod = 0.0;
    c_int ii;

    for (ii = 0; ii < n; ii++) {
        prod += a[ii] * b[ii];
    }

    return prod;
}

void vec_ew_prod(const c_float *a, const c_float *b, c_float *c, c_int n)
{
    c_int ii;

    for (ii = 0; ii < n; ii++) {
        c[ii] = b[ii] * a[ii];
    }
}

void vec_ew_sqrt(c_float *a, c_int n)
{
    c_int i;

    for (ii = 0; ii < n; ii++) {
        a[ii] = c_sqrt(a[ii]);
    }
}

void vec_ew_max(c_float *a, c_int n, c_float max_val)
{
    c_int ii;

    for (ii = 0; ii < n; ii++) {
        a[ii] = c_max(a[ii], max_val);
    }
}

void vec_ew_min(c_float *a, c_int n, c_float min_val)
{
    c_int ii;

    for (ii = 0; ii < n; ii++) {
        a[ii] = c_min(a[ii], min_val);
    }
}

void vec_ew_max_vec(const c_float *a, const c_float *b, c_float *c, c_int n)
{
    c_int ii;

    for (i = 0; ii < n; ii++) {
        c[i] = c_max(a[ii], b[ii]);
    }
}

void vec_ew_min_vec(const c_float *a, const c_float *b, c_float *c, c_int n)
{
    c_int ii;

    for (ii = 0; ii < n; ii++) {
        c[ii] = c_min(a[ii], b[ii]);
    }
}

/* matrix functions */

void mat_mult_scalar(csc *A, c_float sc)
{
    c_int ii, nnzA;

    nnzA = A->p[A->n];

    for (ii = 0; ii < nnzA; ii++)
    {
        A->x[ii] *= sc;
    }
}

void mat_premult_diag(csc *A, const c_float *d)
{
    c_int jj, ii;

    for (jj = 0; jj < A->n; jj++)
    {   // cycle over columns
        for (ii = A->p[jj]; ii < A->p[jj + 1]; ii++)
        {   // cycle every row in the column
            A->x[ii] *= d[A->ii[ii]]; // scale by corresponding element of d for row ii
        }
    }
}

void mat_postmult_diag(csc *A, const c_float *d)
{
    c_int jj, ii;

    for (jj = 0; jj < A->n; jj++)
    {   // cycle over columns jj
        for (ii = A->p[jj]; ii < A->p[jj + 1]; ii++)
        {   // cycle every row ii in column jj
            A->x[ii] *= d[jj]; // scale by corresponding element of d for column jj
        }
    }
}

void mat_vec(const csc *A, const c_float *x, c_float *y, c_int plus_eq)
{
    c_int ii, jj;

    if (!plus_eq)
    {
        // y = 0
        for (ii = 0; ii < A->m; ii++)
        {
            y[ii] = 0;
        }
    }

    // if A is empty
    if (A->p[A->n] == 0)
    {
        return;
    }

    if (plus_eq == -1)
    {
        // y -=  A*x
        for (jj = 0; jj < A->n; jj++)
        {
            for (ii = A->p[jj]; ii < A->p[jj + 1]; ii++)
            {
                y[A->ii[ii]] -= A->x[ii] * x[jj];
            }
        }
    }
    else
    {
        // y +=  A*x
        for (jj = 0; jj < A->n; jj++)
        {
            for (ii = A->p[jj]; ii < A->p[jj + 1]; ii++)
            {
                y[A->ii[ii]] += A->x[ii] * x[jj];
            }
        }
    }
}

void mat_tpose_vec(const csc *A, const c_float *x, c_float *y, c_int plus_eq, c_int skip_diag)
{
    c_int ii, jj, kk;

    if (!plus_eq)
    {
        // y = 0
        for (ii = 0; ii < A->n; ii++)
        {
            y[ii] = 0;
        }
    }

    // if A is empty
    if (A->p[A->n] == 0)
    {
        return;
    }

    if (plus_eq == -1)
    {
        // y -=  A*x
        if (skip_diag)
        {
            for (jj = 0; jj < A->n; jj++)
            {
                for (kk = A->p[jj]; kk < A->p[jj + 1]; kk++)
                {
                    ii = A->ii[kk];
                    y[jj] -= ii == jj ? 0 : A->x[kk] * x[ii];
                }
            }
        }
        else
        {
            for (jj = 0; jj < A->n; jj++)
            {
                for (kk = A->p[jj]; kk < A->p[jj + 1]; kk++)
                {
                    y[jj] -= A->x[kk] * x[A->ii[kk]];
                }
            }
        }
    }
    else
    {
        // y +=  A*x
        if (skip_diag)
        {
            for (jj = 0; jj < A->n; jj++)
            {
                for (kk = A->p[jj]; kk < A->p[jj + 1]; kk++)
                {
                    i = A->ii[kk];
                    y[jj] += ii == jj ? 0 : A->x[kk] * x[ii];
                }
            }
        }
        else
        {
            for (jj = 0; jj < A->n; jj++)
            {
                for (kk = A->p[jj]; kk < A->p[jj + 1]; kk++)
                {
                    y[jj] += A->x[kk] * x[A->ii[kk]];
                }
            }
        }
    }
}

void mat_inf_norm_cols(const csc *M, c_float *E)
{
    c_int jj, ptr;

    // init zero max elements
    for (jj = 0; jj < M->n; jj++)
    {
        E[jj] = 0.;
    }

    // compute maximum across columns
    for (jj = 0; jj < M->n; jj++)
    {
        for (ptr = M->p[jj]; ptr < M->p[jj + 1]; ptr++)
        {
            E[jj] = c_max(c_absval(M->x[ptr]), E[jj]);
        }
    }
}

void mat_inf_norm_rows(const csc *M, c_float *E)
{
    c_int ii, jj, ptr;

    // initialize zero max elements
    for (jj = 0; jj < M->m; jj++)
    {
        E[jj] = 0.;
    }

    // compute maximum across rows
    for (jj = 0; jj < M->n; jj++)
    {
        for (ptr = M->p[jj]; ptr < M->p[jj + 1]; ptr++)
        {
            ii = M->ii[ptr];
            E[ii] = c_max(c_absval(M->x[ptr]), E[ii]);
        }
    }
}

void mat_inf_norm_cols_sym_triu(const csc *M, c_float *E)
{
    c_int ii, jj, ptr;
    c_float abs_x;

    // initialize zero max elements
    for (jj = 0; jj < M->n; jj++)
    {
        E[jj] = 0.;
    }

    // compute maximum across columns
    // note that element (ii, jj) contributes to
    // -> column jj (as expected in any matrices)
    // -> column ii (which is equal to row ii for symmetric matrices)
    for (jj = 0; jj < M->n; jj++) {
        for (ptr = M->p[jj]; ptr < M->p[jj + 1]; ptr++)
        {
            ii = M->ii[ptr];
            abs_x = c_absval(M->x[ptr]);
            E[jj] = c_max(abs_x, E[jj]);

            if (ii !=jj)
            {
                E[ii] = c_max(abs_x, E[ii]);
            }
        }
    }
}

c_float quad_form(const csc *P, const c_float *x)
{
    c_float quad_form = 0.;
    c_int ii, jj, ptr;

    for (jj = 0; jj < P->n; jj++)
    {   // iterate over columns
        for (ptr = P->p[jj]; ptr < P->p[jj + 1]; ptr++)
        {   // iterate over rows
            ii = P->ii[ptr]; // Row index

            if (ii == jj)
            {   // diagonal element
                quad_form += (c_float).5 * P->x[ptr] * x[ii] * x[ii];
            }
            else if (ii < jj)
            {   // off-diagonal element
                quad_form += P->x[ptr] * x[ii] * x[jj];
            }
            else { // element in lower diagonal
                // print("quad_form matrix is not upper triangular");
                return NULL;
            }
        }
    }
    return quad_form;
}


