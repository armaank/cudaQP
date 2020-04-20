/* linalg.c */

/* a set of linear algebra operations
 * needed for qp solvers
 */

#include "../include/linalg.h"
/* vector functions */

void vec_add_scaled(float *c, const float *a, const float *b, int n, float sc)
{
    int ii;

    for (ii = 0; ii < n; ii++)
    {
        c[ii] =  a[ii] + sc * b[ii];
    }
}

float vec_scaled_norm_inf(const float *S, const float *v, int l)
{
    int ii;
    float abs_Sv_i;
    float max = 0.0;

    for (ii = 0; ii < l; ii++)
    {
        abs_Sv_i = c_absval(S[ii] * v[ii]);

        if (abs_Sv_i > max) max = abs_Sv_i;
    }
    return max;
}

float vec_norm_inf(const float *v, int l)
{
    int ii;
    float abs_v_i;
    float max = 0.0;

    for (ii = 0; ii < l; ii++)
    {
        abs_v_i = c_absval(v[ii]);

        if (abs_v_i > max) max = abs_v_i;
    }
    return max;
}

float vec_norm_inf_diff(const float *a, const float *b, int l)
{
    float nmDiff = 0.0, tmp;
    int ii;

    for (ii = 0; ii < l; ii++)
    {
        tmp = c_absval(a[ii] - b[ii]);

        if (tmp > nmDiff) nmDiff = tmp;
    }
    return nmDiff;
}

float vec_mean(const float *a, int n)
{
    float mean = 0.0;
    int ii;

    for (ii = 0; ii < n; ii++)
    {
        mean += a[ii];
    }
    mean /= (float)n;

    return mean;
}

void int_vec_set_scalar(int *a, int sc, int n)
{
    int ii;

    for (ii = 0; ii < n; ii++)
    {
        a[ii] = sc;
    }
}

void vec_set_scalar(float *a, float sc, int n)
{
    int ii;

    for (ii = 0; ii < n; ii++)
    {
        a[ii] = sc;
    }
}

void vec_add_scalar(float *a, float sc, int n)
{
    int ii;

    for (ii = 0; ii < n; ii++)
    {
        a[ii] += sc;
    }
}

void vec_mult_scalar(float *a, float sc, int n)
{
    int ii;

    for (ii = 0; ii < n; ii++)
    {
        a[ii] *= sc;
    }
}

float* vec_copy(float *a, int n)
{
    float *b;
    int ii;

    b = malloc(n * sizeof(float));
    if (!b) return 0;

    for (ii = 0; ii < n; ii++)
    {
        b[ii] = a[ii];
    }

    return b;
}

void prea_int_vec_copy(const int *a, int *b, int n)
{
    int ii;

    for (ii = 0; ii < n; ii++)
    {
        b[ii] = a[ii];
    }
}

void prea_vec_copy(const float *a, float *b, int n)
{
    int ii;

    for (ii = 0; ii < n; ii++)
    {
        b[ii] = a[ii];
    }
}

void vec_ew_recipr(const float *a, float *b, int n)
{
    int ii;

    for (ii = 0; ii < n; ii++)
    {
        b[ii] = (float)1.0 / a[ii];
    }
}

float vec_prod(const float *a, const float *b, int n)
{
    float prod = 0.0;
    int ii;

    for (ii = 0; ii < n; ii++)
    {
        prod += a[ii] * b[ii];
    }

    return prod;
}

void vec_ew_prod(const float *a, const float *b, float *c, int n)
{
    int ii;

    for (ii = 0; ii < n; ii++)
    {
        c[ii] = b[ii] * a[ii];
    }
}

void vec_ew_sqrt(float *a, int n)
{
    int ii;

    for (ii = 0; ii < n; ii++)
    {
        a[ii] = c_sqrt(a[ii]);
    }
}

void vec_ew_max(float *a, int n, float max_val)
{
    int ii;

    for (ii = 0; ii < n; ii++)
    {
        a[ii] = c_max(a[ii], max_val);
    }
}

void vec_ew_min(float *a, int n, float min_val)
{
    int ii;

    for (ii = 0; ii < n; ii++)
    {
        a[ii] = c_min(a[ii], min_val);
    }
}

void vec_ew_max_vec(const float *a, const float *b, float *c, int n)
{
    int ii;

    for (ii = 0; ii < n; ii++)
    {
        c[ii] = c_max(a[ii], b[ii]);
    }
}

void vec_ew_min_vec(const float *a, const float *b, float *c, int n)
{
    int ii;

    for (ii = 0; ii < n; ii++)
    {
        c[ii] = c_min(a[ii], b[ii]);
    }
}

/* matrix functions */

void mat_mult_scalar(csc *A, float sc)
{
    int ii, nnzA;

    nnzA = A->p[A->n];

    for (ii = 0; ii < nnzA; ii++)
    {
        A->x[ii] *= sc;
    }
}

void mat_premult_diag(csc *A, const float *d)
{
    int jj, ii;

    for (jj = 0; jj < A->n; jj++)
    {   // cycle over columns
        for (ii = A->p[jj]; ii < A->p[jj + 1]; ii++)
        {   // cycle every row in the column
            A->x[ii] *= d[A->i[ii]]; // scale by corresponding element of d for row ii
        }
    }
}

void mat_postmult_diag(csc *A, const float *d)
{
    int jj, ii;

    for (jj = 0; jj < A->n; jj++)
    {   // cycle over columns jj
        for (ii = A->p[jj]; ii < A->p[jj + 1]; ii++)
        {   // cycle every row ii in column jj
            A->x[ii] *= d[jj]; // scale by corresponding element of d for column jj
        }
    }
}

void mat_vec(const csc *A, const float *x, float *y, int plus_eq)
{
    int ii, jj;

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
                y[A->i[ii]] -= A->x[ii] * x[jj];
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
                y[A->i[ii]] += A->x[ii] * x[jj];
            }
        }
    }
}

void mat_tpose_vec(const csc *A, const float *x, float *y, int plus_eq, int skip_diag)
{
    int ii, jj, kk;

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
                    ii = A->i[kk];
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
                    y[jj] -= A->x[kk] * x[A->i[kk]];
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
                    ii = A->i[kk];
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
                    y[jj] += A->x[kk] * x[A->i[kk]];
                }
            }
        }
    }
}

void mat_inf_norm_cols(const csc *M, float *E)
{
    int jj, ptr;

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

void mat_inf_norm_rows(const csc *M, float *E)
{
    int ii, jj, ptr;

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
            ii = M->i[ptr];
            E[ii] = c_max(c_absval(M->x[ptr]), E[ii]);
        }
    }
}

void mat_inf_norm_cols_sym_triu(const csc *M, float *E)
{
    int ii, jj, ptr;
    float abs_x;

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
            ii = M->i[ptr];
            abs_x = c_absval(M->x[ptr]);
            E[jj] = c_max(abs_x, E[jj]);

            if (ii !=jj)
            {
                E[ii] = c_max(abs_x, E[ii]);
            }
        }
    }
}

float quad_form(const csc *P, const float *x)
{
    float quad_form = 0.;
    int ii, jj, ptr;

    for (jj = 0; jj < P->n; jj++)
    {   // iterate over columns
        for (ptr = P->p[jj]; ptr < P->p[jj + 1]; ptr++)
        {   // iterate over rows
            ii = P->i[ptr]; // Row index

            if (ii == jj)
            {   // diagonal element
                quad_form += (float).5 * P->x[ptr] * x[ii] * x[ii];
            }
            else if (ii < jj)
            {   // off-diagonal element
                quad_form += P->x[ptr] * x[ii] * x[jj];
            }
            else { // element in lower diagonal
                // print("quad_form matrix is not upper triangular");
                return 0;
            }
        }
    }
    return quad_form;
}


