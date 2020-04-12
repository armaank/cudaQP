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




