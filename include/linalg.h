#include <stdlib.h>
#include <math.h>
#include "ops.h"
#include "types.h"

/* vector functions */
/* uses malloc */
float* vec_copy(float *a,
                  int    n);
/* copy vector a into preallocated vector b */
void prea_vec_copy(const float *a, float *b, int n);

/* copy integer vector a into preallocated vector b */
void prea_int_vec_copy(const int *a, int *b, int n);

/* set float vector to scalar */
void vec_set_scalar(float *a, float sc, int n);

/* set integer vector to scalar */
void int_vec_set_scalar(int *a, int sc, int n);

/* add scalar to vector*/
void vec_add_scalar(float *a, float sc, int n);

/* multiply scalar to vector */
void vec_mult_scalar(float *a, float sc, int n);

/* c = a + sc*b */
void vec_add_scaled(float *c, const float *a, const float *b, int n, float sc);

/* ||v||_inf */
float vec_norm_inf(const float *v, int l);

/* ||Sv||_inf */
float vec_scaled_norm_inf(const float *S, const float *v, int l);

/* ||a - b||_inf */
float vec_norm_inf_diff(const float *a, const float *b, int l);

/* mean of vector elements */
float vec_mean(const float *a, int n);

/* vector elementwise reciprocal b = 1./a (needed for scaling)*/
void vec_ew_recipr(const float *a, float *b, int n);

/* inner product a'b */
float vec_prod(const float *a, const float *b, int n);

/* elementwise product a.*b stored in c*/
void vec_ew_prod(const float *a, const float *b, float *c, int n);

/* elementwise sqrt of the vector elements */
void vec_ew_sqrt(float *a, int n);

/* elementwise max between each vector component and max_val */
void vec_ew_max(float *a, int n, float max_val);

/* elementwise min between each vector component and max_val */
void vec_ew_min(float *a, int n, float min_val);

/* elementwise maximum between vectors c = max(a, b) */
void vec_ew_max_vec(const float *a, const float *b, float *c, int n);

/* elementwise minimum between vectors c = min(a, b) */
void vec_ew_min_vec(const float *a, const float *b, float *c, int n);


/* matrix functons */

/* multiply matrix by a scalar */
void mat_mult_scalar(csc *A, float sc);

/* premultiply matrix A by diagonal matrix with diagonal d,
   i.e. scale the rows of A by d (thanks fred)
 */
void mat_premult_diag(csc *A, const float *d);

/* premultiply matrix A by diagonal matrix with diagonal d,
   i.e. scale the columns of A by d (thanks fred)
 */
void mat_postmult_diag(csc *A, const float *d);

/* matrix-vector multiplication
 *    y  =  A*x  (if plus_eq == 0)
 *    y +=  A*x  (if plus_eq == 1)
 *    y -=  A*x  (if plus_eq == -1)
 */
void mat_vec(const csc *A, const float *x, float *y, int plus_eq);


/* matrix-transpose-vector multiplication
 *    y  =  A'*x  (if plus_eq == 0)
 *    y +=  A'*x  (if plus_eq == 1)
 *    y -=  A'*x  (if plus_eq == -1)
 * If skip_diag == 1, then diagonal elements of A are assumed to be zero.
 */
void mat_tpose_vec(const csc *A, const float *x, float *y, int plus_eq, int skip_diag);

/* infinity norm of each matrix column */
void mat_inf_norm_cols(const csc *M, float *E);

/* infinity norm of each matrix row */
void mat_inf_norm_rows(const csc *M, float *E);

/* infinity norm of each matrix column
 * matrix M is symmetric upper-triangular
 */
void mat_inf_norm_cols_sym_triu(const csc *M, float *E);

/**
 * Compute quadratic form f(x) = 1/2 x' P x
 * matrix P is in CSC form (only upper triangular)
 */
float quad_form(const csc *P, const float *x);


