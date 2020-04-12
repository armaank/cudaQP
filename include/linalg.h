#ifndef LINALG_H
# define LINALG_H

# inlcude "types.h"

/* vector functions */

/* copy vector a into preallocated vector b */
void prea_vec_copy(const c_float *a, c_float *b, c_int n);

/* copy integer vector a into preallocated vector b */
void prea_int_vec_copy(const c_int *a, c_int *b, c_int n);

/* set float vector to scalar */
void vec_set_scalar(c_float *a, c_float sc, c_int n);

/* set integer vector to scalar */
void int_vec_set_scalar(c_int *a, c_int sc, c_int n);

/* add scalar to vector*/
void vec_add_scalar(c_float *a, c_float sc, c_int n);

/* multiply scalar to vector */
void vec_mult_scalar(c_float *a, c_float sc, c_int n);

/* c = a + sc*b */
void vec_add_scaled(c_float *c, const c_float *a, const c_float *b, c_int n, c_float sc);

/* ||v||_inf */
c_float vec_norm_inf(const c_float *v, c_int l);

/* ||Sv||_inf */
c_float vec_scaled_norm_inf(const c_float *S, const c_float *v, c_int l);

/* ||a - b||_inf */
c_float vec_norm_inf_diff(const c_float *a, const c_float *b, c_int l);

/* mean of vector elements */
c_float vec_mean(const c_float *a, c_int n);

/* vector elementwise reciprocal b = 1./a (needed for scaling)*/
void vec_ew_recipr(const c_float *a, c_float *b, c_int n);

/* inner product a'b */
c_float vec_prod(const c_float *a, const c_float *b, c_int n);

/* elementwise product a.*b stored in c*/
void vec_ew_prod(const c_float *a, const c_float *b, c_float *c, c_int n);

/* elementwise sqrt of the vector elements */
void vec_ew_sqrt(c_float *a, c_int n);

/* elementwise max between each vector component and max_val */
void vec_ew_max(c_float *a, c_int n, c_float max_val);

/* elementwise min between each vector component and max_val */
void vec_ew_min(c_float *a, c_int n, c_float min_val);

/* elementwise maximum between vectors c = max(a, b) */
void vec_ew_max_vec(const c_float *a, const c_float *b, c_float *c, c_int n);

/* elementwise minimum between vectors c = min(a, b) */
void vec_ew_min_vec(const c_float *a, const c_float *b, c_float *c, c_int n);


/* matrix functons */

/* multiply matrix by a scalar */
void mat_mult_scalar(csc *A, c_float sc);

/* premultiply matrix A by diagonal matrix with diagonal d,
   i.e. scale the rows of A by d (thanks fred)
 */
void mat_premult_diag(csc *A, const c_float *d);

/* premultiply matrix A by diagonal matrix with diagonal d,
   i.e. scale the columns of A by d (thanks fred)
 */
void mat_postmult_diag(csc *A, const c_float *d);

/* matrix-vector multiplication
 *    y  =  A*x  (if plus_eq == 0)
 *    y +=  A*x  (if plus_eq == 1)
 *    y -=  A*x  (if plus_eq == -1)
 */
void mat_vec(const csc *A, const c_float *x, c_float *y, c_int plus_eq);


/* matrix-transpose-vector multiplication
 *    y  =  A'*x  (if plus_eq == 0)
 *    y +=  A'*x  (if plus_eq == 1)
 *    y -=  A'*x  (if plus_eq == -1)
 * If skip_diag == 1, then diagonal elements of A are assumed to be zero.
 */
void mat_tpose_vec(const csc *A, const c_float *x, c_float *y, c_int plus_eq, c_int skip_diag);

/* infinity norm of each matrix column */
void mat_inf_norm_cols(const csc *M, c_float *E);

/* infinity norm of each matrix row */
void mat_inf_norm_rows(const csc *M, c_float *E);

/* infinity norm of each matrix column
 * matrix M is symmetric upper-triangular
 */
void mat_inf_norm_cols_sym_triu(const csc *M, c_float *E);

/**
 * Compute quadratic form f(x) = 1/2 x' P x
 * matrix P is in CSC form (only upper triangular)
 */
c_float quad_form(const csc *P, const c_float *x);


#endif /* ifndef LINALG_H */

