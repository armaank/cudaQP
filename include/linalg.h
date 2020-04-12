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


#endif /* ifndef LINALG_H */

