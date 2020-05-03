#ifndef CUDA_LIN_ALG_H
# define CUDA_LIN_ALG_H

// #include "algebra_types.h" /
// todo figure out inlcudes, need csr


/*******************************************************************************
 *                           Vector Functions                                  *
 *******************************************************************************/

/*
 * d_y[i] = d_x[i] for i in [0,n-1]
*/
void cuda_vec_copy_d2d(float       *d_y,
                       const float *d_x,
                       int          n);

/*
 * d_y[i] = h_x[i] for i in [0,n-1]
*/
void cuda_vec_copy_h2d(float       *d_y,
                       const float *h_x,
                       int          n);

/*
 * h_y[i] = d_x[i] for i in [0,n-1]
*/
void cuda_vec_copy_d2h(float       *h_y,
                       const float *d_x,
                       int          n);

/*
 * d_y[i] = h_x[i] for i in [0,n-1] (integers)
*/
void cuda_veint_copy_h2d(int       *d_y,
                           const int *h_x,
                           int        n);

/*
 * h_y[i] = d_x[i] for i in [0,n-1] (integers)
*/
void cuda_veint_copy_d2h(int       *h_y,
                           const int *d_x,
                           int        n);

/**
 * d_a[i] = sc for i in [0,n-1]
 */
void cuda_vec_set_sc(float *d_a,
                     float  sc,
                     int    n);

/**
 *           | sc_if_neg   d_test[i]  < 0
 * d_a[i] = <  sc_if_zero  d_test[i] == 0   for i in [0,n-1]
 *           | sc_if_pos   d_test[i]  > 0
 */
void cuda_vec_set_sc_cond(float     *d_a,
                          const int *d_test,
                          float      sc_if_neg,
                          float      sc_if_zero,
                          float      sc_if_pos,
                          float      n);

/**
 * d_a[i] *= sc for i in [0,n-1]
 */
void cuda_vec_mult_sc(float *d_a,
                      float  sc,
                      int    n);

/**
 * d_x[i] = sca * d_a[i] + scb * d_b[i] for i in [0,n-1]
 */
void cuda_vec_add_scaled(float       *d_x,
                         const float *d_a,
                         const float *d_b,
                         float        sca,
                         float        scb,
                         int          n);

/**
 * d_x[i] = sca * d_a[i] + scb * d_b[i] + scc * d_c[i] for i in [0,n-1]
 */
void cuda_vec_add_scaled3(float       *d_x,
                          const float *d_a,
                          const float *d_b,
                          const float *d_c,
                          float        sca,
                          float        scb,
                          float        scc,
                          int          n);

/**
 * h_res = |d_x|_inf
 */
void cuda_vec_norm_inf(const float *d_x,
                       int          n,
                       float       *h_res);

/**
 * h_res = |d_x|_1
 */
void cuda_vec_norm_1(const float *d_x,
                     int          n,
                     float       *h_res);

/**
 * res = |d_x|_2
 */
void cuda_vec_norm_2(const float *d_x,
                     int          n,
                     float       *h_res);

/**
 * h_res = |S*v|_inf
 */
void cuda_vec_scaled_norm_inf(const float *d_S,
                              const float *d_v,
                              int          n,
                              float       *h_res);

/**
 * h_res = |d_a - d_b|_inf
 */
void cuda_vec_diff_norm_inf(const float *d_a,
                            const float *d_b,
                            int          n,
                            float       *h_res);

/**
 * h_res = sum(d_x) / n
 */
void cuda_vec_mean(const float *d_x,
                   int          n,
                   float       *h_res);

/**
 * h_res = d_a' * d_b
 */
void cuda_vec_prod(const float *d_a,
                   const float *d_b,
                   int          n,
                   float       *h_res);

/**
 *          | d_a' * max(d_b, 0)  sign ==  1
 * h_res = <  d_a' * min(d_b, 0)  sign == -1
 *          | d_a' * d_b          otherwise
 */
void cuda_vec_prod_signed(const float *d_a,
                          const float *d_b,
                          int          sign,
                          int          n,
                          float       *h_res);

/**
 * d_c[i] = d_a[i] * d_b[i] for i in [0,n-1]
 */
void cuda_vec_ew_prod(float       *d_c,
                      const float *d_a,
                      const float *d_b,
                      int          n);

/**
 * h_res = all(d_l <= d_u)
 */
void cuda_vec_leq(const float *d_l,
                  const float *d_u,
                  int          n,
                  int         *h_res);

/**
 * d_x[i] = min( max(d_z[i], d_l[i]), d_u[i] ) for i in [0,n-1]
 */
void cuda_vec_bound(float       *d_x,
                    const float *d_z,
                    const float *d_l,
                    const float *d_u,
                    int          n);

/**
 *           | 0.0               d_l < -infval AND d_u > +infval
 * d_y[i] = <  min(d_y[i], 0.0)  d_u > +infval
 *           | max(d_y[i], 0.0)  d_l < -infval
 */
void cuda_vec_project_polar_reccone(float       *d_y,
                                    const float *d_l,
                                    const float *d_u,
                                    float        infval,
                                    int          n);

/**
 *          | d_y[i] \in [-tol,tol]  d_l[i] > -infval AND d_u[i] < +infval
 * h_res = <  d_y[i] < +tol          d_l[i] < -infval AND d_u[i] < +infval
 *          | d_y[i] > -tol          d_l[i] > -infval AND d_u[i] > +infval
 */
void cuda_vec_in_reccone(const float *d_y,
                         const float *d_l,
                         const float *d_u,
                         float        infval,
                         float        tol,
                         int          n,
                         int         *h_res);

/**
 * d_b[i] = 1 / d_a[i] for i in [0,n-1]
 */
void cuda_vec_reciprocal(float       *d_b,
                         const float *d_a,
                         int          n);

/**
 * d_a[i] = sqrt(d_a[i]) for i in [0,n-1]
 */
void cuda_vec_sqrt(float *d_a,
                      int    n);

/**
 * d_c[i] = max(d_a[i], d_b[i]) for i in [0,n-1]
 */
void cuda_vec_max(float       *d_c,
                  const float *d_a,
                  const float *d_b,
                  int          n);

/**
 * d_c[i] = min(d_a[i], d_b[i]) for i in [0,n-1]
 */
void cuda_vec_min(float       *d_c,
                  const float *d_a,
                  const float *d_b,
                  int          n);

void cuda_vec_bounds_type(int         *d_iseq,
                          const float *d_l,
                          const float *d_u,
                          float        infval,
                          float        tol,
                          int          n,
                          int         *h_has_changed);

void cuda_vec_set_sc_if_lt(float       *d_x,
                           const float *d_z,
                           float        testval,
                           float        newval,
                           int          n);

void cuda_vec_set_sc_if_gt(float       *d_x,
                           const float *d_z,
                           float        testval,
                           float        newval,
                           int          n);

void cuda_vec_segmented_sum(const float *d_values,
                            const int   *d_keys,
                            float       *d_res,
                            void          *d_buffer,
                            int          num_segments,
                            int          num_elements);


/*******************************************************************************
 *                           Matrix Functions                                  *
 *******************************************************************************/

/**
 * S = sc * S
 */
void cuda_mat_mult_sc(csr     *S,
                      csr     *At,
                      int    symmetric,
                      float  sc);

/**
 * S = D * S
 */
void cuda_mat_lmult_diag(csr           *S,
                         csr           *At,
                         int          symmetric,
                         const float *d_diag);

/**
 * S = S * D
 */
void cuda_mat_rmult_diag(csr           *S,
                         csr           *At,
                         int          symmetric,
                         const float *d_diag);

/**
 * X = S * D
 * X->val values are stored in d_buffer.
 */
void cuda_mat_rmult_diag_new(const csr     *S,
                             float       *d_buffer,
                             const float *d_diag);

/**
 * d_y = alpha * A*d_x + beta*d_y
 */
void cuda_mat_Axpy(const csr     *A,
                   const float *d_x,
                   float       *d_y,
                   float        alpha,
                   float        beta);

/**
 * h_res = (1/2) d_x' * P * d_x
 */
void cuda_mat_quad_form(const csr     *P,
                        const float *d_x,
                        float       *h_res);

/**
 * d_res[i] = |S_i|_inf where S_i is i-th row of S
 */
void cuda_mat_row_norm_inf(const csr *S,
                           float   *d_res);



#endif /* ifndef CUDA_LIN_ALG_H */