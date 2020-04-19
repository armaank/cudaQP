#include <stdio.h>
#include <stdlib.h>
#include "../../../include/unittest.h"
#include "../../../include/csc.h"

#include "data.h"
#include "../qptest.h"
#define TESTS_TOL 1e-4 // Define tests tolerance


static const char* test_constr_sparse_mat() {
    float *Adns; // Conversion to dense matrix

    lin_alg_sols_data *data = generate_problem_lin_alg_sols_data();

    // Convert sparse to dense
    Adns = csc_to_dns(data->test_sp_matrix_A);

    // Compute norm of the elementwise difference with
    ut_assert("Linear algebra tests: error in constructing sparse/dense matrix!",
              vec_norm_inf_diff(Adns, data->test_sp_matrix_Adns,
                                data->test_sp_matrix_A->m *
                                data->test_sp_matrix_A->n) < TESTS_TOL);

    // Free memory
    free(Adns); // because of vars from file matrices.h
    clean_problem_lin_alg_sols_data(data);

    return 0;
}

static const char* test_vec_operations() {
    float  norm_inf, vecprod; // normInf;
    float *ew_reciprocal;
    float *add_scaled;
    float *vec_ew_max_vec_test, *vec_ew_min_vec_test;

    lin_alg_sols_data *data = generate_problem_lin_alg_sols_data();


    // Add scaled
    add_scaled = vec_copy(data->test_vec_ops_v1, data->test_vec_ops_n);
    vec_add_scaled(add_scaled,
                   add_scaled,
                   data->test_vec_ops_v2,
                   data->test_vec_ops_n,
                   data->test_vec_ops_sc);
    ut_assert(
        "Linear algebra tests: error in vector operation, adding scaled vector",
        vec_norm_inf_diff(add_scaled, data->test_vec_ops_add_scaled,
                          data->test_vec_ops_n) < TESTS_TOL);

    // Norm_inf of the difference
    ut_assert(
        "Linear algebra tests: error in vector operation, norm_inf of difference",
        c_absval(vec_norm_inf_diff(data->test_vec_ops_v1,
                                   data->test_vec_ops_v2,
                                   data->test_vec_ops_n) -
                 data->test_vec_ops_norm_inf_diff) <
        TESTS_TOL);

    // norm_inf
    norm_inf = vec_norm_inf(data->test_vec_ops_v1, data->test_vec_ops_n);
    ut_assert("Linear algebra tests: error in vector operation, norm_inf",
              c_absval(norm_inf - data->test_vec_ops_norm_inf) < TESTS_TOL);

    // Elementwise reciprocal
    ew_reciprocal = (float *)malloc(data->test_vec_ops_n * sizeof(float));
    vec_ew_recipr(data->test_vec_ops_v1, ew_reciprocal, data->test_vec_ops_n);
    ut_assert(
        "Linear algebra tests: error in vector operation, elementwise reciprocal",
        vec_norm_inf_diff(ew_reciprocal, data->test_vec_ops_ew_reciprocal,
                          data->test_vec_ops_n) < TESTS_TOL);


    // Vector product
    vecprod = vec_prod(data->test_vec_ops_v1,
                       data->test_vec_ops_v2,
                       data->test_vec_ops_n);
    ut_assert("Linear algebra tests: error in vector operation, vector product",
              c_absval(vecprod - data->test_vec_ops_vec_prod) < TESTS_TOL);

    // Elementwise maximum between two vectors
    vec_ew_max_vec_test =
        (float *)malloc(data->test_vec_ops_n * sizeof(float));
    vec_ew_max_vec(data->test_vec_ops_v1,
                   data->test_vec_ops_v2,
                   vec_ew_max_vec_test,
                   data->test_vec_ops_n);
    ut_assert(
        "Linear algebra tests: error in vector operation, elementwise maximum between vectors",
        vec_norm_inf_diff(vec_ew_max_vec_test, data->test_vec_ops_ew_max_vec,
                          data
                          ->test_vec_ops_n) < TESTS_TOL);

    // Elementwise minimum between two vectors
    vec_ew_min_vec_test =
        (float *)malloc(data->test_vec_ops_n * sizeof(float));
    vec_ew_min_vec(data->test_vec_ops_v1,
                   data->test_vec_ops_v2,
                   vec_ew_min_vec_test,
                   data->test_vec_ops_n);
    ut_assert(
        "Linear algebra tests: error in vector operation, elementwise minimum between vectors",
        vec_norm_inf_diff(vec_ew_min_vec_test, data->test_vec_ops_ew_min_vec,
                          data
                          ->test_vec_ops_n) < TESTS_TOL);

    // cleanup
    free(add_scaled);
    free(ew_reciprocal);
    free(vec_ew_min_vec_test);
    free(vec_ew_max_vec_test);
    clean_problem_lin_alg_sols_data(data);

    return 0;
}

static const char* test_mat_operations() {
    csc *Ad, *dA; // Matrices used for tests
    // csc *A_ewsq, *A_ewabs;     // Matrices used for tests
    int exitflag = 0;

    // float trace, fro_sq;
    float *inf_norm_cols_rows_test;


    lin_alg_sols_data *data = generate_problem_lin_alg_sols_data();


    // Copy matrices
    Ad = copy_csc_mat(data->test_mat_ops_A);
    dA = copy_csc_mat(data->test_mat_ops_A);


    // Premultiply matrix A
    mat_premult_diag(dA, data->test_mat_ops_d);
    ut_assert(
        "Linear algebra tests: error in matrix operation, premultiply diagonal",
        is_eq_csc(dA, data->test_mat_ops_prem_diag, TESTS_TOL));


    // Postmultiply matrix A
    mat_postmult_diag(Ad, data->test_mat_ops_d);
    ut_assert(
        "Linear algebra tests: error in matrix operation, postmultiply diagonal",
        is_eq_csc(Ad, data->test_mat_ops_postm_diag, TESTS_TOL));

    // Maximum norm over columns
    inf_norm_cols_rows_test =
        (float *)malloc(data->test_mat_ops_n * sizeof(float));
    mat_inf_norm_cols(data->test_mat_ops_A, inf_norm_cols_rows_test);
    ut_assert(
        "Linear algebra tests: error in matrix operation, max norm over columns",
        vec_norm_inf_diff(inf_norm_cols_rows_test, data->test_mat_ops_inf_norm_cols,
                          data
                          ->test_mat_ops_n) < TESTS_TOL);

    // Maximum norm over rows
    mat_inf_norm_rows(data->test_mat_ops_A, inf_norm_cols_rows_test);
    ut_assert("Linear algebra tests: error in matrix operation, max norm over rows",
              vec_norm_inf_diff(inf_norm_cols_rows_test,
                                data->test_mat_ops_inf_norm_rows,
                                data
                                ->test_mat_ops_n) < TESTS_TOL);


    // cleanup
    free(inf_norm_cols_rows_test);
    csc_spfree(Ad);
    csc_spfree(dA);
    clean_problem_lin_alg_sols_data(data);

    return 0;
}

static const char* test_mat_vec_multiplication() {
    float *Ax, *ATy, *Px, *Ax_cum, *ATy_cum, *Px_cum;

    lin_alg_sols_data *data = generate_problem_lin_alg_sols_data();


    // Allocate vectors
    Ax  = (float *)malloc(data->test_mat_vec_m * sizeof(float));
    ATy = (float *)malloc(data->test_mat_vec_n * sizeof(float));
    Px  = (float *)malloc(data->test_mat_vec_n * sizeof(float));


    // Matrix-vector multiplication:  y = Ax
    mat_vec(data->test_mat_vec_A, data->test_mat_vec_x, Ax, 0);
    ut_assert(
        "Linear algebra tests: error in matrix-vector operation, matrix-vector multiplication",
        vec_norm_inf_diff(Ax, data->test_mat_vec_Ax,
                          data->test_mat_vec_m) < TESTS_TOL);

    // Cumulative matrix-vector multiplication:  y += Ax
    Ax_cum = vec_copy(data->test_mat_vec_y, data->test_mat_vec_m);
    mat_vec(data->test_mat_vec_A, data->test_mat_vec_x, Ax_cum, 1);
    ut_assert(
        "Linear algebra tests: error in matrix-vector operation, cumulative matrix-vector multiplication",
        vec_norm_inf_diff(Ax_cum, data->test_mat_vec_Ax_cum,
                          data->test_mat_vec_m) < TESTS_TOL);

    // Matrix-transpose-vector multiplication:  x = A'*y
    mat_tpose_vec(data->test_mat_vec_A, data->test_mat_vec_y, ATy, 0, 0);
    ut_assert(
        "Linear algebra tests: error in matrix-vector operation, matrix-transpose-vector multiplication",
        vec_norm_inf_diff(ATy, data->test_mat_vec_ATy,
                          data->test_mat_vec_n) < TESTS_TOL);

    // Cumulative matrix-transpose-vector multiplication:  x += A'*y
    ATy_cum = vec_copy(data->test_mat_vec_x, data->test_mat_vec_n);
    mat_tpose_vec(data->test_mat_vec_A, data->test_mat_vec_y, ATy_cum, 1, 0);
    ut_assert(
        "Linear algebra tests: error in matrix-vector operation, cumulative matrix-transpose-vector multiplication",
        vec_norm_inf_diff(ATy_cum, data->test_mat_vec_ATy_cum,
                          data->test_mat_vec_n) < TESTS_TOL);

    // Symmetric-matrix-vector multiplication (only upper part is stored)
    mat_vec(data->test_mat_vec_Pu, data->test_mat_vec_x, Px, 0);          // upper
    // traingular
    // part
    mat_tpose_vec(data->test_mat_vec_Pu, data->test_mat_vec_x, Px, 1, 1); // lower
    // traingular
    // part
    // (without
    // diagonal)
    ut_assert(
        "Linear algebra tests: error in matrix-vector operation, symmetric matrix-vector multiplication",
        vec_norm_inf_diff(Px, data->test_mat_vec_Px,
                          data->test_mat_vec_n) < TESTS_TOL);


    // Cumulative symmetric-matrix-vector multiplication
    Px_cum = vec_copy(data->test_mat_vec_x, data->test_mat_vec_n);
    mat_vec(data->test_mat_vec_Pu, data->test_mat_vec_x, Px_cum, 1);          // upper
    // traingular
    // part
    mat_tpose_vec(data->test_mat_vec_Pu, data->test_mat_vec_x, Px_cum, 1, 1); // lower
    // traingular
    // part
    // (without
    // diagonal)
    ut_assert(
        "Linear algebra tests: error in matrix-vector operation, cumulative symmetric matrix-vector multiplication",
        vec_norm_inf_diff(Px_cum, data->test_mat_vec_Px_cum,
                          data->test_mat_vec_n) < TESTS_TOL);


    // cleanup
    free(Ax);
    free(ATy);
    free(Px);
    free(Ax_cum);
    free(ATy_cum);
    free(Px_cum);
    clean_problem_lin_alg_sols_data(data);

    return 0;
}

static const char* test_extract_upper_triangular() {
    float *inf_norm_cols_test;
    lin_alg_sols_data *data = generate_problem_lin_alg_sols_data();

    // Extract upper triangular part
    csc *Ptriu = csc_to_triu(data->test_mat_extr_triu_P);

    ut_assert("Linear algebra tests: error in forming upper triangular matrix!",
              is_eq_csc(data->test_mat_extr_triu_Pu, Ptriu, TESTS_TOL));

    // Compute infinity norm over columns of the original matrix by using the
    // upper triangular part only
    inf_norm_cols_test = (float *)malloc(data->test_mat_extr_triu_n
                                         * sizeof(float));
    mat_inf_norm_cols_sym_triu(Ptriu, inf_norm_cols_test);
    ut_assert(
        "Linear algebra tests: error in forming upper triangular matrix, infinity norm over columns",
        vec_norm_inf_diff(inf_norm_cols_test,
                          data->test_mat_extr_triu_P_inf_norm_cols,
                          data->test_mat_extr_triu_n) < TESTS_TOL);

    // Cleanup
    free(inf_norm_cols_test);
    csc_spfree(Ptriu);
    clean_problem_lin_alg_sols_data(data);

    return 0;
}

static const char* test_quad_form_upper_triang() {
    float quad_form_t;

    lin_alg_sols_data *data = generate_problem_lin_alg_sols_data();

    // Compute quadratic form
    quad_form_t = quad_form(data->test_qpform_Pu, data->test_qpform_x);

    ut_assert(
        "Linear algebra tests: error in computing quadratic form using upper triangular matrix!",
        (c_absval(quad_form_t - data->test_qpform_value) < TESTS_TOL));

    // cleanup
    clean_problem_lin_alg_sols_data(data);

    return 0;
}

static const char* test_lin_alg()
{
    ut_run_test(test_constr_sparse_mat);
    ut_run_test(test_vec_operations);
    ut_run_test(test_mat_operations);
    ut_run_test(test_mat_vec_multiplication);
    ut_run_test(test_extract_upper_triangular);
    ut_run_test(test_quad_form_upper_triang);

    return 0;
}
