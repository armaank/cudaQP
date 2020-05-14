#ifndef BASIC_SVM_DATA_H
#define BASIC_SVM_DATA_H
#include "osqp.h"


/* create additional data and solutions structure */
typedef struct {
c_float * x_test;
c_float obj_value_test;
c_float * q_new;
c_float obj_value_test_new;
c_float * y_test_new;
c_int status_test_new;
c_int status_test;
c_float * x_test_new;
c_float * y_test;
c_float * u_new;
} basic_svm_sols_data;

/* function prototypes */
OSQPTestData * generate_problem_basic_svm();
void clean_problem_basic_svm(OSQPTestData * data);
basic_svm_sols_data *  generate_problem_basic_svm_sols_data();
void clean_problem_basic_svm_sols_data(basic_svm_sols_data * data);


/* function to generate QP problem data */
OSQPTestData * generate_problem_basic_svm(){

OSQPTestData * data = (OSQPTestData *)c_malloc(sizeof(OSQPTestData));

// Problem dimensions
data->n = 2;
data->m = 5;

// Problem vectors
data->l = (c_float*) c_malloc(5 * sizeof(c_float));
data->l[0] = -OSQP_INFTY;
data->l[1] = -OSQP_INFTY;
data->l[2] = -OSQP_INFTY;
data->l[3] = -OSQP_INFTY;
data->l[4] = -OSQP_INFTY;
data->u = (c_float*) c_malloc(5 * sizeof(c_float));
data->u[0] = 0.00000000000000000000;
data->u[1] = 0.00000000000000000000;
data->u[2] = -15.00000000000000000000;
data->u[3] = 100.00000000000000000000;
data->u[4] = 80.00000000000000000000;
data->q = (c_float*) c_malloc(2 * sizeof(c_float));
data->q[0] = 3.00000000000000000000;
data->q[1] = 4.00000000000000000000;


// Matrix A
//---------
data->A = (csc*) c_malloc(sizeof(csc));
data->A->m = 5;
data->A->n = 2;
data->A->nz = -1;
data->A->nzmax = 8;
data->A->x = (c_float*) c_malloc(8 * sizeof(c_float));
data->A->x[0] = -1.00000000000000000000;
data->A->x[1] = -1.00000000000000000000;
data->A->x[2] = 2.00000000000000000000;
data->A->x[3] = 3.00000000000000000000;
data->A->x[4] = -1.00000000000000000000;
data->A->x[5] = 3.00000000000000000000;
data->A->x[6] = 5.00000000000000000000;
data->A->x[7] = 4.00000000000000000000;
data->A->i = (c_int*) c_malloc(8 * sizeof(c_int));
data->A->i[0] = 0;
data->A->i[1] = 2;
data->A->i[2] = 3;
data->A->i[3] = 4;
data->A->i[4] = 1;
data->A->i[5] = 2;
data->A->i[6] = 3;
data->A->i[7] = 4;
data->A->p = (c_int*) c_malloc((2 + 1) * sizeof(c_int));
data->A->p[0] = 0;
data->A->p[1] = 4;
data->A->p[2] = 8;


// Matrix P
//---------
data->P = (csc*) c_malloc(sizeof(csc));
data->P->m = 2;
data->P->n = 2;
data->P->nz = -1;
data->P->nzmax = 1;
data->P->x = (c_float*) c_malloc(1 * sizeof(c_float));
data->P->x[0] = 11.00000000000000000000;
data->P->i = (c_int*) c_malloc(1 * sizeof(c_int));
data->P->i[0] = 0;
data->P->p = (c_int*) c_malloc((2 + 1) * sizeof(c_int));
data->P->p[0] = 0;
data->P->p[1] = 1;
data->P->p[2] = 1;

return data;

}

/* function to clean problem data structure */
void clean_problem_basic_svm(OSQPTestData * data){

// Clean vectors
c_free(data->l);
c_free(data->u);
c_free(data->q);

//Clean Matrices
c_free(data->A->x);
c_free(data->A->i);
c_free(data->A->p);
c_free(data->A);
c_free(data->P->x);
c_free(data->P->i);
c_free(data->P->p);
c_free(data->P);

c_free(data);

}

/* function to define solutions and additional data struct */
basic_svm_sols_data *  generate_problem_basic_svm_sols_data(){

basic_svm_sols_data * data = (basic_svm_sols_data *)c_malloc(sizeof(basic_svm_sols_data));

data->x_test = (c_float*) c_malloc(2 * sizeof(c_float));
data->x_test[0] = 15.00000000000000000000;
data->x_test[1] = -0.00000000000000000000;
data->obj_value_test = 1282.50000000000000000000;
data->q_new = (c_float*) c_malloc(2 * sizeof(c_float));
data->q_new[0] = 1.00000000000000000000;
data->q_new[1] = 1.00000000000000000000;
data->obj_value_test_new = 2220.00000000000000000000;
data->y_test_new = (c_float*) c_malloc(5 * sizeof(c_float));
data->y_test_new[0] = 0.00000000000000000000;
data->y_test_new[1] = 664.00000000000000000000;
data->y_test_new[2] = 221.00000000000000000000;
data->y_test_new[3] = 0.00000000000000000000;
data->y_test_new[4] = 0.00000000000000000000;
data->status_test_new = OSQP_SOLVED;
data->status_test = OSQP_SOLVED;
data->x_test_new = (c_float*) c_malloc(2 * sizeof(c_float));
data->x_test_new[0] = 20.00000000000000000000;
data->x_test_new[1] = -0.00000000000000000000;
data->y_test = (c_float*) c_malloc(5 * sizeof(c_float));
data->y_test[0] = 0.00000000000000000000;
data->y_test[1] = 508.00000000000000000000;
data->y_test[2] = 168.00000000000000000000;
data->y_test[3] = 0.00000000000000000000;
data->y_test[4] = 0.00000000000000000000;
data->u_new = (c_float*) c_malloc(5 * sizeof(c_float));
data->u_new[0] = -2.00000000000000000000;
data->u_new[1] = 0.00000000000000000000;
data->u_new[2] = -20.00000000000000000000;
data->u_new[3] = 100.00000000000000000000;
data->u_new[4] = 80.00000000000000000000;

return data;

}

/* function to clean solutions and additional data struct */
void clean_problem_basic_svm_sols_data(basic_svm_sols_data * data){

c_free(data->x_test);
c_free(data->q_new);
c_free(data->y_test_new);
c_free(data->x_test_new);
c_free(data->y_test);
c_free(data->u_new);

c_free(data);

}

#endif
