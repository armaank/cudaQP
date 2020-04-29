// #ifndef BASIC_QP_DATA_H
// #define BASIC_QP_DATA_H
// #include "osqp.h"


/* create additional data and solutions structure */
typedef struct {
float * x_test;
float * y_test;
float obj_value_test;
int status_test;
float * q_new;
float * l_new;
float * u_new;
} basic_qp_sols_data;

/* function prototypes */
qpData * generate_problem_basic_qp();
void clean_problem_basic_qp(qpData * data);
basic_qp_sols_data *  generate_problem_basic_qp_sols_data();
void clean_problem_basic_qp_sols_data(basic_qp_sols_data * data);


/* function to generate QP problem data */
qpData * generate_problem_basic_qp(){

qpData * data = (qpData *)malloc(sizeof(qpData));

// Problem dimensions
data->n = 2;
data->m = 4;

// Problem vectors
data->l = (float*) malloc(4 * sizeof(float));
data->l[0] = 1.00000000000000000000;
data->l[1] = 0.00000000000000000000;
data->l[2] = 0.00000000000000000000;
data->l[3] = -qp_INFTY;
data->u = (float*) malloc(4 * sizeof(float));
data->u[0] = 1.00000000000000000000;
data->u[1] = 0.69999999999999995559;
data->u[2] = 0.69999999999999995559;
data->u[3] = qp_INFTY;
data->q = (float*) malloc(2 * sizeof(float));
data->q[0] = 1.00000000000000000000;
data->q[1] = 1.00000000000000000000;


// Matrix A
//---------
data->A = (csc*) malloc(sizeof(csc));
data->A->m = 4;
data->A->n = 2;
data->A->nz = -1;
data->A->nzmax = 5;
data->A->x = (float*) malloc(5 * sizeof(float));
data->A->x[0] = 1.00000000000000000000;
data->A->x[1] = 1.00000000000000000000;
data->A->x[2] = 1.00000000000000000000;
data->A->x[3] = 1.00000000000000000000;
data->A->x[4] = 1.00000000000000000000;
data->A->i = (int*) malloc(5 * sizeof(int));
data->A->i[0] = 0;
data->A->i[1] = 1;
data->A->i[2] = 0;
data->A->i[3] = 2;
data->A->i[4] = 3;
data->A->p = (int*) malloc((2 + 1) * sizeof(int));
data->A->p[0] = 0;
data->A->p[1] = 2;
data->A->p[2] = 5;


// Matrix P
//---------
data->P = (csc*) malloc(sizeof(csc));
data->P->m = 2;
data->P->n = 2;
data->P->nz = -1;
data->P->nzmax = 3;
data->P->x = (float*) malloc(3 * sizeof(float));
data->P->x[0] = 4.00000000000000000000;
data->P->x[1] = 1.00000000000000000000;
data->P->x[2] = 2.00000000000000000000;
data->P->i = (int*) malloc(3 * sizeof(int));
data->P->i[0] = 0;
data->P->i[1] = 0;
data->P->i[2] = 1;
data->P->p = (int*) malloc((2 + 1) * sizeof(int));
data->P->p[0] = 0;
data->P->p[1] = 1;
data->P->p[2] = 3;

return data;

}

/* function to clean problem data structure */
void clean_problem_basic_qp(qpData * data){

// Clean vectors
free(data->l);
free(data->u);
free(data->q);

//Clean Matrices
free(data->A->x);
free(data->A->i);
free(data->A->p);
free(data->A);
free(data->P->x);
free(data->P->i);
free(data->P->p);
free(data->P);

free(data);

}

/* function to define solutions and additional data struct */
basic_qp_sols_data *  generate_problem_basic_qp_sols_data(){

basic_qp_sols_data * data = (basic_qp_sols_data *)malloc(sizeof(basic_qp_sols_data));

data->x_test = (float*) malloc(2 * sizeof(float));
data->x_test[0] = 0.29999999999999998890;
data->x_test[1] = 0.69999999999999995559;
data->y_test = (float*) malloc(4 * sizeof(float));
data->y_test[0] = -2.89999999999999991118;
data->y_test[1] = 0.00000000000000000000;
data->y_test[2] = 0.20000000000000001110;
data->y_test[3] = 0.00000000000000000000;
data->obj_value_test = 1.87999999999999989342;
data->status_test = qp_SOLVED;
data->q_new = (float*) malloc(2 * sizeof(float));
data->q_new[0] = 2.50000000000000000000;
data->q_new[1] = 3.20000000000000017764;
data->l_new = (float*) malloc(4 * sizeof(float));
data->l_new[0] = 0.80000000000000004441;
data->l_new[1] = -3.39999999999999991118;
data->l_new[2] = -qp_INFTY;
data->l_new[3] = 0.50000000000000000000;
data->u_new = (float*) malloc(4 * sizeof(float));
data->u_new[0] = 1.60000000000000008882;
data->u_new[1] = 1.00000000000000000000;
data->u_new[2] = qp_INFTY;
data->u_new[3] = 0.50000000000000000000;

return data;

}

/* function to clean solutions and additional data struct */
void clean_problem_basic_qp_sols_data(basic_qp_sols_data * data){

free(data->x_test);
free(data->y_test);
free(data->q_new);
free(data->l_new);
free(data->u_new);

free(data);

}

