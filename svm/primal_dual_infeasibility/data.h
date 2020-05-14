#ifndef PRIMAL_DUAL_INFEASIBILITY_DATA_H
#define PRIMAL_DUAL_INFEASIBILITY_DATA_H
#include "osqp.h"


/* create data and solutions structure */
typedef struct {
c_int status4;
c_int status1;
c_int status3;
csc * A12;
c_float * u4;
csc * A34;
c_float * u1;
c_float * l;
c_float * u3;
c_float * u2;
c_float * q;
csc * P;
c_float obj_value1;
c_float * y1;
c_float * x1;
c_int status2;
} primal_dual_infeasibility_sols_data;

/* function prototypes */
primal_dual_infeasibility_sols_data *  generate_problem_primal_dual_infeasibility_sols_data();
void clean_problem_primal_dual_infeasibility_sols_data(primal_dual_infeasibility_sols_data * data);


/* function to define problem data */
primal_dual_infeasibility_sols_data *  generate_problem_primal_dual_infeasibility_sols_data(){

primal_dual_infeasibility_sols_data * data = (primal_dual_infeasibility_sols_data *)c_malloc(sizeof(primal_dual_infeasibility_sols_data));

data->status4 = OSQP_PRIMAL_INFEASIBLE;
data->status1 = OSQP_SOLVED;
data->status3 = OSQP_DUAL_INFEASIBLE;

// Matrix A12
//-----------
data->A12 = (csc*) c_malloc(sizeof(csc));
data->A12->m = 3;
data->A12->n = 2;
data->A12->nz = -1;
data->A12->nzmax = 4;
data->A12->x = (c_float*) c_malloc(4 * sizeof(c_float));
data->A12->x[0] = 1.00000000000000000000;
data->A12->x[1] = 1.00000000000000000000;
data->A12->x[2] = 1.00000000000000000000;
data->A12->x[3] = 1.00000000000000000000;
data->A12->i = (c_int*) c_malloc(4 * sizeof(c_int));
data->A12->i[0] = 0;
data->A12->i[1] = 1;
data->A12->i[2] = 0;
data->A12->i[3] = 2;
data->A12->p = (c_int*) c_malloc((2 + 1) * sizeof(c_int));
data->A12->p[0] = 0;
data->A12->p[1] = 2;
data->A12->p[2] = 4;

data->u4 = (c_float*) c_malloc(3 * sizeof(c_float));
data->u4[0] = 0.00000000000000000000;
data->u4[1] = 3.00000000000000000000;
data->u4[2] = OSQP_INFTY;

// Matrix A34
//-----------
data->A34 = (csc*) c_malloc(sizeof(csc));
data->A34->m = 3;
data->A34->n = 2;
data->A34->nz = -1;
data->A34->nzmax = 3;
data->A34->x = (c_float*) c_malloc(3 * sizeof(c_float));
data->A34->x[0] = 1.00000000000000000000;
data->A34->x[1] = 1.00000000000000000000;
data->A34->x[2] = 1.00000000000000000000;
data->A34->i = (c_int*) c_malloc(3 * sizeof(c_int));
data->A34->i[0] = 0;
data->A34->i[1] = 1;
data->A34->i[2] = 2;
data->A34->p = (c_int*) c_malloc((2 + 1) * sizeof(c_int));
data->A34->p[0] = 0;
data->A34->p[1] = 2;
data->A34->p[2] = 3;

data->u1 = (c_float*) c_malloc(3 * sizeof(c_float));
data->u1[0] = 5.00000000000000000000;
data->u1[1] = 3.00000000000000000000;
data->u1[2] = 3.00000000000000000000;
data->l = (c_float*) c_malloc(3 * sizeof(c_float));
data->l[0] = 0.00000000000000000000;
data->l[1] = 1.00000000000000000000;
data->l[2] = 1.00000000000000000000;
data->u3 = (c_float*) c_malloc(3 * sizeof(c_float));
data->u3[0] = 2.00000000000000000000;
data->u3[1] = 3.00000000000000000000;
data->u3[2] = OSQP_INFTY;
data->u2 = (c_float*) c_malloc(3 * sizeof(c_float));
data->u2[0] = 0.00000000000000000000;
data->u2[1] = 3.00000000000000000000;
data->u2[2] = 3.00000000000000000000;
data->q = (c_float*) c_malloc(2 * sizeof(c_float));
data->q[0] = 1.00000000000000000000;
data->q[1] = -1.00000000000000000000;

// Matrix P
//---------
data->P = (csc*) c_malloc(sizeof(csc));
data->P->m = 2;
data->P->n = 2;
data->P->nz = -1;
data->P->nzmax = 1;
data->P->x = (c_float*) c_malloc(1 * sizeof(c_float));
data->P->x[0] = 1.00000000000000000000;
data->P->i = (c_int*) c_malloc(1 * sizeof(c_int));
data->P->i[0] = 0;
data->P->p = (c_int*) c_malloc((2 + 1) * sizeof(c_int));
data->P->p[0] = 0;
data->P->p[1] = 1;
data->P->p[2] = 1;

data->obj_value1 = -1.50000000000000000000;
data->y1 = (c_float*) c_malloc(3 * sizeof(c_float));
data->y1[0] = 0.00000000000000000000;
data->y1[1] = -2.00000000000000000000;
data->y1[2] = 1.00000000000000000000;
data->x1 = (c_float*) c_malloc(2 * sizeof(c_float));
data->x1[0] = 1.00000000000000000000;
data->x1[1] = 3.00000000000000000000;
data->status2 = OSQP_PRIMAL_INFEASIBLE;

return data;

}

/* function to clean data struct */
void clean_problem_primal_dual_infeasibility_sols_data(primal_dual_infeasibility_sols_data * data){

c_free(data->A12->x);
c_free(data->A12->i);
c_free(data->A12->p);
c_free(data->A12);
c_free(data->u4);
c_free(data->A34->x);
c_free(data->A34->i);
c_free(data->A34->p);
c_free(data->A34);
c_free(data->u1);
c_free(data->l);
c_free(data->u3);
c_free(data->u2);
c_free(data->q);
c_free(data->P->x);
c_free(data->P->i);
c_free(data->P->p);
c_free(data->P);
c_free(data->y1);
c_free(data->x1);

c_free(data);

}

#endif
