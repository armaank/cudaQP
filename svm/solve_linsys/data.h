#ifndef SOLVE_LINSYS_DATA_H
#define SOLVE_LINSYS_DATA_H
#include "osqp.h"


/* create data and solutions structure */
typedef struct {
c_float * test_solve_KKT_rhs;
c_float * test_solve_KKT_x;
c_int test_solve_KKT_n;
c_int test_solve_KKT_m;
csc * test_solve_KKT_Pu;
csc * test_solve_KKT_A;
c_float test_solve_KKT_rho;
csc * test_solve_KKT_KKT;
c_float test_solve_KKT_sigma;
} solve_linsys_sols_data;

/* function prototypes */
solve_linsys_sols_data *  generate_problem_solve_linsys_sols_data();
void clean_problem_solve_linsys_sols_data(solve_linsys_sols_data * data);


/* function to define problem data */
solve_linsys_sols_data *  generate_problem_solve_linsys_sols_data(){

solve_linsys_sols_data * data = (solve_linsys_sols_data *)c_malloc(sizeof(solve_linsys_sols_data));

data->test_solve_KKT_rhs = (c_float*) c_malloc(7 * sizeof(c_float));
data->test_solve_KKT_rhs[0] = -1.05436800185014400988;
data->test_solve_KKT_rhs[1] = -1.00889150318170695009;
data->test_solve_KKT_rhs[2] = -0.06752199193321305193;
data->test_solve_KKT_rhs[3] = -0.08142504953145184021;
data->test_solve_KKT_rhs[4] = 2.83521598124578755318;
data->test_solve_KKT_rhs[5] = -0.24617516933298103088;
data->test_solve_KKT_rhs[6] = 0.70551086502222482011;
data->test_solve_KKT_x = (c_float*) c_malloc(7 * sizeof(c_float));
data->test_solve_KKT_x[0] = 1.59146823014194582768;
data->test_solve_KKT_x[1] = -0.39979419212238193060;
data->test_solve_KKT_x[2] = -0.51903227424323850059;
data->test_solve_KKT_x[3] = 0.55635674324897865795;
data->test_solve_KKT_x[4] = 0.86927159669526288255;
data->test_solve_KKT_x[5] = -0.38765340434320433305;
data->test_solve_KKT_x[6] = 0.00000000000000000000;
data->test_solve_KKT_n = 3;
data->test_solve_KKT_m = 4;

// Matrix test_solve_KKT_Pu
//-------------------------
data->test_solve_KKT_Pu = (csc*) c_malloc(sizeof(csc));
data->test_solve_KKT_Pu->m = 3;
data->test_solve_KKT_Pu->n = 3;
data->test_solve_KKT_Pu->nz = -1;
data->test_solve_KKT_Pu->nzmax = 2;
data->test_solve_KKT_Pu->x = (c_float*) c_malloc(2 * sizeof(c_float));
data->test_solve_KKT_Pu->x[0] = 0.38349652976488485256;
data->test_solve_KKT_Pu->x[1] = 0.15100215738121633424;
data->test_solve_KKT_Pu->i = (c_int*) c_malloc(2 * sizeof(c_int));
data->test_solve_KKT_Pu->i[0] = 0;
data->test_solve_KKT_Pu->i[1] = 1;
data->test_solve_KKT_Pu->p = (c_int*) c_malloc((3 + 1) * sizeof(c_int));
data->test_solve_KKT_Pu->p[0] = 0;
data->test_solve_KKT_Pu->p[1] = 1;
data->test_solve_KKT_Pu->p[2] = 2;
data->test_solve_KKT_Pu->p[3] = 2;


// Matrix test_solve_KKT_A
//------------------------
data->test_solve_KKT_A = (csc*) c_malloc(sizeof(csc));
data->test_solve_KKT_A->m = 4;
data->test_solve_KKT_A->n = 3;
data->test_solve_KKT_A->nz = -1;
data->test_solve_KKT_A->nzmax = 4;
data->test_solve_KKT_A->x = (c_float*) c_malloc(4 * sizeof(c_float));
data->test_solve_KKT_A->x[0] = 0.40730783228994515976;
data->test_solve_KKT_A->x[1] = 0.54620731990215798390;
data->test_solve_KKT_A->x[2] = 0.96963240582679277590;
data->test_solve_KKT_A->x[3] = 0.17698462366793654699;
data->test_solve_KKT_A->i = (c_int*) c_malloc(4 * sizeof(c_int));
data->test_solve_KKT_A->i[0] = 0;
data->test_solve_KKT_A->i[1] = 1;
data->test_solve_KKT_A->i[2] = 2;
data->test_solve_KKT_A->i[3] = 0;
data->test_solve_KKT_A->p = (c_int*) c_malloc((3 + 1) * sizeof(c_int));
data->test_solve_KKT_A->p[0] = 0;
data->test_solve_KKT_A->p[1] = 2;
data->test_solve_KKT_A->p[2] = 3;
data->test_solve_KKT_A->p[3] = 4;

data->test_solve_KKT_rho = 4.00000000000000000000;

// Matrix test_solve_KKT_KKT
//--------------------------
data->test_solve_KKT_KKT = (csc*) c_malloc(sizeof(csc));
data->test_solve_KKT_KKT->m = 7;
data->test_solve_KKT_KKT->n = 7;
data->test_solve_KKT_KKT->nz = -1;
data->test_solve_KKT_KKT->nzmax = 15;
data->test_solve_KKT_KKT->x = (c_float*) c_malloc(15 * sizeof(c_float));
data->test_solve_KKT_KKT->x[0] = 1.38349652976488490808;
data->test_solve_KKT_KKT->x[1] = 0.40730783228994515976;
data->test_solve_KKT_KKT->x[2] = 0.54620731990215798390;
data->test_solve_KKT_KKT->x[3] = 1.15100215738121636200;
data->test_solve_KKT_KKT->x[4] = 0.96963240582679277590;
data->test_solve_KKT_KKT->x[5] = 1.00000000000000000000;
data->test_solve_KKT_KKT->x[6] = 0.17698462366793654699;
data->test_solve_KKT_KKT->x[7] = 0.40730783228994515976;
data->test_solve_KKT_KKT->x[8] = 0.17698462366793654699;
data->test_solve_KKT_KKT->x[9] = -0.25000000000000000000;
data->test_solve_KKT_KKT->x[10] = 0.54620731990215798390;
data->test_solve_KKT_KKT->x[11] = -0.25000000000000000000;
data->test_solve_KKT_KKT->x[12] = 0.96963240582679277590;
data->test_solve_KKT_KKT->x[13] = -0.25000000000000000000;
data->test_solve_KKT_KKT->x[14] = -0.25000000000000000000;
data->test_solve_KKT_KKT->i = (c_int*) c_malloc(15 * sizeof(c_int));
data->test_solve_KKT_KKT->i[0] = 0;
data->test_solve_KKT_KKT->i[1] = 3;
data->test_solve_KKT_KKT->i[2] = 4;
data->test_solve_KKT_KKT->i[3] = 1;
data->test_solve_KKT_KKT->i[4] = 5;
data->test_solve_KKT_KKT->i[5] = 2;
data->test_solve_KKT_KKT->i[6] = 3;
data->test_solve_KKT_KKT->i[7] = 0;
data->test_solve_KKT_KKT->i[8] = 2;
data->test_solve_KKT_KKT->i[9] = 3;
data->test_solve_KKT_KKT->i[10] = 0;
data->test_solve_KKT_KKT->i[11] = 4;
data->test_solve_KKT_KKT->i[12] = 1;
data->test_solve_KKT_KKT->i[13] = 5;
data->test_solve_KKT_KKT->i[14] = 6;
data->test_solve_KKT_KKT->p = (c_int*) c_malloc((7 + 1) * sizeof(c_int));
data->test_solve_KKT_KKT->p[0] = 0;
data->test_solve_KKT_KKT->p[1] = 3;
data->test_solve_KKT_KKT->p[2] = 5;
data->test_solve_KKT_KKT->p[3] = 7;
data->test_solve_KKT_KKT->p[4] = 10;
data->test_solve_KKT_KKT->p[5] = 12;
data->test_solve_KKT_KKT->p[6] = 14;
data->test_solve_KKT_KKT->p[7] = 15;

data->test_solve_KKT_sigma = 1.00000000000000000000;

return data;

}

/* function to clean data struct */
void clean_problem_solve_linsys_sols_data(solve_linsys_sols_data * data){

c_free(data->test_solve_KKT_rhs);
c_free(data->test_solve_KKT_x);
c_free(data->test_solve_KKT_Pu->x);
c_free(data->test_solve_KKT_Pu->i);
c_free(data->test_solve_KKT_Pu->p);
c_free(data->test_solve_KKT_Pu);
c_free(data->test_solve_KKT_A->x);
c_free(data->test_solve_KKT_A->i);
c_free(data->test_solve_KKT_A->p);
c_free(data->test_solve_KKT_A);
c_free(data->test_solve_KKT_KKT->x);
c_free(data->test_solve_KKT_KKT->i);
c_free(data->test_solve_KKT_KKT->p);
c_free(data->test_solve_KKT_KKT);

c_free(data);

}

#endif