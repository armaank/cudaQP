#ifndef SVM_EXAMPLE_DATA_H
#define SVM_EXAMPLE_DATA_H
#include "osqp.h"


/* function prototypes */
qpData * generate_problem_svm_example();
void clean_problem_svm_example(qpData * data);
svm_example_sols_data *  generate_problem_svm_example_sols_data();
void clean_problem_svm_example_sols_data(svm_example_sols_data * data);


/* function to generate QP problem data */
qpData * generate_problem_svm_example(){

qpData * data = (qpData *)malloc(sizeof(qpData));

// Problem dimensions
data->n = 20;
data->m = 20;

// Problem vectors
data->l = (float*) c_malloc(20 * sizeof(float));
data->l[0] = -qpINFTY;
data->l[1] = -qpINFTY;
data->l[2] = -qpINFTY;
data->l[3] = -qpINFTY;
data->l[4] = -qpINFTY;
data->l[5] = -qpINFTY;
data->l[6] = -qpINFTY;
data->l[7] = -qpINFTY;
data->l[8] = -qpINFTY;
data->l[9] = -qpINFTY;
data->l[10] = 0.00000000000000000000;
data->l[11] = 0.00000000000000000000;
data->l[12] = 0.00000000000000000000;
data->l[13] = 0.00000000000000000000;
data->l[14] = 0.00000000000000000000;
data->l[15] = 0.00000000000000000000;
data->l[16] = 0.00000000000000000000;
data->l[17] = 0.00000000000000000000;
data->l[18] = 0.00000000000000000000;
data->l[19] = 0.00000000000000000000;
data->u = (float*) c_malloc(20 * sizeof(float));
data->u[0] = -1.00000000000000000000;
data->u[1] = -1.00000000000000000000;
data->u[2] = -1.00000000000000000000;
data->u[3] = -1.00000000000000000000;
data->u[4] = -1.00000000000000000000;
data->u[5] = -1.00000000000000000000;
data->u[6] = -1.00000000000000000000;
data->u[7] = -1.00000000000000000000;
data->u[8] = -1.00000000000000000000;
data->u[9] = -1.00000000000000000000;
data->u[10] = qpINFTY;
data->u[11] = qpINFTY;
data->u[12] = qpINFTY;
data->u[13] = qpINFTY;
data->u[14] = qpINFTY;
data->u[15] = qpINFTY;
data->u[16] = qpINFTY;
data->u[17] = qpINFTY;
data->u[18] = qpINFTY;
data->u[19] = qpINFTY;
data->q = (float*) c_malloc(20 * sizeof(float));
data->q[0] = 0.00000000000000000000;
data->q[1] = 0.00000000000000000000;
data->q[2] = 0.00000000000000000000;
data->q[3] = 0.00000000000000000000;
data->q[4] = 0.00000000000000000000;
data->q[5] = 0.00000000000000000000;
data->q[6] = 0.00000000000000000000;
data->q[7] = 0.00000000000000000000;
data->q[8] = 0.00000000000000000000;
data->q[9] = 0.00000000000000000000;
data->q[10] = 1.00000000000000000000;
data->q[11] = 1.00000000000000000000;
data->q[12] = 1.00000000000000000000;
data->q[13] = 1.00000000000000000000;
data->q[14] = 1.00000000000000000000;
data->q[15] = 1.00000000000000000000;
data->q[16] = 1.00000000000000000000;
data->q[17] = 1.00000000000000000000;
data->q[18] = 1.00000000000000000000;
data->q[19] = 1.00000000000000000000;


// Matrix A
//---------
data->A = (csc*) malloc(sizeof(csc));
data->A->m = 20;
data->A->n = 20;
data->A->nz = -1;
data->A->nzmax = 70;
data->A->x = (float*) malloc(70 * sizeof(float));
data->A->x[0] = 0.19158893420206205005;
data->A->x[1] = 0.14494233294725600292;
data->A->x[2] = -0.16300428943430142481;
data->A->x[3] = 0.07666491708488196166;
data->A->x[4] = -0.17870475788235587467;
data->A->x[5] = 0.26595613455476962983;
data->A->x[6] = -0.19783499425640985181;
data->A->x[7] = 0.38993701156555271581;
data->A->x[8] = 0.12319978813583790100;
data->A->x[9] = -0.05751591832126282111;
data->A->x[10] = -0.12421893274074219393;
data->A->x[11] = -0.01235900214958017551;
data->A->x[12] = -0.14885590673501225556;
data->A->x[13] = 0.00470144377080186637;
data->A->x[14] = 0.24108011820263472447;
data->A->x[15] = 0.13778620374984859587;
data->A->x[16] = 0.28480006613549951888;
data->A->x[17] = 0.19453000157924121849;
data->A->x[18] = -0.07434425450439319238;
data->A->x[19] = -0.04780254404420808401;
data->A->x[20] = 0.38567072963187121193;
data->A->x[21] = 0.28754378395301266602;
data->A->x[22] = 0.23045623624367633786;
data->A->x[23] = -0.02010434490071751068;
data->A->x[24] = 0.29733560855125101829;
data->A->x[25] = 0.30891136112406603065;
data->A->x[26] = -0.04882969120344768377;
data->A->x[27] = -0.14178540040022505342;
data->A->x[28] = -0.08298011525809939615;
data->A->x[29] = 0.07974013141059665966;
data->A->x[30] = 0.31239567476385721179;
data->A->x[31] = 0.17024194944461251700;
data->A->x[32] = 0.16247107914434760767;
data->A->x[33] = 0.19439250022629289694;
data->A->x[34] = 0.05163934385638805497;
data->A->x[35] = 0.09630987695107322277;
data->A->x[36] = 0.23021876103919566847;
data->A->x[37] = 0.12644391591552198162;
data->A->x[38] = 0.24838625094456695530;
data->A->x[39] = -0.14399976250133683653;
data->A->x[40] = -0.14125651047491069590;
data->A->x[41] = 0.09715142124472293805;
data->A->x[42] = -0.04940972576560687113;
data->A->x[43] = 0.38790905784807760970;
data->A->x[44] = 0.11080591739816407493;
data->A->x[45] = 0.34770578578896127464;
data->A->x[46] = 0.13041233907345786691;
data->A->x[47] = -0.14436157409923061623;
data->A->x[48] = 0.01493439850455613449;
data->A->x[49] = 0.07783828507048379253;
data->A->x[50] = -1.00000000000000000000;
data->A->x[51] = 1.00000000000000000000;
data->A->x[52] = -1.00000000000000000000;
data->A->x[53] = 1.00000000000000000000;
data->A->x[54] = -1.00000000000000000000;
data->A->x[55] = 1.00000000000000000000;
data->A->x[56] = -1.00000000000000000000;
data->A->x[57] = 1.00000000000000000000;
data->A->x[58] = -1.00000000000000000000;
data->A->x[59] = 1.00000000000000000000;
data->A->x[60] = -1.00000000000000000000;
data->A->x[61] = 1.00000000000000000000;
data->A->x[62] = -1.00000000000000000000;
data->A->x[63] = 1.00000000000000000000;
data->A->x[64] = -1.00000000000000000000;
data->A->x[65] = 1.00000000000000000000;
data->A->x[66] = -1.00000000000000000000;
data->A->x[67] = 1.00000000000000000000;
data->A->x[68] = -1.00000000000000000000;
data->A->x[69] = 1.00000000000000000000;
data->A->i = (int*) malloc(70 * sizeof(int));
data->A->i[0] = 2;
data->A->i[1] = 3;
data->A->i[2] = 9;
data->A->i[3] = 5;
data->A->i[4] = 9;
data->A->i[5] = 3;
data->A->i[6] = 9;
data->A->i[7] = 2;
data->A->i[8] = 4;
data->A->i[9] = 5;
data->A->i[10] = 6;
data->A->i[11] = 7;
data->A->i[12] = 8;
data->A->i[13] = 9;
data->A->i[14] = 1;
data->A->i[15] = 2;
data->A->i[16] = 3;
data->A->i[17] = 4;
data->A->i[18] = 5;
data->A->i[19] = 7;
data->A->i[20] = 1;
data->A->i[21] = 2;
data->A->i[22] = 4;
data->A->i[23] = 6;
data->A->i[24] = 1;
data->A->i[25] = 2;
data->A->i[26] = 5;
data->A->i[27] = 6;
data->A->i[28] = 8;
data->A->i[29] = 9;
data->A->i[30] = 0;
data->A->i[31] = 1;
data->A->i[32] = 3;
data->A->i[33] = 4;
data->A->i[34] = 8;
data->A->i[35] = 9;
data->A->i[36] = 0;
data->A->i[37] = 1;
data->A->i[38] = 2;
data->A->i[39] = 6;
data->A->i[40] = 7;
data->A->i[41] = 8;
data->A->i[42] = 9;
data->A->i[43] = 0;
data->A->i[44] = 1;
data->A->i[45] = 3;
data->A->i[46] = 4;
data->A->i[47] = 5;
data->A->i[48] = 6;
data->A->i[49] = 8;
data->A->i[50] = 0;
data->A->i[51] = 10;
data->A->i[52] = 1;
data->A->i[53] = 11;
data->A->i[54] = 2;
data->A->i[55] = 12;
data->A->i[56] = 3;
data->A->i[57] = 13;
data->A->i[58] = 4;
data->A->i[59] = 14;
data->A->i[60] = 5;
data->A->i[61] = 15;
data->A->i[62] = 6;
data->A->i[63] = 16;
data->A->i[64] = 7;
data->A->i[65] = 17;
data->A->i[66] = 8;
data->A->i[67] = 18;
data->A->i[68] = 9;
data->A->i[69] = 19;
data->A->p = (int*) malloc((20 + 1) * sizeof(int));
data->A->p[0] = 0;
data->A->p[1] = 3;
data->A->p[2] = 5;
data->A->p[3] = 7;
data->A->p[4] = 14;
data->A->p[5] = 20;
data->A->p[6] = 24;
data->A->p[7] = 30;
data->A->p[8] = 36;
data->A->p[9] = 43;
data->A->p[10] = 50;
data->A->p[11] = 52;
data->A->p[12] = 54;
data->A->p[13] = 56;
data->A->p[14] = 58;
data->A->p[15] = 60;
data->A->p[16] = 62;
data->A->p[17] = 64;
data->A->p[18] = 66;
data->A->p[19] = 68;
data->A->p[20] = 70;


// Matrix P
//---------
data->P = (csc*) malloc(sizeof(csc));
data->P->m = 20;
data->P->n = 20;
data->P->nz = -1;
data->P->nzmax = 10;
data->P->x = (float*) malloc(10 * sizeof(float));
data->P->x[0] = 1.00000000000000000000;
data->P->x[1] = 1.00000000000000000000;
data->P->x[2] = 1.00000000000000000000;
data->P->x[3] = 1.00000000000000000000;
data->P->x[4] = 1.00000000000000000000;
data->P->x[5] = 1.00000000000000000000;
data->P->x[6] = 1.00000000000000000000;
data->P->x[7] = 1.00000000000000000000;
data->P->x[8] = 1.00000000000000000000;
data->P->x[9] = 1.00000000000000000000;
data->P->i = (int*) malloc(10 * sizeof(int));
data->P->i[0] = 0;
data->P->i[1] = 1;
data->P->i[2] = 2;
data->P->i[3] = 3;
data->P->i[4] = 4;
data->P->i[5] = 5;
data->P->i[6] = 6;
data->P->i[7] = 7;
data->P->i[8] = 8;
data->P->i[9] = 9;
data->P->p = (int*) malloc((20 + 1) * sizeof(int));
data->P->p[0] = 0;
data->P->p[1] = 1;
data->P->p[2] = 2;
data->P->p[3] = 3;
data->P->p[4] = 4;
data->P->p[5] = 5;
data->P->p[6] = 6;
data->P->p[7] = 7;
data->P->p[8] = 8;
data->P->p[9] = 9;
data->P->p[10] = 10;
data->P->p[11] = 10;
data->P->p[12] = 10;
data->P->p[13] = 10;
data->P->p[14] = 10;
data->P->p[15] = 10;
data->P->p[16] = 10;
data->P->p[17] = 10;
data->P->p[18] = 10;
data->P->p[19] = 10;
data->P->p[20] = 10;

return data;

}

/* function to clean problem data structure */
void clean_problem_svm_example(qpData * data){

// Clean vectors
c_free(data->l);
c_free(data->u);
c_free(data->q);

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

