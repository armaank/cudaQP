/* basic tests of ldl factorization routines */

#include <stdio.h>
#include <stdlib.h>
#include "../ldl.h"
#include "../../include/unittest.h"
#define ERR 1e-4

float vec_diff_norm(float* x, float* y, int len);
int ldl_factor_solve(int An, int* Ap, int* Ai, float* Ax, float* b);
static char* basic();
static char* kkt();
static char* rank();

int tests_run = 0;

static char* run_tests()
{
    ut_run_test(basic);
    ut_run_test(kkt);
    ut_run_test(rank);

    return 0;
}

int main(void) {
    char *result = run_tests();

    if (result != 0)
    {
        printf("%s\n", result);
    }
    else
    {
        printf("tests passed\n");
    }

    printf("tests run: %d\n", tests_run);

    return result != 0;
}

float vec_diff_norm(float* x, float* y, int len) {

    float maxDiff = 0.0;
    float elDiff  = 0.0;
    int ii;

    for(ii = 0; ii < len; ii++)
    {
        elDiff  = x[ii] - y[ii];
        maxDiff = (elDiff > maxDiff) ? elDiff : ((-elDiff > maxDiff) ? -elDiff : maxDiff);
    }
    return maxDiff;

}

int ldl_factor_solve(int An, int* Ap, int* Ai, float* Ax, float* b)
{

    /* data for L and D factors */
    int Ln = An;
    int *Lp;
    int *Li;
    float *Lx;
    float *D;
    float *Dinv;

    /* data for elim tree */
    int *etree;
    int *Lnz;
    int  sumLnz;

    /* data for factorisation */
    int   *iwork;
    bool  *bwork;
    float *fwork;
    int   factorStatus;

    /*--------------------------------
     * pre-factorisation memory allocations
     *---------------------------------*/

    //These can happen *before* the etree is calculated
    //since the sizes are not sparsity pattern specific

    //For the elimination tree
    etree = (int*)malloc(sizeof(int)*An);
    Lnz   = (int*)malloc(sizeof(int)*An);

    //For the L factors.   Li and Lx are sparsity dependent
    //so must be done after the etree is constructed
    Lp    = (int*)malloc(sizeof(int)*(An+1));
    D     = (float*)malloc(sizeof(float)*An);
    Dinv  = (float*)malloc(sizeof(float)*An);

    //Working memory.  Note that both the etree and factor
    //calls requires a working vector of QDLDL_int, with
    //the factor function requiring 3*An elements and the
    //etree only An elements.   Just allocate the larger
    //amount here and use it in both places
    iwork = (int*)malloc(sizeof(int)*(3*An));
    bwork = (bool*)malloc(sizeof(bool)*An);
    fwork = (float*)malloc(sizeof(float)*An);


    /*--------------------------------
     * elimination tree calculation
     *---------------------------------*/
    sumLnz = LDL_etree(An,Ap,Ai,iwork,Lnz,etree);

    //not perfect triu A = bomb
    if(sumLnz < 0) {
        free(Lp);
        free(D);
        free(Dinv);
        free(etree);
        free(Lnz);
        free(iwork);
        free(bwork);
        free(fwork);
        return sumLnz;
    }

    /*--------------------------------
     * LDL factorisation
     *---------------------------------*/

    //First allocate memory for Li and Lx
    Li    = (int*)malloc(sizeof(int)*sumLnz);
    Lx    = (float*)malloc(sizeof(float)*sumLnz);

    //now factor
    factorStatus = LDL_factor(An,Ap,Ai,Ax,Lp,Li,Lx,D,Dinv,Lnz,etree,bwork,iwork,fwork);

    //Zero on the diagonal = bomb
    if(factorStatus < 0) {
        free(Lp);
        free(D);
        free(Dinv);
        free(etree);
        free(Lnz);
        free(iwork);
        free(bwork);
        free(fwork);
        free(Li);
        free(Lx);
        return factorStatus;
    }

    /*--------------------------------
     * solve
     *---------------------------------*/
    LDL_solve(Ln,Lp,Li,Lx,Dinv,b);


    /*--------------------------------
     * clean up
     *---------------------------------*/
    free(Lp);
    free(D);
    free(Dinv);
    free(etree);
    free(Lnz);
    free(iwork);
    free(bwork);
    free(fwork);
    free(Li);
    free(Lx);

    return 0 ;

}


static char* basic()
{
    //A matrix data
    int Ap[]  = {0, 1, 2, 4, 5, 6, 8, 10, 12, 14, 17};
    int Ai[]  = {0, 1, 1, 2, 3, 4, 1, 5, 0, 6, 3, 7, 6, 8, 1, 2, 9};
    float Ax[] = {1.0, 0.460641, -0.121189, 0.417928, 0.177828,
                  0.1, -0.0290058, -1.0, 0.350321, -0.441092, -0.0845395,
                  -0.316228, 0.178663, -0.299077, 0.182452, -1.56506, -0.1
                 };
    int An = 10;

    // RHS and solution to Ax = b
    float b[]    = {1,2,3,4,5,6,7,8,9,10};
    float sol[] = {10.2171, 3.9416, -5.69096, 9.28661, 50.0, -6.11433,
                   -26.3104, -27.7809, -45.8099, -3.74178
                  };

    //x replaces b during solve
    int status = ldl_factor_solve(An,Ap,Ai,Ax,b);

    ut_assert("Factorisation failed", status >= 0);
    ut_assert("Solve accuracy failed", vec_diff_norm(b,sol,An) < ERR);

    return 0;
}

static char* kkt()
{
    // Unordered A
    int Ap[]  = {0, 1, 2, 5, 6, 7, 8, 12};
    int Ai[]  = {0, 1, 2, 1, 0, 3, 4, 5, 5, 6, 4, 3};
    float Ax[] = {-0.25,  -0.25,   1.0,   0.513578,   0.529142,  -0.25,  -0.25,   1.10274,   0.15538,   1.25883,   0.13458,   0.621134};

    int An = 7;

    // RHS and solution to Ax = b
    float b[]    = {-0.595598, -0.0193715, -0.576156, -0.168746, 0.61543, 0.419073, 1.31087};
    float xsol[] = {1.13141, -1.1367, -0.591044, 1.68867, -2.24209, 0.32254, 0.407998};

    //x replaces b during solve
    int status = ldl_factor_solve(An,Ap,Ai,Ax,b);

    ut_assert("Factorisation failed", status >= 0);
    ut_assert("Solve accuracy failed", vec_diff_norm(b,xsol,An) < ERR);

    return 0;
}

static char* rank()
{
    //A matrix data
    int Ap[]  = {0, 1, 3};
    int Ai[]  = {0, 0, 1};
    float Ax[] = {1.0, 1.0, 1.0};
    int An = 2;

    // RHS for Ax = b (should fail to solve)
    float b[]    = {1,1};

    //x replaces b during solve
    int status = ldl_factor_solve(An,Ap,Ai,Ax,b);

    ut_assert("rank deficiency undetected", status < 0);
    return 0;
}
