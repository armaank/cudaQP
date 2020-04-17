/* basic tests of ldl factorization routines */

#include <stdio.h>
#include <stdlib.h>
#include "../include/ldl.h"
#include "../include/ldl_unittest.h"
#define ERR 1e-4

float vec_diff_norm(float* x, float* y, int len);
int ldl_factor_solve(int An, int* Ap, int* Ai, float* Ax, float* b);
static char* basic();

int tests_run = 0;

static char* run_tests()
{
    ut_run_test(basic);

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
        printf("test passed\n");
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



