#include "../include/kkt.h"



void update_KKT_P(csc *KKT, const csc *P, const c_int *PtoKKT, const c_float param1, const c_int *Pdiag_idx, const c_int Pdiag_n)
{
    c_int ii, jj;
    s

    // update elements of KKT using P
    for (ii = 0; ii < P->p[P->n]; ii++)
    {
        KKT->x[PtoKKT[ii]] = P->x[ii];
    }

    // update diagonal elements of KKT by adding sigma
    for (ii = 0; ii < Pdiag_n; ii++)
    {
        jj = Pdiag_idx[ii]; // extract index of the element on the diagonal
        KKT->x[PtoKKT[jj]] += param1;
    }
}

void update_KKT_A(csc *KKT, const csc *A, const c_int *AtoKKT)
{
    c_int ii;

    // update elements of KKT using A
    for (ii = 0; ii < A->p[A->n]; ii++)
    {
        KKT->x[AtoKKT[ii]] = A->x[ii];
    }
}

void update_KKT_param2(csc *KKT, const c_float *param2, const c_int *param2toKKT, const c_int m)
{
    c_int ii;

    // update elements of KKT using param2
    for (ii = 0; ii < m; ii++) {
        KKT->x[param2toKKT[ii]] = -param2[ii];
    }
}


