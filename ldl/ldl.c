#include <stdbool.h>
#include "../include/ldl.h"
#include "limits.h"

#define LDL_UNKNOWN (-1)
#define LDL_USED (1)
#define LDL_UNUSED (0)

/* compute the elimination tree for a in compressed sparse column form */

int LDL_etree(const int n, const int* Ap, const int* Ai, int* work, int* Lnz, int* etree)
{

    int sumLnz;
    int ii,jj,pp;


    for(ii = 0; ii < n; ii++)
    {
        /* zero out Lnz and work.  Set all etree values to unknown */
        work[ii]  = 0;
        Lnz[ii]   = 0;
        etree[ii] = LDL_UNKNOWN;

        /* quit if A doesn't have at least one entry in every column */
        if(Ap[ii] == Ap[ii+1])
        {
            return -1;
        }

    }

    for(jj = 0; jj < n; jj++)
    {
        work[jj] = jj;
        for(pp = Ap[jj]; pp < Ap[jj+1]; pp++)
        {
            ii = Ai[pp];
            if(ii > jj)
            {
                return -1;
            }; // quit if entries on lower triangle
            while(work[ii] != jj)
            {
                if(etree[ii] == LDL_UNKNOWN)
                {
                    etree[ii] = jj;
                }
                Lnz[ii]++; //nonzeros in this column
                work[ii] = jj;
                ii = etree[ii];
            }
        }
    }

    /* compute the total nonzeros in L */
    sumLnz  = 0;
    for(ii = 0; ii < n; ii++)
    {
        if(sumLnz > INT_MAX - Lnz[ii])
        {
            sumLnz = -2; // return -2 if the nz count overflows uint type
            break;
        }
        else
        {
            sumLnz += Lnz[ii];
        }
    }

    return sumLnz;
}

int LDL_factor(const int n, const int* Ap, const int* Ai, const float* Ax, int* Lp, int* Li, float* Lx, float* D, float* Dinv, const int* Lnz, const int* etree, bool*  bwork, int* iwork, float* fwork)
{

    int ii,jj,kk,nnzY, bidx, cidx, nextIdx, nnzE, tmpIdx;
    int *yIdx, *elimBuffer, *LNextSpaceInCol;
    float *yVals;
    float yVals_cidx;
    bool *yMarkers;
    int positiveValuesInD = 0;

    /* partition working memory into pieces */
    yMarkers = bwork;
    yIdx = iwork;
    elimBuffer = iwork + n;
    LNextSpaceInCol = iwork + n*2;
    yVals = fwork;

    Lp[0] = 0;

    for(ii = 0; ii < n; ii++)
    {

        /* compute L column indices */
        Lp[ii+1] = Lp[ii] + Lnz[ii];

        /* set all Yidx to be 'unused' initially */
        yMarkers[ii] = LDL_UNUSED;
        yVals[ii] = 0.0;
        D[ii] = 0.0;
        LNextSpaceInCol[ii] = Lp[ii];
    }

    /* First element of the diagonal D */
    D[0] = Ax[0];
    if(D[0] == 0.0)
    {
        return -1;
    }
    if(D[0]  > 0.0)
    {
        positiveValuesInD++;
    }
    Dinv[0] = 1/D[0];

    /* start from 1 here. The upper LH corner is trivially 0 */
    for(kk = 1; kk < n; kk++)
    {

        /*
        for each k, we compute a solution to
        y = L(0:(k-1),0:k-1))\b, where b is the kth
        column of A that sits above the diagonal.
        the solution y is then the kth row of L,
        with an implied '1' at the diagonal entry.
        */

        nnzY = 0; // number of nonzeros in this row of L

        /* this loop determines where nonzeros will go in the kth row of L */
        tmpIdx = Ap[kk+1];

        for(ii = Ap[kk]; ii < tmpIdx; ii++)
        {

            bidx = Ai[ii];

            /*
            initialize D[k] as the element of this column
            corresponding to the diagonal place.  Don't use
            this element as part of the elimination step
            that computes the k^th row of L
            */

            if(bidx == kk)
            {
                D[kk] = Ax[ii];
                continue;
            }

            yVals[bidx] = Ax[ii];

            /*
            use the forward elimination tree to figure
            out which elements must be eliminated after
            this element of b
            */
            nextIdx = bidx;

            if(yMarkers[nextIdx] == LDL_UNUSED)
            {   //this y term not already visited

                yMarkers[nextIdx] = LDL_USED;     // mark as visited
                elimBuffer[0] = nextIdx;  // it goes at the start of the current list
                nnzE = 1;         //length of unvisited elimination path from here

                nextIdx = etree[bidx];

                while(nextIdx != LDL_UNKNOWN && nextIdx < kk)
                {
                    if(yMarkers[nextIdx] == LDL_USED)
                        break;

                    yMarkers[nextIdx] = LDL_USED;
                    elimBuffer[nnzE] = nextIdx;
                    nnzE++;
                    nextIdx = etree[nextIdx];
                }


                while(nnzE)
                {
                    yIdx[nnzY++] = elimBuffer[--nnzE];
                }
            }

        }

        /* this for loop places nonzeros values in the kkth row */
        for(ii = (nnzY-1); ii >=0; ii--)
        {

            // current col
            cidx = yIdx[ii];

            tmpIdx = LNextSpaceInCol[cidx];
            yVals_cidx = yVals[cidx];
            for(jj = Lp[cidx]; jj < tmpIdx; jj++)
            {
                yVals[Li[jj]] -= Lx[jj]*yVals_cidx;
            }

           // the cidx^th element of y = L\b.
           //so compute the corresponding element of
           //this row of L and put it into the right place
            Li[tmpIdx] = kk;
            Lx[tmpIdx] = yVals_cidx *Dinv[cidx];

            D[kk] -= yVals_cidx*Lx[tmpIdx];
            LNextSpaceInCol[cidx]++;

            //reset the yvalues and indices back to zero and LDL_UNUSEDs
            yVals[cidx] = 0.0;
            yMarkers[cidx] = LDL_UNUSED;

        }

        /*
        maintain a count of the positive entries
        in D.  If we hit a zero, we can't factor
        this matrix, so abort
        */
        if(D[kk] == 0.0)
        {
            return -1;
        }
        if(D[kk]  > 0.0)
        {
            positiveValuesInD++;
        }

        //compute the inverse of the diagonal
        Dinv[kk]= 1/D[kk];

    } //end for k

    return positiveValuesInD;

}

// Solves (L+I)x = b
void LDL_Lsolve(const int n, const int* Lp, const int* Li, const float* Lx, float* x)
{

    int ii,jj;

    for(ii = 0; ii < n; ii++)
    {
        for(jj = Lp[ii]; jj < Lp[ii+1]; jj++)
        {
            x[Li[jj]] -= Lx[jj]*x[ii];
        }
    }

}

// Solves (L+I)'x = b
void LDL_Ltsolve(const int n, const int* Lp, const int* Li, const float* Lx, float* x)
{

    int ii,jj;

    for(ii = n-1; ii>=0; ii--)
    {
        for(jj = Lp[ii]; jj < Lp[ii+1]; jj++)
        {
            x[ii] -= Lx[jj]*x[Li[jj]];
        }
    }
}

// Solves Ax = b where A has given LDL factors
void LDL_solve(const int n, const int* Lp, const int* Li, const float* Lx, const float* Dinv, float* x)
{
    int ii;

    LDL_Lsolve(n,Lp,Li,Lx,x);

    for(ii = 0; ii < n; ii++)
        x[ii] *= Dinv[ii];

    LDL_Ltsolve(n,Lp,Li,Lx,x);

}
