#include "../include/kkt.h"

/* todo: figure out triplet format, check nulls */

/* form kkt matrix */
csc* form_KKT(const csc *P, const csc *A, int format, float param1, float *param2, int *PtoKKT, int *AtoKKT, int **Pdiag_idx, int *Pdiag_n, int *param2toKKT)
{
    /* size, number of nonzeros and max number of nonzeros */
    int nKKT, nnzKKTmax;
    /* ktt matrix in triplet and csc format */
    csc *KKT_trip, *KKT;
    /* counts and index pointer */
    int ptr, ii, jj;
    /* counter for num of elements in P */
    int zKKT = 0;
    /* pointer to vector from KKT in triplet to csc */
    int *KKT_TtoC;

    /* get matrix dimensions */
    nKKT = P->m + A->m;

    /* get maximum number of nonzero elements in the UT part of the matrix */
    nnzKKTmax = P->p[P->n] + // number of elements in P
                P->m +       // number of elements in param1 * I
                A->p[A->n] + // number of nonzeros in A
                A->m;        // number of elements in - diag(param2)

    /* preallocate KKT matrix in triplet format */
    KKT_trip = csc_spalloc(nKKT, nKKT, nnzKKTmax, 1, 1);

    if (!KKT_trip)
        return 0;  // failed to preallocate matrix

    /* allocate vector of indices on the diagonal. worst case it has m elements */
    if (Pdiag_idx != 0)
    {
        (*Pdiag_idx) = malloc(P->m * sizeof(int));
        *Pdiag_n = 0; // set 0 diagonal elements to start
    }

    /* allocate Triplet matrices */
    for (jj = 0; jj < P->n; jj++)
    {   // cycle over columns
        // no elements in column jj => add diagonal element param1
        if (P->p[jj] == P->p[jj + 1])
        {
            KKT_trip->i[zKKT] = jj;
            KKT_trip->p[zKKT] = jj;
            KKT_trip->x[zKKT] = param1;
            zKKT++;
        }

        for (ptr = P->p[jj]; ptr < P->p[jj + 1]; ptr++)
        {   // cycle over rows
            // get current row
            ii = P->i[ptr];

            // add element of P
            KKT_trip->i[zKKT] = ii;
            KKT_trip->p[zKKT] = jj;
            KKT_trip->x[zKKT] = P->x[ptr];

            if (PtoKKT != 0)
                PtoKKT[ptr] = zKKT;  // update index from P to KTTtrip

            if (ii == jj)
            {   // P has a diagonal element,
                // add param1
                KKT_trip->x[zKKT] += param1;

                // If index vector pointer supplied -> Store the index
                if (Pdiag_idx != 0)
                {
                    (*Pdiag_idx)[*Pdiag_n] = ptr;
                    (*Pdiag_n)++;
                }
            }
            zKKT++;

            // Add diagonal param1 in case
            if ((ii < jj) && (ptr + 1 == P->p[jj+1]))
            {   // Diagonal element not reached, last element of col jj
                // Add diagonal element param1
                KKT_trip->i[zKKT] = jj;
                KKT_trip->p[zKKT] = jj;
                KKT_trip->x[zKKT] = param1;
                zKKT++;
            }
        }
    }

    if (Pdiag_idx != 0)
    {
        /* realloc Pdiag_idx so that it contains exactly *Pdiag_n diagonal elements */
        (*Pdiag_idx) = realloc((*Pdiag_idx), (*Pdiag_n) * sizeof(int));
    }


    /* A' at top right */
    for (jj = 0; jj < A->n; jj++)
    {   // cycle over columns of A
        for (ptr = A->p[jj]; ptr < A->p[jj + 1]; ptr++)
        {
            KKT_trip->p[zKKT] = P->m + A->i[ptr]; // assign column index from
            // row index of A
            KKT_trip->i[zKKT] = jj; // assign row index from
            // column index of A
            KKT_trip->x[zKKT] = A->x[ptr]; // assign A value element

            if (AtoKKT != 0)
                AtoKKT[ptr] = zKKT; // Update index from A to
            // KKTtrip
            zKKT++;
        }
    }

    /* - diag(param2) at bottom right */
    for (jj = 0; jj < A->m; jj++)
    {
        KKT_trip->i[zKKT] = jj + P->n;
        KKT_trip->p[zKKT] = jj + P->n;
        KKT_trip->x[zKKT] = -param2[jj];

        if (param2toKKT != 0)
            param2toKKT[jj] = zKKT;  // update index from
        // param2 to KKTtrip
        zKKT++;
    }

    /* allocate number of nonzeros */
    KKT_trip->nz = zKKT;

    /* convert triplet matrix to csc format */
    if (!PtoKKT && !AtoKKT && !param2toKKT)
    {
        // if no index vectors passed, do not store KKT mapping from Trip to CSC/CSR
        if (format == 0)
            KKT = triplet_to_csc(KKT_trip, 0);
        else
            KKT = triplet_to_csr(KKT_trip, 0);
    }
    else
    {
        // allocate vector of indices from triplet to csc
        KKT_TtoC = malloc((zKKT) * sizeof(int));

        if (!KKT_TtoC)
        {
            // error in allocating KKT_TtoC vector
            csc_spfree(KKT_trip);
            free(*Pdiag_idx);
            return 0;
        }

        // store KKT mapping from Trip to CSC/CSR
        if (format == 0)
            KKT = triplet_to_csc(KKT_trip, KKT_TtoC);
        else
            KKT = triplet_to_csr(KKT_trip, KKT_TtoC);

        // update vectors of indices from P, A, param2 to KKT (now in CSC format)
        if (PtoKKT != 0)
        {
            for (ii = 0; ii < P->p[P->n]; ii++)
            {
                PtoKKT[ii] = KKT_TtoC[PtoKKT[ii]];
            }
        }

        if (AtoKKT != 0)
        {
            for (ii = 0; ii < A->p[A->n]; ii++)
            {
                AtoKKT[ii] = KKT_TtoC[AtoKKT[ii]];
            }
        }

        if (param2toKKT != 0)
        {
            for (ii = 0; ii < A->m; ii++)
            {
                param2toKKT[ii] = KKT_TtoC[param2toKKT[ii]];
            }
        }

        // Free mapping
        free(KKT_TtoC);
    }

    // free matrix in triplet format and return result
    csc_spfree(KKT_trip);

    return KKT;
}

void update_KKT_P(csc *KKT, const csc *P, const int *PtoKKT, const float param1, const int *Pdiag_idx, const int Pdiag_n)
{
    int ii, jj;

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

void update_KKT_A(csc *KKT, const csc *A, const int *AtoKKT)
{
    int ii;

    // update elements of KKT using A
    for (ii = 0; ii < A->p[A->n]; ii++)
    {
        KKT->x[AtoKKT[ii]] = A->x[ii];
    }
}

void update_KKT_param2(csc *KKT, const float *param2, const int *param2toKKT, const int m)
{
    int ii;

    // update elements of KKT using param2
    for (ii = 0; ii < m; ii++) {
        KKT->x[param2toKKT[ii]] = -param2[ii];
    }
}
