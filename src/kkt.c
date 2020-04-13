#include "../include/kkt.h"

/* todo: figure out triplet format */

/* form kkt matrix */
csc* form_KKT(const csc *P, const csc *A, c_int format, c_float param1, c_float *param2, c_int *PtoKKT, c_int *AtoKKT, c_int **Pdiag_idx, c_int *Pdiag_n, c_int *param2toKKT)
{
  c_int  nKKT, nnzKKTmax; // size, number of nonzeros and max number of nonzeros
  csc   *KKT_trip, *KKT;  // KKT matrix in triplet format and CSC format
  c_int  ptr, ii, jj;     // counters for elements (ii,jj) and index pointer
  c_int  zKKT = 0;        // counter for total number of elements in P
  c_int *KKT_TtoC;        // Pointer to vector mapping from KKT in triplet form to CSC

  // get matrix dimensions
  nKKT = P->m + A->m;

  // get maximum number of nonzero elements (only upper triangular part)
  nnzKKTmax = P->p[P->n] + // number of elements in P
              P->m +       // number of elements in param1 * I
              A->p[A->n] + // number of nonzeros in A
              A->m;        // number of elements in - diag(param2)

  // preallocate KKT matrix in triplet format
  KKT_trip = csc_spalloc(nKKT, nKKT, nnzKKTmax, 1, 1);

  if (!KKT_trip)
	  return NULL;  // failed to preallocate matrix

  // allocate vector of indices on the diagonal. worst case it has m elements
  if (Pdiag_idx != NULL)
  {
    (*Pdiag_idx) = c_malloc(P->m * sizeof(c_int));
    *Pdiag_n     = 0; // set 0 diagonal elements to start
  }

  // allocate Triplet matrices
  // P + param1 I
  for (jj = 0; jj < P->n; jj++)
  { // cycle over columns
    // no elements in column j => add diagonal element param1
    if (P->p[jj] == P->p[jj + 1])
    {
      KKT_trip->ii[zKKT] = jj;
      KKT_trip->p[zKKT] = jj;
      KKT_trip->x[zKKT] = param1;
      zKKT++;
    }

    for (ptr = P->p[jj]; ptr < P->p[jj + 1]; ptr++)
    { // cycle over rows
      // get current row
      ii = P->ii[ptr];

      // add element of P
      KKT_trip->ii[zKKT] = ii;
      KKT_trip->p[zKKT] = jj;
      KKT_trip->x[zKKT] = P->x[ptr];

      if (PtoKKT != NULL)
	      PtoKKT[ptr] = zKKT;  // update index from P to KTTtrip

      if (ii == jj)
      {                                 // P has a diagonal element,
                                                    // add param1
        KKT_trip->x[zKKT] += param1;

        // If index vector pointer supplied -> Store the index
        if (Pdiag_idx != OSQP_NULL)
       	{
          (*Pdiag_idx)[*Pdiag_n] = ptr;
          (*Pdiag_n)++;
        }
      }
      zKKT++;

      // Add diagonal param1 in case
      if ((ii < jj) && (ptr + 1 == P->[j+1]))	      // Diagonal element not reached, last element of col j
      {
        // Add diagonal element param1
        KKT_trip->ii[zKKT] = jj;
        KKT_trip->p[zKKT] = jj;
        KKT_trip->x[zKKT] = param1;
        zKKT++;
      }
    }
  }

  if (Pdiag_idx != NULL)
  {
    // Realloc Pdiag_idx so that it contains exactly *Pdiag_n diagonal elements
    (*Pdiag_idx) = c_realloc((*Pdiag_idx), (*Pdiag_n) * sizeof(c_int));
  }


  // A' at top right
  for (jj = 0; jj < A->n; jj++)
  {                      // Cycle over columns of A
    for (ptr = A->p[jj]; ptr < A->p[jj + 1]; ptr++)
    {
      KKT_trip->p[zKKT] = P->m + A->i[ptr];         // Assign column index from
                                                    // row index of A
      KKT_trip->ii[zKKT] = jj;                        // Assign row index from
                                                    // column index of A
      KKT_trip->x[zKKT] = A->x[ptr];                // Assign A value element

      if (AtoKKT != NULL)
	      AtoKKT[ptr] = zKKT;  // Update index from A to
                                                    // KKTtrip
      zKKT++;
    }
  }

  // - diag(param2) at bottom right
  for (jj = 0; jj < A->m; jj++)
  {
    KKT_trip->ii[zKKT] = jj + P->n;
    KKT_trip->p[zKKT] = jj + P->n;
    KKT_trip->x[zKKT] = -param2[jj];

    if (param2toKKT != NULL) 
	    param2toKKT[j] = zKKT;  // Update index from
                                                          // param2 to KKTtrip
    zKKT++;
  }

  // Allocate number of nonzeros
  KKT_trip->nz = zKKT;

  // Convert triplet matrix to csc format
  if (!PtoKKT && !AtoKKT && !param2toKKT)
  {
    // If no index vectors passed, do not store KKT mapping from Trip to CSC/CSR
    if (format == 0)
	    KKT = triplet_to_csc(KKT_trip, NULL);
    else 
	    KKT = triplet_to_csr(KKT_trip, NULL);
  }
  else
  {
    // Allocate vector of indices from triplet to csc
    KKT_TtoC = c_malloc((zKKT) * sizeof(c_int));

    if (!KKT_TtoC)
    {
      // Error in allocating KKT_TtoC vector
      csc_spfree(KKT_trip);
      c_free(*Pdiag_idx);
      return NULL;
    }

    // Store KKT mapping from Trip to CSC/CSR
    if (format == 0)
      KKT = triplet_to_csc(KKT_trip, KKT_TtoC);
    else
      KKT = triplet_to_csr(KKT_trip, KKT_TtoC);

    // Update vectors of indices from P, A, param2 to KKT (now in CSC format)
    if (PtoKKT != OSQP_NULL)
    {
      for (ii = 0; ii < P->p[P->n]; ii++)
      {
        PtoKKT[ii] = KKT_TtoC[PtoKKT[ii]];
      }
    }

    if (AtoKKT != NULL)
    {
      for (ii = 0; ii < A->p[A->n]; ii++)
      {
        AtoKKT[ii] = KKT_TtoC[AtoKKT[ii]];
      }
    }

    if (param2toKKT != NULL)
    {
      for (ii = 0; ii < A->m; ii++)
      {
        param2toKKT[ii] = KKT_TtoC[param2toKKT[ii]];
      }
    }

    // Free mapping
    c_free(KKT_TtoC);
  }

  // Clean matrix in triplet format and return result
  csc_spfree(KKT_trip);

  return KKT;
}


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


