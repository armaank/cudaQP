#include "../../include/ops.h"

float* csc_to_dns(csc *M)
{
    int i, j = 0; // Predefine row index and column index
    int idx;

    // Initialize matrix of zeros
    float *A = (float *)calloc(M->m * M->n, sizeof(float));

    // Allocate elements
    for (idx = 0; idx < M->p[M->n]; idx++)
    {
        // Get row index i (starting from 1)
        i = M->i[idx];

        // Get column index j (increase if necessary) (starting from 1)
        while (M->p[j + 1] <= idx) {
            j++;
        }

        // Assign values to A
        A[j * (M->m) + i] = M->x[idx];
    }
    return A;
}

int is_eq_csc(csc *A, csc *B, float tol) {
    int j, i;

    // If number of columns does not coincide, they are not equal.
    if (A->n != B->n) return 0;

    for (j = 0; j < A->n; j++) { // Cycle over columns j
        // if column pointer does not coincide, they are not equal
        if (A->p[j] != B->p[j]) return 0;

        for (i = A->p[j]; i < A->p[j + 1]; i++) { // Cycle rows i in column j
            if ((A->i[i] != B->i[i]) ||             // Different row indices
                    (c_absval(A->x[i] - B->x[i]) > tol)) {
                return 0;
            }
        }
    }
    return 1;
}
