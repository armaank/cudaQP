/* defining matrix/vector types for cuda interface */

#ifndef ALGEBRA_TYPES_H
# define ALGEBRA_TYPES_H

# ifdef __cplusplus
extern "C" {
# endif


#include "osqp_api_types.h"

struct OSQPVectori_ {
    c_int *d_val;
    c_int  length;
};

struct OSQPVectorf_ {
    c_float *d_val;
    c_int    length;
};


/* matrix in CSR format stored in GPU memory */
typedef struct csr_t csr;

struct OSQPMatrix_ {
    csr     *S;   /* P or A */
    csr     *At;
    c_int   *d_A_to_At_ind;
    c_float *d_P_triu_val;
    c_int   *d_P_triu_to_full_ind;
    c_int   *d_P_diag_ind;
    c_int    P_triu_nnz;
    c_int    symmetric;
};


# ifdef __cplusplus
}
# endif

#endif /* ifndef ALGEBRA_TYPES_H */
