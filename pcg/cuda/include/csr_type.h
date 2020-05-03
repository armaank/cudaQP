
#ifndef CSR_TYPE_H
# define CSR_TYPE_H


#include <cusparse_v2.h>
// #include "osqp_api_types.h"   /* --> c_int, c_float */


/* CSR matrix structure */
struct csr_t {
  int               m;          ///< number of rows
  int               n;          ///< number of columns
  int              *row_ptr;    ///< row pointers (size m+1)
  int              *row_ind;    ///< uncompressed row indices (size nnz), NULL if not needed 
  int              *col_ind;    ///< column indices (size nnz)
  float            *val;        ///< numerical values (size nnz)
  int               nnz;        ///< number of non-zero entries in matrix

  void               *buffer;
  size_t              bufferSizeInBytes;
  cusparseAlgMode_t   alg;
  cusparseMatDescr_t  MatDescription;
};


#endif /* ifndef CSR_TYPE_H */
