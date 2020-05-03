/**
 *  Copyright (c) 2019 ETH Zurich, Automatic Control Lab, Michel Schubiger, Goran Banjac.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

/********************************************************************
 *       Wrapper functions to abstract floating point type          *
 *                                                                  *
 *  They make the code work when either single or double precision  *
 *  floating-point type is used.                                    *
 ********************************************************************/
// note - might not need this
#ifndef CUDA_WRAPPER_H
# define CUDA_WRAPPER_H

#include <cusparse.h>
#include <cublas_v2.h>

// #include "osqp_api_types.h"


static cublasStatus_t cublasTaxpy(cublasHandle_t  handle,
                                  int           n,
                                  const float  *alpha,
                                  const float  *x,
                                  int           incx,
                                  float        *y,
                                  int           incy) {

// #ifdef DFLOAT
//   return cublasSaxpy(handle, n, alpha, x, incx, y, incy);
// #else
  return cublasDaxpy(handle, n, alpha, x, incx, y, incy);
// #endif
}


static cublasStatus_t cublasTscal(cublasHandle_t  handle,
                                  int           n,
                                  const float  *alpha,
                                  float        *x,
                                  int           incx) {

// #ifdef DFLOAT
//   return cublasSscal(handle, n, alpha, x, incx);
// #else
  return cublasDscal(handle, n, alpha, x, incx);
// #endif
}


static cublasStatus_t cublasTdot(cublasHandle_t  handle,
                                 int           n,
                                 const float  *x,
                                 int           incx,
                                 const float  *y,
                                 int           incy,
                                 float        *result) {

// #ifdef DFLOAT
//   return cublasSdot (handle, n, x, incx, y, incy, result);
// #else
  return cublasDdot (handle, n, x, incx, y, incy, result);
// #endif
}


static cublasStatus_t cublasITamax(cublasHandle_t  handle,
                                   int           n,
                                   const float  *x,
                                   int           incx,
                                   int          *result) {

// #ifdef DFLOAT
//   return cublasIsamax(handle, n, x, incx, result);
// #else
  return cublasIdamax(handle, n, x, incx, result);
// #endif
}


static cublasStatus_t cublasTasum(cublasHandle_t  handle,
                                  int           n,
                                  const float  *x,
                                  int           incx,
                                  float        *result) {

// #ifdef DFLOAT
//   return cublasSasum(handle, n, x, incx, result);
// #else
  return cublasDasum(handle, n, x, incx, result);
// #endif
}


static cusparseStatus_t cusparseTgthr(cusparseHandle_t     handle,
                                      int                nnz,
                                      const float       *y,
                                      float             *xVal,
                                      const int         *xInd,
                                      cusparseIndexBase_t  idxBase) {

// #ifdef DFLOAT
//   return cusparseSgthr(handle, nnz, y, xVal, xInd, idxBase);
// #else
  return cusparseDgthr(handle, nnz, y, xVal, xInd, idxBase);
// #endif
}


static cublasStatus_t cublasTnrm2(cublasHandle_t  handle,
                                  int           n,
                                  const float  *x,
                                  int           incx,
                                  float        *result) {

// #ifdef DFLOAT
//   return cublasSnrm2(handle, n, x, incx, result);
// #else
  return cublasDnrm2(handle, n, x, incx, result);
// #endif
}


static cusparseStatus_t cusparseCsrmv(cusparseHandle_t          handle,
                                      cusparseAlgMode_t         alg,
                                      int                     m,
                                      int                     n,
                                      int                     nnz,
                                      const float            *alpha,
                                      const cusparseMatDescr_t  descrA,
                                      const float            *csrValA,
                                      const int              *csrRowPtrA,
                                      const int              *csrColIndA,
                                      const float            *x,
                                      const float            *beta,
                                      float                  *y,
                                      void                     *buffer) {

// #ifdef DFLOAT
//   return cusparseCsrmvEx(handle, alg, CUSPARSE_OPERATION_NON_TRANSPOSE, m, n, nnz, alpha,
//                          CUDA_R_32F, descrA, csrValA, CUDA_R_32F, csrRowPtrA, csrColIndA, x,
//                          CUDA_R_32F, beta, CUDA_R_32F, y, CUDA_R_32F, CUDA_R_32F, buffer);
// #else
  return cusparseCsrmvEx(handle, alg, CUSPARSE_OPERATION_NON_TRANSPOSE, m, n, nnz, alpha,
                         CUDA_R_64F, descrA, csrValA, CUDA_R_64F, csrRowPtrA, csrColIndA, x,
                         CUDA_R_64F, beta, CUDA_R_64F, y, CUDA_R_64F, CUDA_R_64F, buffer);
// #endif
} 


static cusparseStatus_t cusparseCsrmv_bufferSize(cusparseHandle_t          handle,
                                                 cusparseAlgMode_t         alg,
                                                 int                     m,
                                                 int                     n,
                                                 int                     nnz,
                                                 const float            *alpha,
                                                 const cusparseMatDescr_t  descrA,
                                                 const float            *csrValA,
                                                 const int              *csrRowPtrA,
                                                 const int              *csrColIndA,
                                                 const float            *x,
                                                 const float            *beta,
                                                 float                  *y,
                                                 size_t                   *bufferSizeInBytes) {

// #ifdef DFLOAT
//   return cusparseCsrmvEx_bufferSize(handle, alg, CUSPARSE_OPERATION_NON_TRANSPOSE, m, n, nnz, alpha,
                                    CUDA_R_32F, descrA, csrValA, CUDA_R_32F, csrRowPtrA, csrColIndA, x,
                                    CUDA_R_32F, beta, CUDA_R_32F, y, CUDA_R_32F, CUDA_R_32F, bufferSizeInBytes);
// #else
  return cusparseCsrmvEx_bufferSize(handle, alg, CUSPARSE_OPERATION_NON_TRANSPOSE, m, n, nnz, alpha,
                                    CUDA_R_64F, descrA, csrValA, CUDA_R_64F, csrRowPtrA, csrColIndA, x,
                                    CUDA_R_64F, beta, CUDA_R_64F, y, CUDA_R_64F, CUDA_R_64F, bufferSizeInBytes);
#endif 
}


#endif /* ifndef CUDA_WRAPPER */