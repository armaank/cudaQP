/* library handler for cuda */
#ifndef CUDA_HANDLER_H
# define CUDA_HANDLER_H

#include <cusparse_v2.h>
#include <cublas_v2.h>

# ifdef __cplusplus
extern "C" {
# endif

typedef struct {
    cublasHandle_t    cublasHandle;
    cusparseHandle_t  cusparseHandle;
    int              *d_index;
} CUDA_Handle_t;


CUDA_Handle_t* cuda_init_libs(void);


void cuda_free_libs(CUDA_Handle_t *CUDA_handle);


# ifdef __cplusplus
}
# endif

#endif /* ifndef CUDA_HANDLER_H */
