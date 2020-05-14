/* constants for cuda_pgh */

#ifndef CUDA_PCG_CONSTANTS_H
# define CUDA_PCG_CONSTANTS_H

#ifdef __cplusplus
extern "C" {
#endif

/* PCG parameters */
#define CUDA_PCG_PRECONDITION        (1)
#define CUDA_PCG_MAX_ITER            (20)
#define CUDA_PCG_WARM_START          (1)
#define CUDA_PCG_NORM                (0)     /* 0: inf;  2: Euclidean */
#define CUDA_PCG_EPS_MIN             (1e-7)

/* tolerance parameters */
#define CUDA_PCG_START_TOL           (50)
#define CUDA_PCG_DECAY_RATE          (2.75)
#define CUDA_PCG_REDUCTION_FACTOR    (0.15)
#define CUDA_PCG_REDUCTION_THRESHOLD (10)

/* polish parameters */
#define CUDA_PCG_POLISH_ACCURACY     (1e-5)
#define CUDA_PCG_POLISH_MAX_ITER     (1e3)

/* pcg tolerance stragegies*/
enum pcg_eps_strategy { SCS_STRATEGY, RESIDUAL_STRATEGY };


#ifdef __cplusplus
}
#endif

#endif /* #ifndef CUDA_PCG_CONSTANTS_H */
