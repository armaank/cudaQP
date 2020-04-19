#ifndef GLOB_OPTS_H
# define GLOB_OPTS_H


/*
   Define OSQP compiler flags
 */



/* DATA CUSTOMIZATIONS (depending on memory manager)-----------------------   */

/* If no custom memory allocator defined, use
 * standard linux functions. Custom memory allocator definitions
 * appear in the osqp_configure.h generated file. */
    #  include <stdlib.h>
    #  define c_malloc  malloc
    #  define c_calloc  calloc
    #  define c_free    free
    #  define c_realloc realloc




/* Use customized operations */

# ifndef c_absval
#  define c_absval(x) (((x) < 0) ? -(x) : (x))
# endif /* ifndef c_absval */

# ifndef c_max
#  define c_max(a, b) (((a) > (b)) ? (a) : (b))
# endif /* ifndef c_max */

# ifndef c_min
#  define c_min(a, b) (((a) < (b)) ? (a) : (b))
# endif /* ifndef c_min */

// Round x to the nearest multiple of N
# ifndef c_roundmultiple
#  define c_roundmultiple(x, N) ((x) + .5 * (N)-c_fmod((x) + .5 * (N), (N)))
# endif /* ifndef c_roundmultiple */


/* Use customized functions -----------------------------------------------   */


#  include <math.h>
#  ifndef DFLOAT // Doubles
#   define c_sqrt sqrt
#   define c_fmod fmod
#  else          // Floats
#   define c_sqrt sqrtf
#   define c_fmod fmodf
#  endif /* ifndef DFLOAT */


#  include <stdio.h>
#  include <string.h>



/* Print error macro */
#  define c_eprint(...) c_print("ERROR in %s: ", __FUNCTION__); c_print(\
    __VA_ARGS__); c_print("\n");




#endif /* ifndef GLOB_OPTS_H */
