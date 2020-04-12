#ifndef TYPES_H
# define TYPES_H

/* Use customized number representation */
# ifdef DLONG            // long integers
typedef long long c_int; /* for indices */
# else // standard integers
typedef int c_int;       /* for indices */
# endif /* ifdef DLONG */


# ifndef DFLOAT         // Doubles
typedef double c_float; /* for numerical values  */
# else                  // Floats
typedef float c_float;  /* for numerical values  */
# endif /* ifndef DFLOAT */

/*
 * matrix in compressed-column form.
 */
typedef struct {
    c_int    nzmax; ///< maximum number of entries
    c_int    m;     ///< number of rows
    c_int    n;     ///< number of columns
    c_int   *p;     ///< column pointers (size n+1); col indices (size nzmax) start from 0 when using triplet format (direct KKT matrix formation)
    c_int   *i;     ///< row indices, size nzmax starting from 0
    c_float *x;     ///< numerical values, size nzmax
    c_int    nz;    ///< number of entries in triplet matrix, -1 for csc
} csc;

# endif /* ifndef TYPES_H */
