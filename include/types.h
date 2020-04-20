typedef struct {
  int    nzmax; ///< maximum number of entries
  int    m;     ///< number of rows
  int    n;     ///< number of columns
  int   *p;     ///< column pointers (size n+1); col indices (size nzmax) start from 0 when using triplet format (direct KKT matrix formation)
  int   *i;     ///< row indices, size nzmax starting from 0
  float *x;     ///< numerical values, size nzmax
  int    nz;    ///< number of entries in triplet matrix, -1 for csc
} csc;
