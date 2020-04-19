#include "./include/lin_sys.h"

#include "./include/ldl.h" // Include only this solver in the same directory

const char *LINSYS_SOLVER_NAME[] = {
  "qdldl"
};


// Load linear system solver shared library
c_int load_linsys_solver(enum linsys_solver_type linsys_solver) {
  switch (linsys_solver) {
  case QDLDL_SOLVER:

    return 0;

  default: // QDLDL
    return 0;
  }
}

// Unload linear system solver shared library
c_int unload_linsys_solver(enum linsys_solver_type linsys_solver) {
  switch (linsys_solver) {
  case QDLDL_SOLVER:

    // We do not load QDLDL solver. We have the source.
    return 0;

  default: //  QDLDL
    return 0;
  }
}

// Initialize linear system solver structure
// NB: Only the upper triangular part of P is stuffed!
c_int init_linsys_solver(LinSysSolver          **s,
                         const csc              *P,
                         const csc              *A,
                         c_float                 sigma,
                         const c_float          *rho_vec,
                         enum linsys_solver_type linsys_solver,
                         c_int                   polish) {
  switch (linsys_solver) {
  case QDLDL_SOLVER:
    return init_linsys_solver_qdldl((qdldl_solver **)s, P, A, sigma, rho_vec, polish);

  default: // QDLDL
    return init_linsys_solver_qdldl((qdldl_solver **)s, P, A, sigma, rho_vec, polish);
  }
}