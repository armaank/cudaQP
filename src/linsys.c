#include "../include/linsys.h"

const char *LINSYS_SOLVER_NAME[] = {
  "ldl"
};

// Load linear system solver shared library
int load_linsys_solver(enum linsys_solver_type linsys_solver) {
  switch (linsys_solver) {
  case LDL_SOLVER:

    // We do not load  QDLDL solver. We have the source.
    return 0;



  default: // QDLDL
    return 0;
  }
}

// Unload linear system solver shared library
int unload_linsys_solver(enum linsys_solver_type linsys_solver) {
  switch (linsys_solver) {
  case LDL_SOLVER:

    // We do not load QDLDL solver. We have the source.
    return 0;


  default: //  QDLDL
    return 0;
  }
}

// Initialize linear system solver structure
// NB: Only the upper triangular part of P is stuffed!
int init_linsys_solver(LinSysSolver          **s,
                         const csc              *P,
                         const csc              *A,
                         float                 sigma,
                         const float          *rho_vec,
                         enum linsys_solver_type linsys_solver,
                         int                   polish) {
  switch (linsys_solver) {
  case LDL_SOLVER:
    return init_linsys_solver_ldl((ldl_solver **)s, P, A, sigma, rho_vec, polish);

  default: // QDLDL
    return init_linsys_solver_ldl((ldl_solver **)s, P, A, sigma, rho_vec, polish);
  }
}