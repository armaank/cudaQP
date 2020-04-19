#ifndef LDL_INTERFACE_H
#define LDL_INTERFACE_H
#include <stdbool.h>
#include "ldl.h"
#include "csc.h"
/**
 * QDLDL solver structure
 */
typedef struct ldl ldl_solver;

struct ldl {
    enum linsys_solver_type type;

    /**
     * @name Functions
     * @{
     */
    int (*solve)(struct ldl * self, float * b);

    void (*free)(struct ldl * self); ///< Free workspace (only if desktop)

    // This used only in non embedded or embedded 2 version
    int (*update_matrices)(struct ldl * self, const csc *P, const csc *A);  ///< Update solver matrices
    int (*update_rho_vec)(struct ldl * self, const float * rho_vec);      ///< Update rho_vec parameter

    int nthreads;
    /** @} */

    /**
     * @name Attributes
     * @{
     */
    csc *L;                 ///< lower triangular matrix in LDL factorization
    float *Dinv;          ///< inverse of diag matrix in LDL (as a vector)
    int   *P;             ///< permutation of KKT matrix for factorization
    float *bp;            ///< workspace memory for solves
    float *sol;           ///< solution to the KKT system
    float *rho_inv_vec;   ///< parameter vector
    float sigma;          ///< scalar parameter

    int polish;           ///< polishing flag

    int n;                ///< number of QP variables
    int m;                ///< number of QP constraints


    // These are required for matrix updates
    int * Pdiag_idx, Pdiag_n;  ///< index and number of diagonal elements in P
    csc   * KKT;                 ///< Permuted KKT matrix in sparse form (used to update P and A matrices)
    int * PtoKKT, * AtoKKT;    ///< Index of elements from P and A to KKT matrix
    int * rhotoKKT;            ///< Index of rho places in KKT matrix
    // QDLDL Numeric workspace
    float *D;
    int   *etree;
    int   *Lnz;
    int   *iwork;
    bool  *bwork;
    float *fwork;

    /** @} */
};



/**
 * Initialize QDLDL Solver
 *
 * @param  s         Pointer to a private structure
 * @param  P         Cost function matrix (upper triangular form)
 * @param  A         Constraints matrix
 * @param  sigma     Algorithm parameter. If polish, then sigma = delta.
 * @param  rho_vec   Algorithm parameter. If polish, then rho_vec = OSQP_NULL.
 * @param  polish    Flag whether we are initializing for polish or not
 * @return           Exitflag for error (0 if no errors)
 */
int init_linsys_solver_ldl(ldl_solver ** sp, const csc * P, const csc * A, float sigma, const float * rho_vec, int polish);

/**
 * Solve linear system and store result in b
 * @param  s        Linear system solver structure
 * @param  b        Right-hand side
 * @return          Exitflag
 */
int solve_linsys_ldl(ldl_solver * s, float * b);


/**
 * Update linear system solver matrices
 * @param  s        Linear system solver structure
 * @param  P        Matrix P
 * @param  A        Matrix A
 * @return          Exitflag
 */
int update_linsys_solver_matrices_ldl(ldl_solver * s, const csc *P, const csc *A);




/**
 * Update rho_vec parameter in linear system solver structure
 * @param  s        Linear system solver structure
 * @param  rho_vec  new rho_vec value
 * @return          exitflag
 */
int update_linsys_solver_rho_vec_ldl(ldl_solver * s, const float * rho_vec);

/**
 * Free linear system solver
 * @param s linear system solver object
 */
void free_linsys_solver_ldl(ldl_solver * s);

