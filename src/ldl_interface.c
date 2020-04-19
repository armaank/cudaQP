#include "../include/glob_opts.h"

#include "../include/ldl.h"
#include "../include/ldl_interface.h"
#include "../include/kkt.h"


// Free LDL Factorization structure
void free_linsys_solver_ldl(ldl_solver *s) {
    if (s) {
        if (s->L)           csc_spfree(s->L);
        if (s->P)           c_free(s->P);
        if (s->Dinv)        c_free(s->Dinv);
        if (s->bp)          c_free(s->bp);
        if (s->sol)         c_free(s->sol);
        if (s->rho_inv_vec) c_free(s->rho_inv_vec);

        // These are required for matrix updates
        if (s->Pdiag_idx) c_free(s->Pdiag_idx);
        if (s->KKT)       csc_spfree(s->KKT);
        if (s->PtoKKT)    c_free(s->PtoKKT);
        if (s->AtoKKT)    c_free(s->AtoKKT);
        if (s->rhotoKKT)  c_free(s->rhotoKKT);

        // QDLDL workspace
        if (s->D)         c_free(s->D);
        if (s->etree)     c_free(s->etree);
        if (s->Lnz)       c_free(s->Lnz);
        if (s->iwork)     c_free(s->iwork);
        if (s->bwork)     c_free(s->bwork);
        if (s->fwork)     c_free(s->fwork);
        c_free(s);

    }
}


/**
 * Compute LDL factorization of matrix A
 * @param  A    Matrix to be factorized
 * @param  p    Private workspace
 * @param  nvar Number of QP variables
 * @return      exitstatus (0 is good)
 */
static int LDL_factor(csc *A,  ldl_solver * p, int nvar){

    int sum_Lnz;
    int factor_status;

    // Compute elimination tree
    sum_Lnz = QDLDL_etree(A->n, A->p, A->i, p->iwork, p->Lnz, p->etree);

    if (sum_Lnz < 0){
      // Error
#ifdef PRINTING
      c_eprint("Error in KKT matrix LDL factorization when computing the elimination tree. A is not perfectly upper triangular");
#endif
      return sum_Lnz;
    }

    // Allocate memory for Li and Lx
    p->L->i = (int *)c_malloc(sizeof(int)*sum_Lnz);
    p->L->x = (float *)c_malloc(sizeof(float)*sum_Lnz);

    // Factor matrix
    factor_status = QDLDL_factor(A->n, A->p, A->i, A->x,
                                 p->L->p, p->L->i, p->L->x,
                                 p->D, p->Dinv, p->Lnz,
                                 p->etree, p->bwork, p->iwork, p->fwork);


    if (factor_status < 0){
      // Error
#ifdef PRINTING
      c_eprint("Error in KKT matrix LDL factorization when computing the nonzero elements. There are zeros in the diagonal matrix");
#endif
      return factor_status;
    } else if (factor_status < nvar) {
      // Error: Number of positive elements of D should be equal to nvar
#ifdef PRINTING
      c_eprint("Error in KKT matrix LDL factorization when computing the nonzero elements. The problem seems to be non-convex");
#endif
      return -2;
    }

    return 0;

}


static int permute_KKT(csc ** KKT, ldl_solver * p, int Pnz, int Anz, int m, int * PtoKKT, int * AtoKKT, int * rhotoKKT){
    float *info;
    int amd_status;
    int * Pinv;
    csc *KKT_temp;
    int * KtoPKPt;
    int i; // Indexing

    info = (float *)c_malloc(AMD_INFO * sizeof(float));

    // Compute permutation matrix P using AMD
#ifdef DLONG
    amd_status = amd_l_order((*KKT)->n, (*KKT)->p, (*KKT)->i, p->P, (float *)OSQP_NULL, info);
#else
    amd_status = amd_order((*KKT)->n, (*KKT)->p, (*KKT)->i, p->P, (float *)OSQP_NULL, info);
#endif
    if (amd_status < 0) {
        // Free Amd info and return an error
        c_free(info);
        return amd_status;
    }


    // Inverse of the permutation vector
    Pinv = csc_pinv(p->P, (*KKT)->n);

    // Permute KKT matrix
    if (!PtoKKT && !AtoKKT && !rhotoKKT){  // No vectors to be stored
        // Assign values of mapping
        KKT_temp = csc_symperm((*KKT), Pinv, OSQP_NULL, 1);
    }
    else {
        // Allocate vector of mappings from unpermuted to permuted
        KtoPKPt = c_malloc((*KKT)->p[(*KKT)->n] * sizeof(int));
        KKT_temp = csc_symperm((*KKT), Pinv, KtoPKPt, 1);

        // Update vectors PtoKKT, AtoKKT and rhotoKKT
        if (PtoKKT){
            for (i = 0; i < Pnz; i++){
                PtoKKT[i] = KtoPKPt[PtoKKT[i]];
            }
        }
        if (AtoKKT){
            for (i = 0; i < Anz; i++){
                AtoKKT[i] = KtoPKPt[AtoKKT[i]];
            }
        }
        if (rhotoKKT){
            for (i = 0; i < m; i++){
                rhotoKKT[i] = KtoPKPt[rhotoKKT[i]];
            }
        }

        // Cleanup vector of mapping
        c_free(KtoPKPt);
    }

    // Cleanup
    // Free previous KKT matrix and assign pointer to new one
    csc_spfree((*KKT));
    (*KKT) = KKT_temp;
    // Free Pinv
    c_free(Pinv);
    // Free Amd info
    c_free(info);

    return 0;
}


// Initialize LDL Factorization structure
int init_linsys_solver_ldl(ldl_solver ** sp, const csc * P, const csc * A, float sigma, const float * rho_vec, int polish){

    // Define Variables
    csc * KKT_temp;     // Temporary KKT pointer
    int i;            // Loop counter
    int n_plus_m;     // Define n_plus_m dimension

    // Allocate private structure to store KKT factorization
    ldl_solver *s;
    s = c_calloc(1, sizeof(ldl_solver));
    *sp = s;

    // Size of KKT
    s->n = P->n;
    s->m = A->m;
    n_plus_m = s->n + s->m;

    // Sigma parameter
    s->sigma = sigma;

    // Polishing flag
    s->polish = polish;

    // Link Functions
    s->solve = &solve_linsys_ldl;

#ifndef EMBEDDED
    s->free = &free_linsys_solver_ldl;
#endif

#if EMBEDDED != 1
    s->update_matrices = &update_linsys_solver_matrices_ldl;
    s->update_rho_vec = &update_linsys_solver_rho_vec_ldl;
#endif

    // Assign type
    s->type = QDLDL_SOLVER;

    // Set number of threads to 1 (single threaded)
    s->nthreads = 1;

    // Sparse matrix L (lower triangular)
    // NB: We don not allocate L completely (CSC elements)
    //      L will be allocated during the factorization depending on the
    //      resulting number of elements.
    s->L = c_malloc(sizeof(csc));
    s->L->m = n_plus_m;
    s->L->n = n_plus_m;
    s->L->nz = -1;

    // Diagonal matrix stored as a vector D
    s->Dinv = (QDLDL_float *)c_malloc(sizeof(QDLDL_float) * n_plus_m);
    s->D    = (QDLDL_float *)c_malloc(sizeof(QDLDL_float) * n_plus_m);

    // Permutation vector P
    s->P    = (QDLDL_int *)c_malloc(sizeof(QDLDL_int) * n_plus_m);

    // Working vector
    s->bp   = (QDLDL_float *)c_malloc(sizeof(QDLDL_float) * n_plus_m);

    // Solution vector
    s->sol  = (QDLDL_float *)c_malloc(sizeof(QDLDL_float) * n_plus_m);

    // Parameter vector
    s->rho_inv_vec = (float *)c_malloc(sizeof(float) * s->m);

    // Elimination tree workspace
    s->etree = (QDLDL_int *)c_malloc(n_plus_m * sizeof(QDLDL_int));
    s->Lnz   = (QDLDL_int *)c_malloc(n_plus_m * sizeof(QDLDL_int));

    // Preallocate L matrix (Lx and Li are sparsity dependent)
    s->L->p = (int *)c_malloc((n_plus_m+1) * sizeof(QDLDL_int));

    // Lx and Li are sparsity dependent, so set them to
    // null initially so we don't try to free them prematurely
    s->L->i = OSQP_NULL;
    s->L->x = OSQP_NULL;

    // Preallocate workspace
    s->iwork = (QDLDL_int *)c_malloc(sizeof(QDLDL_int)*(3*n_plus_m));
    s->bwork = (QDLDL_bool *)c_malloc(sizeof(QDLDL_bool)*n_plus_m);
    s->fwork = (QDLDL_float *)c_malloc(sizeof(QDLDL_float)*n_plus_m);

    // Form and permute KKT matrix
    if (polish){ // Called from polish()
        // Use s->rho_inv_vec for storing param2 = vec(delta)
        for (i = 0; i < A->m; i++){
            s->rho_inv_vec[i] = sigma;
        }

        KKT_temp = form_KKT(P, A, 0, sigma, s->rho_inv_vec, OSQP_NULL, OSQP_NULL, OSQP_NULL, OSQP_NULL, OSQP_NULL);

        // Permute matrix
        if (KKT_temp)
            permute_KKT(&KKT_temp, s, OSQP_NULL, OSQP_NULL, OSQP_NULL, OSQP_NULL, OSQP_NULL, OSQP_NULL);
    }
    else { // Called from ADMM algorithm

        // Allocate vectors of indices
        s->PtoKKT = c_malloc((P->p[P->n]) * sizeof(int));
        s->AtoKKT = c_malloc((A->p[A->n]) * sizeof(int));
        s->rhotoKKT = c_malloc((A->m) * sizeof(int));

        // Use p->rho_inv_vec for storing param2 = rho_inv_vec
        for (i = 0; i < A->m; i++){
            s->rho_inv_vec[i] = 1. / rho_vec[i];
        }

        KKT_temp = form_KKT(P, A, 0, sigma, s->rho_inv_vec,
                            s->PtoKKT, s->AtoKKT,
                            &(s->Pdiag_idx), &(s->Pdiag_n), s->rhotoKKT);

        // Permute matrix
        if (KKT_temp)
            permute_KKT(&KKT_temp, s, P->p[P->n], A->p[A->n], A->m, s->PtoKKT, s->AtoKKT, s->rhotoKKT);
    }

    // Check if matrix has been created
    if (!KKT_temp){
#ifdef PRINTING
        c_eprint("Error forming and permuting KKT matrix");
#endif
        free_linsys_solver_ldl(s);
        *sp = OSQP_NULL;
        return OSQP_LINSYS_SOLVER_INIT_ERROR;
    }

    // Factorize the KKT matrix
    if (LDL_factor(KKT_temp, s, P->n) < 0) {
        csc_spfree(KKT_temp);
        free_linsys_solver_ldl(s);
        *sp = OSQP_NULL;
        return OSQP_NONCVX_ERROR;
    }

    if (polish){ // If KKT passed, assign it to KKT_temp
        // Polish, no need for KKT_temp
        csc_spfree(KKT_temp);
    }
    else { // If not embedded option 1 copy pointer to KKT_temp. Do not free it.
        s->KKT = KKT_temp;
    }


    // No error
    return 0;
}

#endif  // EMBEDDED


// Permute x = P*b using P
void permute_x(int n, float * x, const float * b, const int * P) {
    int j;
    for (j = 0 ; j < n ; j++) x[j] = b[P[j]];
}

// Permute x = P'*b using P
void permutet_x(int n, float * x, const float * b, const int * P) {
    int j;
    for (j = 0 ; j < n ; j++) x[P[j]] = b[j];
}


static void LDLSolve(float *x, float *b, const csc *L, const float *Dinv, const int *P, float *bp) {
    /* solves P'LDL'P x = b for x */
    permute_x(L->n, bp, b, P);
    QDLDL_solve(L->n, L->p, L->i, L->x, Dinv, bp);
    permutet_x(L->n, x, bp, P);

}


int solve_linsys_ldl(ldl_solver * s, float * b) {
    int j;

#ifndef EMBEDDED
    if (s->polish) {
        /* stores solution to the KKT system in b */
        LDLSolve(b, b, s->L, s->Dinv, s->P, s->bp);
    } else {
#endif
        /* stores solution to the KKT system in s->sol */
        LDLSolve(s->sol, b, s->L, s->Dinv, s->P, s->bp);

        /* copy x_tilde from s->sol */
        for (j = 0 ; j < s->n ; j++) {
            b[j] = s->sol[j];
        }

        /* compute z_tilde from b and s->sol */
        for (j = 0 ; j < s->m ; j++) {
            b[j + s->n] += s->rho_inv_vec[j] * s->sol[j + s->n];
        }
#ifndef EMBEDDED
    }
#endif

    return 0;
}


#if EMBEDDED != 1
// Update private structure with new P and A
int update_linsys_solver_matrices_ldl(ldl_solver * s, const csc *P, const csc *A) {

    // Update KKT matrix with new P
    update_KKT_P(s->KKT, P, s->PtoKKT, s->sigma, s->Pdiag_idx, s->Pdiag_n);

    // Update KKT matrix with new A
    update_KKT_A(s->KKT, A, s->AtoKKT);

    return (QDLDL_factor(s->KKT->n, s->KKT->p, s->KKT->i, s->KKT->x,
        s->L->p, s->L->i, s->L->x, s->D, s->Dinv, s->Lnz,
        s->etree, s->bwork, s->iwork, s->fwork) < 0);

}


int update_linsys_solver_rho_vec_ldl(ldl_solver * s, const float * rho_vec){
    int i;

    // Update internal rho_inv_vec
    for (i = 0; i < s->m; i++){
        s->rho_inv_vec[i] = 1. / rho_vec[i];
    }

    // Update KKT matrix with new rho_vec
    update_KKT_param2(s->KKT, s->rho_inv_vec, s->rhotoKKT, s->m);

    return (QDLDL_factor(s->KKT->n, s->KKT->p, s->KKT->i, s->KKT->x,
        s->L->p, s->L->i, s->L->x, s->D, s->Dinv, s->Lnz,
        s->etree, s->bwork, s->iwork, s->fwork) < 0);
}


#endif
