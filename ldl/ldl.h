#include <stdbool.h>

/* compute the elimination tree for a quasidefinite matrix */
int LDL_etree(const int n, const int* Ap, const int* Ai, int* work, int* Lnz, int* etree);

/* compute an LDL decomposition for a quasidefinite matrix */
int LDL_factor(const int n, const int* Ap, const int* Ai, const float* Ax, int* Lp, int* Li, float* Lx, float* D, float* Dinv, const int* Lnz, const int* etree, bool* bwork, int* iwork, float* fwork);

/* solves LDL'x = b */
void LDL_solve(const int n, const int* Lp, const int* Li, const float* Lx, const float* Dinv, float* x);

/* solves (L+I)x = b */
void LDL_Lsolve(const int n, const int* Lp, const int* Li, const float* Lx, float* x);

/* solves (L+I)'x = b */
void LDL_Ltsolve(const int n, const int* Lp, const int* Li, const float* Lx, float* x);

