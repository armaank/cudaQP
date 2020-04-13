#ifndef KKT_H
# define KKT_H

# include "types.h"

#  include "cs.h"

/* form square symmetric KKT matrix */
csc* form_KKT(const csc *P, const csc *A, c_int format, c_float param1, c_float *param2, c_int *PtoKKT, c_int *AtoKKT, c_int **Pdiag_idx, c_int *Pdiag_n, c_int *param2toKKT);

/* update KKT matrix using the elements of P */
void update_KKT_P(csc *KKT, const csc *P, const c_int *PtoKKT, const c_float param1, const c_int *Pdiag_idx, const c_int Pdiag_n);

/* update KKT matrix using the elements of A */
void update_KKT_A(csc *KKT, const csc *A, const c_int *AtoKKT);

/* update KKT matrix with new param2 */
void update_KKT_param2(csc *KKT, const c_float *param2, const c_int *param2toKKT, const c_int m);

#endif // ifndef KKT_H
