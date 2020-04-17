#ifndef KKT_H
# define KKT_H

#  include "cs.h"

/* form square symmetric KKT matrix */
csc* form_KKT(const csc *P, const csc *A, int format, float param1, float *param2, int *PtoKKT, int *AtoKKT, int **Pdiag_idx, int *Pdiag_n, int *param2toKKT);

/* update KKT matrix using the elements of P */
void update_KKT_P(csc *KKT, const csc *P, const int *PtoKKT, const float param1, const int *Pdiag_idx, const int Pdiag_n);

/* update KKT matrix using the elements of A */
void update_KKT_A(csc *KKT, const csc *A, const int *AtoKKT);

/* update KKT matrix with new param2 */
void update_KKT_param2(csc *KKT, const float *param2, const int *param2toKKT, const int m);

#endif // ifndef KKT_H
