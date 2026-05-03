#include "matrix.h"

void matrix_geninv(struct matrix* restrict res, const struct matrix* a)
{
    static struct matrix tmp[3];

    matrix_gram         (&tmp[1], a);
    matrix_cstmCholesky (&tmp[0], &tmp[1]);
    matrix_gram         (&tmp[1], &tmp[0]);
    matrix_cholesky     (&tmp[2], &tmp[1]);
    matrix_forwardSubs  (&tmp[1], &tmp[2]);
    matrix_gram         (&tmp[2], &tmp[1]);
    matrix_multiply     (&tmp[1], &tmp[2], &tmp[2]);
    matrix_multiply     (&tmp[2], &tmp[0], &tmp[1]);
    matrix_multiplyTranB(&tmp[1], &tmp[2], &tmp[0]);
    matrix_multiplyTranB(res,     &tmp[1], a);
}
