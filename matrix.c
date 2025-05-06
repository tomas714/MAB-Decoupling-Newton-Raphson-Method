#include "matrix.h"

#include <string.h>
#include <math.h>

static void matrix_chol(struct matrix* l, const struct matrix* a)
{
    // Full rank Cholesky factorization: A = L * L^T

    unsigned int size = a->rows;
    memset(l, 0, sizeof(struct matrix));
    l->rows = size;
    l->cols = size;

    float tol = INFINITY;
    for (unsigned int i = 0; i < size; i++) {
        float val = a->data[i][i];
        if ((val > 0.0F) && (val < tol)) {
            tol = val;
        }
    }
    tol *= 1e-6F;

    unsigned int r = 0;
    for (unsigned int k = 0; k < size; k++) {
        for (unsigned int i = k; i < size; i++) {
            float sum = 0.0F;
            for (unsigned int j = 0; j < r; j++) {
                sum += l->data[i][j] * l->data[k][j];
            }
            l->data[i][r] = a->data[i][k] - sum;
        }

        if (l->data[k][r] > tol) {
            l->data[k][r] = sqrtf(l->data[k][r]);
            for (unsigned int i = k + 1; i < size; i++) {
                l->data[i][r] = l->data[i][r] / l->data[k][r];
            }
            r++;
        }
    }

    l->cols = r;
}

void matrix_inv(struct matrix* res, struct matrix* a)
{
    static struct matrix tmp;
    matrix_chol(&tmp, a);

    unsigned int size = a->rows;
    memset(res, 0, sizeof(struct matrix));
    res->rows = size;
    res->cols = size;

    float sum;

    // Forward substitution : res = L^(-1)
    for (unsigned int i = 0; i < size; i++) {
        res->data[i][i] = 1.0F / tmp.data[i][i];
        for (unsigned int j = i + 1; j < size; j++) {
            sum = 0.0F;
            for (unsigned int k = i; k < j; k++) {
                sum -= tmp.data[j][k] * res->data[k][i];
            }
            res->data[j][i] = sum / tmp.data[j][j];
        }
    }

    // Inversion : A^(-1) = (L^(-1)') * L^(-1)
    for (int i = 0; i < size; i++) {
       for (int j = i; j < size; j++) {
           sum = 0.0F;
           for (int k = j; k < size; k++) {
               sum += res->data[k][i] * res->data[k][j];
           }
           res->data[i][j] = sum;
           res->data[j][i] = sum;
       }
    }
}

void matrix_multiple(struct matrix* res, struct matrix* a, struct matrix* b)
{
    unsigned int i_max = a->rows;
    unsigned int j_max = b->cols;
    unsigned int k_max = a->cols;

    for (unsigned int i = 0; i < i_max; i++) {
        for (unsigned int j = 0; j < j_max; j++) {
            res->data[i][j] = 0.0F;
            for (unsigned int k = 0; k < k_max; k++) {
                res->data[i][j] += a->data[i][k] * b->data[k][j];
            }
        }
    }
    res->rows = i_max;
    res->cols = j_max;
}

void matrix_transpose(struct matrix* res, struct matrix* a)
{
    unsigned int i_max = a->rows;
    unsigned int j_max = a->cols;

    for (unsigned int i = 0; i < i_max; i++) {
        for (unsigned int j = 0; j < j_max; j++) {
            res->data[j][i] = a->data[i][j];
        }
    }
    res->rows = a->cols;
    res->cols = a->rows;
}

void matrix_geninv(struct matrix* res, struct matrix* a)
{
    static struct matrix  tmp[3];
    struct matrix* ptr[5] = { &tmp[0], &tmp[1], &tmp[2], a, res};

    matrix_transpose(ptr[0], a             );
    matrix_multiple (ptr[1], ptr[0], a     );
    matrix_chol     (ptr[2], ptr[1]        );
    matrix_transpose(ptr[1], ptr[2]        );
    matrix_multiple (ptr[3], ptr[1], ptr[2]);
    matrix_inv      (ptr[4], ptr[3]        );
    matrix_multiple (ptr[3], ptr[2], ptr[4]);
    matrix_multiple (ptr[2], ptr[3], ptr[4]);
    matrix_multiple (ptr[3], ptr[2], ptr[1]);
    matrix_multiple (ptr[4], ptr[3], ptr[0]);
}
