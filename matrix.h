#ifndef MATRIX_H
#define MATRIX_H

#include "config.h"
#include <string.h>
#include <math.h>

struct matrix {
    unsigned int rows;
    unsigned int cols;
    float data[MAB_PORTS][MAB_PORTS];
};

static inline void matrix_multiply(struct matrix* restrict res, const struct matrix* a, const struct matrix* b)
{
    unsigned int i_max = a->rows;
    unsigned int j_max = b->cols;
    unsigned int k_max = a->cols;

    for (unsigned int i = 0; i < i_max; i++) {
        for (unsigned int j = 0; j < j_max; j++) {
            float sum = 0.0F;
            for (unsigned int k = 0; k < k_max; k++) {
                sum += a->data[i][k] * b->data[k][j];
            }
            res->data[i][j] = sum;
        }
    }
    res->rows = i_max;
    res->cols = j_max;
}

static inline void matrix_multiplyTranB(struct matrix* restrict res, const struct matrix* a, const struct matrix* b)
{
    unsigned int i_max = a->rows;
    unsigned int j_max = b->rows;
    unsigned int k_max = a->cols;

    for (unsigned int i = 0; i < i_max; i++) {
        for (unsigned int j = 0; j < j_max; j++) {
            float sum = 0.0F;
            for (unsigned int k = 0; k < k_max; k++) {
                sum += a->data[i][k] * b->data[j][k];
            }
            res->data[i][j] = sum;
        }
    }
    res->rows = i_max;
    res->cols = j_max;
}

static inline void matrix_gram(struct matrix* restrict res, const struct matrix* a)
{
    unsigned int rows = a->rows;
    unsigned int cols = a->cols;

    for (unsigned int i = 0; i < cols; i++) {
        for (unsigned int j = i; j < cols; j++) {
            float sum = 0.0F;
            for (unsigned int k = 0; k < rows; k++) {
                sum += a->data[k][i] * a->data[k][j];
            }
            res->data[i][j] = sum;
            res->data[j][i] = sum;
        }
    }

    res->rows = cols;
    res->cols = cols;
}

static inline void matrix_cstmCholesky(struct matrix* restrict res, const struct matrix* a)
{
    unsigned int size = a->rows;

    memset(res, 0, sizeof(struct matrix));
    for (unsigned int i = 0; i < (size - 1); i++)
    {
        float sum = 0.0f;
        for (unsigned int k = 0; k < i; k++) {
            sum += res->data[i][k] * res->data[i][k];
        }

        float val = a->data[i][i] - sum;
        res->data[i][i] = (val > 1e-6F) ? sqrtf(val) : 1e-6F;

        float invDiag = 1.0f / res->data[i][i];
        for (unsigned int j = i + 1; j < size; j++) {
            sum = 0.0F;
            for (unsigned int k = 0; k < i; k++) {
                sum += res->data[i][k] * res->data[j][k];
            }
            res->data[j][i] = (a->data[i][j] - sum) * invDiag;
        }
    }

    res->rows = size;
    res->cols = size - 1;
}

static inline void matrix_cholesky(struct matrix* restrict res, const struct matrix* a)
{
    unsigned int size = a->rows;

    memset(res, 0, sizeof(struct matrix));
    for (unsigned int i = 0; i < size; i++)
    {
        float sum = 0.0f;
        for (unsigned int k = 0; k < i; k++) {
            sum += res->data[i][k] * res->data[i][k];
        }

        float val = a->data[i][i] - sum;
        res->data[i][i] = (val > 1e-6F) ? sqrtf(val) : 1e-6F;

        float invDiag = 1.0f / res->data[i][i];
        for (unsigned int j = i + 1; j < size; j++) {
            sum = 0.0F;
            for (unsigned int k = 0; k < i; k++) {
                sum += res->data[i][k] * res->data[j][k];
            }
            res->data[j][i] = (a->data[i][j] - sum) * invDiag;
        }
    }

    res->rows = size;
    res->cols = size;
}

static inline void matrix_forwardSubs(struct matrix* restrict res, const struct matrix* a)
{
    unsigned int size = a->rows;

    memset(res, 0, sizeof(struct matrix));
    for (unsigned int i = 0; i < size; i++) {
        res->data[i][i] = 1.0F / a->data[i][i];
        for (unsigned int j = i + 1; j < size; j++) {
            float sum = 0.0F;
            for (unsigned int k = i; k < j; k++) {
                sum -= a->data[j][k] * res->data[k][i];
            }
            res->data[j][i] = sum / a->data[j][j];
        }
    }

    res->rows = size;
    res->cols = size;
}

void matrix_geninv(struct matrix* restrict res, const struct matrix* a);

#endif /* MATRIX_H */
