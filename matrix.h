#ifndef MATRIX_H
#define MATRIX_H

#include "config.h"

struct matrix {
    unsigned int rows;
    unsigned int cols;
    float data[MAB_PORTS][MAB_PORTS];
};

void matrix_transpose(struct matrix* res, struct matrix* a);
void matrix_multiple(struct matrix* res, struct matrix* a, struct matrix* b);
void matrix_geninv(struct matrix* res, struct matrix* a);
void matrix_inv(struct matrix* res, struct matrix* a);

#endif /* MATRIX_H */
