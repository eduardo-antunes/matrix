/* I dedicate any and all copyright interest in this software to the
 * public domain. I make this dedication for the benefit of the public at
 * large and to the detriment of my heirs and successors. I intend this
 * dedication to be an overt act of relinquishment in perpetuity of all
 * present and future rights to this software under copyright law.
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <stdbool.h>

typedef struct {
    int rows, cols;
    double **p;
} Matrix;

// Meta operations:

void matrix_init(Matrix *mat, int rows, int cols);

void matrix_init_id(Matrix *mat, int order);

void matrix_init_with(Matrix *mat, int rows, int cols, double p[rows][cols]);

void matrix_init_copy(Matrix *mat, const Matrix *src);

void matrix_print(const Matrix *mat);

void matrix_free(Matrix *mat);

// Basic arithmetic operations:

void matrix_scale(Matrix *mat, double k);

void matrix_add(const Matrix *a, const Matrix *b, Matrix *res);

void matrix_prod(const Matrix *a, const Matrix *b, Matrix *res);

// Solving linear systems of equations:

bool matrix_solve(Matrix *a, Matrix *b, Matrix *res);

#endif // MATRIX_H
