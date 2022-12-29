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

void matrix_init_with(Matrix *mat, int rows, int cols, double p[rows][cols]);

void matrix_init_id(Matrix *mat, int order);

Matrix matrix_copy(Matrix *mat);

void matrix_free(Matrix *mat);

void matrix_print(const Matrix *mat);

// Basic operations:

bool matrix_square(Matrix *mat);

void matrix_scale(Matrix *mat, double k);

Matrix matrix_add(Matrix *a, Matrix *b);

Matrix matrix_prod(Matrix *a, Matrix *b);

// Gaussian reduction operations:

void matrix_swap_rows(Matrix *mat, int r1, int r2);

void matrix_scale_row(Matrix *mat, int r, double k);

void matrix_divide_row(Matrix *mat, int r, double k);

void matrix_add_rows(Matrix *mat, int rr, int rf, double k);

void matrix_reduce(Matrix *mat);

// System solving operations:

Matrix matrix_augment(Matrix *a, Matrix *b);

Matrix matrix_solve(Matrix *a, Matrix *b, bool *has_solution);

Matrix matrix_solve_aug(Matrix *sys);

Matrix matrix_inverse(Matrix *mat);

#endif // MATRIX_H
