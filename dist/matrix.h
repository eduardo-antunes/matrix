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

// Initializers, print and free:

void matrix_init(Matrix *mat, int rows, int cols);

void matrix_init_id(Matrix *mat, int order);

void matrix_init_with(Matrix *mat, int rows, int cols, double p[rows][cols]);

void matrix_init_with_diag(Matrix *mat, int order, double diag[order]);

void matrix_init_copy(Matrix *mat, const Matrix *src);

void matrix_print(const Matrix *mat);

void matrix_free(Matrix *mat);

// Basic arithmetic operations:

void matrix_scale(Matrix *mat, double k);

Matrix matrix_add(const Matrix *a, const Matrix *b);

Matrix matrix_prod(const Matrix *a, const Matrix *b);

Matrix matrix_transpose(const Matrix *mat);

// Tests and conditions:

bool matrix_is_square(const Matrix *mat);

bool matrix_is_diagonally_dominant(const Matrix *mat);

// Solving linear systems of equations:

void matrix_reduce(Matrix *mat);

Matrix matrix_augment(const Matrix *a, const Matrix *b);

Matrix matrix_inverse(const Matrix *mat);

Matrix matrix_solve(Matrix *a, Matrix *b, bool *solution);

Matrix matrix_solve_const(const Matrix *a, const Matrix *b, bool *solution);

#define DEFAULT_GAUSS_SEIDEL_ITERS 64

/* This implements the Gauss-Seidel algorithm. As such, it requires
 * matrix A to be a square matrix. To guarantee convergence, it also
 * must be strictly diagonally dominant.
 */
Matrix matrix_solve_numerical(const Matrix *a, const Matrix *b, int iters);

#endif // MATRIX_H
