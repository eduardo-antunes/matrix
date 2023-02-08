/* I dedicate any and all copyright interest in this software to the
 * public domain. I make this dedication for the benefit of the public at
 * large and to the detriment of my heirs and successors. I intend this
 * dedication to be an overt act of relinquishment in perpetuity of all
 * present and future rights to this software under copyright law.
 */

#ifndef MATRIX_H
#define MATRIX_H

typedef struct {
    int rows, cols;
    double *p;
} Matrix;

// Initializers, print and free:

Matrix *matrix_alloc(int rows, int cols);

Matrix *matrix_calloc(int rows, int cols);

Matrix *matrix_identity(int order);

Matrix *matrix_copy(const Matrix *src);

void matrix_print(const Matrix *mat);

void matrix_free(Matrix *mat);

// Arithmetic operations:

void matrix_scale(Matrix *mat, double k);

Matrix *matrix_scale_copy(const Matrix *mat, double k);

Matrix *matrix_add(const Matrix *a, const Matrix *b);

Matrix *matrix_prod(const Matrix *a, const Matrix *b);

Matrix *matrix_transpose(const Matrix *mat);

// Tests and conditions:

int matrix_are_equal(const Matrix *a, const Matrix *b);

int matrix_is_square(const Matrix *mat);

int matrix_is_diagonally_dominant(const Matrix *mat);

// Solving linear systems of equations:

void matrix_reduce(Matrix *mat);

Matrix *matrix_reduce_const(const Matrix *mat);

Matrix *matrix_augment(const Matrix *a, const Matrix *b);

Matrix *matrix_inverse(const Matrix *mat);

Matrix *matrix_solve(Matrix *a, Matrix *b, int *solution);

Matrix *matrix_solve_copy(const Matrix *a, const Matrix *b, int *solution);

#define DEFAULT_GAUSS_SEIDEL_ITERS 64

/* This implements the Gauss-Seidel algorithm. As such, it requires
 * matrix A to be a square matrix. To guarantee convergence, it also
 * must be strictly diagonally dominant.
 */
Matrix *matrix_solve_numerical(const Matrix *a, const Matrix *b, int iters);

// Access elements of a matrix:

inline double matrix_get(const Matrix *mat, int i, int j) {
    return mat->p[mat->cols * i + j]; 
}

inline void matrix_set(Matrix *mat, int i, int j, double x) {
    mat->p[mat->cols * i + j] = x; 
}

inline double *matrix_get_ptr(const Matrix *mat, int i, int j) {
    return &mat->p[mat->cols * i + j]; 
}

#endif // MATRIX_H
