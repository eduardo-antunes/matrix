/* I dedicate any and all copyright interest in this software to the
 * public domain. I make this dedication for the benefit of the public at
 * large and to the detriment of my heirs and successors. I intend this
 * dedication to be an overt act of relinquishment in perpetuity of all
 * present and future rights to this software under copyright law.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "matrix.h"

// Initializers, free and print:

Matrix *matrix_alloc(int rows, int cols) {
    Matrix *mat = (Matrix*) malloc(sizeof(Matrix));
    if(mat == NULL) return NULL;

    mat->p = (double**) malloc(rows * sizeof(double*));
    if(mat->p == NULL) return NULL;

    for(int i = 0; i < rows; ++i) {
        mat->p[i] = (double*) malloc(cols * sizeof(double));
        if(mat->p[i] == NULL) return NULL;
    }
    mat->rows = rows;
    mat->cols = cols;
    return mat;
}

Matrix *matrix_calloc(int rows, int cols) {
    Matrix *mat = matrix_alloc(rows, cols);
    if(mat == NULL) return NULL;
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            mat->p[i][j] = 0;
    return mat;
}

Matrix *matrix_identity(int order) {
    Matrix *mat = matrix_alloc(order, order);
    if(mat == NULL) return NULL;
    for(int i = 0; i < order; ++i)
        for(int j = 0; j < order; ++j)
            mat->p[i][j] = (i == j) ? 1 : 0;
    return mat;
}

Matrix *matrix_copy(const Matrix *src) {
    Matrix *mat = matrix_alloc(src->rows, src->cols);
    if(mat == NULL) return NULL;
    for(int i = 0; i < src->rows; ++i)
        for(int j = 0; j < src->cols; ++j)
            mat->p[i][j] = src->p[i][j];
    return mat;
}

void matrix_print(const Matrix *mat) {
    int last_col = mat->cols - 1;
    for(int i = 0; i < mat->rows; ++i) {
        for(int j = 0; j < last_col; ++j)
            printf("%g ", mat->p[i][j]);
        printf("%g\n", mat->p[i][last_col]);
    }
}

void matrix_free(Matrix *mat) {
    for(int i = 0; i < mat->rows; ++i)
        free(mat->p[i]);
    free(mat->p);
    free(mat);
}

// Arithmetic operations:

void matrix_scale(Matrix *mat, double k) {
    for(int i = 0; i < mat->rows; ++i)
        for(int j = 0; j < mat->cols; ++j)
            mat->p[i][j] *= k;
}

Matrix *matrix_scale_const(const Matrix *mat, double k) {
    Matrix *res = matrix_alloc(mat->rows, mat->cols);
    for(int i = 0; i < mat->rows; ++i)
        for(int j = 0; j < mat->cols; ++j)
            res->p[i][j] *= k;
    return res;
}

Matrix *matrix_add(const Matrix *a, const Matrix *b) {
    Matrix *res = matrix_alloc(a->rows, a->cols);
    for(int i = 0; i < res->rows; ++i)
        for(int j = 0; j < res->cols; ++j)
            res->p[i][j] = a->p[i][j] + b->p[i][j];
    return res;
}

Matrix *matrix_prod(const Matrix *a, const Matrix *b) {
    Matrix *res = matrix_calloc(a->rows, b->cols);
    for(int i = 0; i < a->rows; ++i)
        for(int j = 0; j < b->cols; ++j)
            for(int k = 0; k < a->cols; ++k)
                res->p[i][j] += a->p[i][k] * b->p[k][j];
    return res;
}

Matrix *matrix_transpose(const Matrix *mat) {
    Matrix *tr = matrix_alloc(mat->cols, mat->rows);
    for(int i = 0; i < tr->rows; ++i)
        for(int j = 0; j < tr->cols; ++j)
            tr->p[i][j] = mat->p[j][i];
    return tr;
}

// Tests and conditions:

bool matrix_are_equal(const Matrix *a, const Matrix *b) {
    if(a->rows != b->rows || a->cols != b->cols) 
        return false;
    for(int i = 0; i < a->rows; ++i)
        for(int j = 0; j < a->cols; ++j)
            if(a->p[i][j] != b->p[i][j])
                return false;
    return true;
}

bool matrix_is_square(const Matrix *mat) {
    return mat->rows == mat->cols;
}

bool matrix_is_diagonally_dominant(const Matrix *mat) {
    for(int i = 0; i < mat->rows; ++i) {
        double s = 0;
        for(int j = 0; j < mat->cols; ++j) {
            if(i == j) continue;
            s += mat->p[i][j];
        }
        if(s >= mat->p[i][i]) return false;
    }
    return true;
}

// Solving linear systems of equations:

static void swap_rows(Matrix *mat, int r1, int r2) {
    double *aux = mat->p[r1];
    mat->p[r1] = mat->p[r2];
    mat->p[r2] = aux;
}

static void find_pivots(const Matrix *mat, int *pivots) {
    for(int i = 0; i < mat->rows; ++i)
        for(int j = 0; j < mat->cols; ++j) {
            if(mat->p[i][j] != 0) {
                pivots[i] = j;
                break;
            }
            pivots[i] = mat->cols;
        }
}

void matrix_reduce(Matrix *mat) {
    int aux, pivots[mat->rows];
    for(int i = 0; i < mat->rows; ++i) {
        // Step #1: sort rows by insertion, with pivot indexes as keys
        find_pivots(mat, pivots);
        for(int t = 1; t < mat->rows; ++t) {
            aux = pivots[t];
            int u = t - 1;
            while(u >= 0 && aux < pivots[u]) {
                swap_rows(mat, u, u + 1);
                pivots[u + 1] = pivots[u];
                --u;
            }
            pivots[u + 1] = aux;
        }
        // Step #2: divide the i-th row by the value of its pivot
        int p = pivots[i];
        if(p >= mat->cols) break;
        double c = mat->p[i][p];
        for(int j = 0; j < mat->cols; ++j)
            mat->p[i][j] /= c;
        // Step #3: gaussian elimination
        for(int k = 0; k < mat->rows; ++k) {
            if(k == i) continue;
            c = mat->p[k][p];
            for(int j = 0; j < mat->cols; ++j)
                mat->p[k][j] -= mat->p[i][j] * c;
        }
    }
}

Matrix *matrix_reduce_const(const Matrix *mat) {
    Matrix *res = matrix_copy(mat);
    matrix_reduce(res);
    return res;
}

Matrix *matrix_augment(const Matrix *a, const Matrix *b) {
    Matrix *ab = matrix_alloc(a->rows, a->cols + b->cols);
    for(int i = 0; i < ab->rows; ++i)
        for(int j = 0; j < ab->cols; ++j)
            ab->p[i][j] = (j < a->cols) ?
                a->p[i][j] : b->p[i][j - a->cols];
    return ab;
}

Matrix *matrix_inverse(const Matrix *mat) {
    // Textbook method of inversion
    Matrix *inv = matrix_alloc(mat->rows, mat->cols);
    Matrix *id = matrix_identity(mat->rows);

    Matrix *aug = matrix_augment(mat, id);
    // reduce(A|I) = I|inv(A)
    matrix_reduce(aug);
    for(int i = 0; i < mat->rows; ++i)
        for(int j = 0; j < mat->cols; ++j)
            inv->p[i][j] = aug->p[i][j + mat->cols];

    matrix_free(aug);
    matrix_free(id);
    return inv;
}

// Oh my god some code duplication ;-;

Matrix *matrix_solve(Matrix *a, Matrix *b, bool *solution) {
    if(solution != NULL) *solution = true;
    int aux, pivots[a->rows];
    for(int i = 0; i < a->rows; ++i) {
        // Step #1: sort rows by insertion, with pivot indexes as keys
        find_pivots(a, pivots);
        for(int t = 1; t < a->rows; ++t) {
            aux = pivots[t];
            int u = t - 1;
            while(u >= 0 && aux < pivots[u]) {
                swap_rows(a, u, u + 1);
                swap_rows(b, u, u + 1);
                pivots[u + 1] = pivots[u];
                --u;
            }
            pivots[u + 1] = aux;
        }
        // Step #2: divide the i-th row by the value of its pivot
        int p = pivots[i];
        if(p >= a->cols) break;
        double c = a->p[i][p];
        for(int j = 0; j < a->cols; ++j)
            a->p[i][j] /= c;
        b->p[i][0] /= c;
        // Step #3: gaussian elimination
        for(int k = 0; k < a->rows; ++k) {
            if(k == i) continue;
            c = a->p[k][p];
            for(int j = 0; j < a->cols; ++j)
                a->p[k][j] -= a->p[i][j] * c;
            b->p[k][0] -= b->p[i][0] * c;
        }
    }
    // Step #4: finish things off
    Matrix *res = matrix_alloc(a->rows, 1);
    for(int i = 0; i < a->rows; ++i)
        if(pivots[i] >= a->rows && b->p[i] != 0) {
            // rank(A) < rank(A|B) => no solution
            if(solution != NULL) *solution = false;
            return res;
        }
        else
            res->p[i][0] = b->p[i][0];
    return res;
}

Matrix *matrix_solve_const(const Matrix *a, const Matrix *b, bool *solution) {
    Matrix *_a = matrix_copy(a);
    Matrix *_b = matrix_copy(b);

    Matrix *res = matrix_solve(_a, _b, solution);
    matrix_free(_a);
    matrix_free(_b);
    return res;
}

Matrix *matrix_solve_numerical(const Matrix *a, const Matrix *b, int iters) {
    Matrix *res = matrix_calloc(a->rows, 1);
    if(iters <= 0) iters = DEFAULT_GAUSS_SEIDEL_ITERS;

    for(int k = 0; k < iters; ++k) {
        for(int i = 0; i < res->rows; ++i) {
            res->p[i][0] = b->p[i][0];
            for(int j = 0; j < res->rows; ++j) {
                if(i == j) continue;
                res->p[i][0] -= a->p[i][j] * res->p[j][0];
            }
            res->p[i][0] /= a->p[i][i];
        }
    }
    return res;
}

// Access elements of a matrix: (external declarations)

extern double matrix_get(const Matrix*, const int, const int);

extern double *matrix_get_ptr(const Matrix*, const int, const int);

extern void matrix_set(Matrix*, const int, const int, double);
