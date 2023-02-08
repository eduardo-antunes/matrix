/* I dedicate any and all copyright interest in this software to the
 * public domain. I make this dedication for the benefit of the public at
 * large and to the detriment of my heirs and successors. I intend this
 * dedication to be an overt act of relinquishment in perpetuity of all
 * present and future rights to this software under copyright law.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.h"

// L-value access to matrix elements; this abstracts the
// specific storage strategy used for the contents of the
// matrix and simplifies the implementation
#define MSUB(m, i, j) ((m)->p[(((m)->cols * i) + (j))])

// Access elements of a matrix: (defined in matrix.h)

double matrix_get(const Matrix*, int, int);

double *matrix_get_ptr(const Matrix*, int, int);

void matrix_set(Matrix*, int, int, double);

// Initializers, free and print:

Matrix *matrix_alloc(int rows, int cols) {
    // Allocate the matrix structure
    Matrix *mat = (Matrix *) malloc(sizeof(Matrix));
    if(mat == NULL) return NULL;
    // Allocate the actual data of the matrix, which will be 
    // stored as a 1D-array in row-major order
    mat->p = (double *) malloc(rows * cols * sizeof(Matrix));
    if(mat->p == NULL) {
        free(mat);
        return NULL;
    }
    // Fill in the dimensions of the matrix
    mat->rows = rows;
    mat->cols = cols;
    return mat;
}

Matrix *matrix_calloc(int rows, int cols) {
    Matrix *mat = matrix_alloc(rows, cols);
    if(mat == NULL) return NULL;
    // Zero-initialize the data
    memset(mat->p, 0, rows * cols * sizeof(double));
    return mat;
}

Matrix *matrix_identity(int order) {
    Matrix *mat = matrix_alloc(order, order);
    if(mat == NULL) return NULL;
    for(int i = 0; i < order; ++i)
        for(int j = 0; j < order; ++j)
            // Data is initialized with the Kronecker delta
            MSUB(mat, i, j) = (i == j) ? 1 : 0;
    return mat;
}

Matrix *matrix_copy(const Matrix *src) {
    Matrix *c = matrix_alloc(src->rows, src->cols);
    if(c == NULL) return NULL;
    // Copy the data from the source matrix
    memcpy(c->p, src->p, src->rows * src->cols * sizeof(double));
    return c;
}

void matrix_print(const Matrix *mat) {
    int last_col = mat->cols - 1;
    for(int i = 0; i < mat->rows; ++i) {
        for(int j = 0; j < last_col; ++j)
            printf("%g ", MSUB(mat, i, j));
        printf("%g\n", MSUB(mat, i, last_col));
    }
}

void matrix_free(Matrix *mat) {
    free(mat->p);
    free(mat);
}

// Arithmetic operations:

void matrix_scale(Matrix *mat, double k) {
    for(int i = 0; i < mat->rows; ++i)
        for(int j = 0; j < mat->cols; ++j)
            MSUB(mat, i, j) *= k;
}

Matrix *matrix_scale_copy(const Matrix *mat, double k) {
    Matrix *c = matrix_copy(mat);
    for(int i = 0; i < c->rows; ++i)
        for(int j = 0; j < c->cols; ++j)
            MSUB(c, i, j) *= k;
    return c;
}

Matrix *matrix_add(const Matrix *a, const Matrix *b) {
    Matrix *res = matrix_alloc(a->rows, a->cols);
    for(int i = 0; i < res->rows; ++i)
        for(int j = 0; j < res->cols; ++j)
            MSUB(res, i, j) = MSUB(a, i, j) + MSUB(b, i, j);
    return res;
}

Matrix *matrix_prod(const Matrix *a, const Matrix *b) {
    Matrix *res = matrix_calloc(a->rows, b->cols);
    for(int i = 0; i < a->rows; ++i)
        for(int j = 0; j < b->cols; ++j)
            for(int k = 0; k < a->cols; ++k)
                MSUB(res, i, j) += MSUB(a, i, k) * MSUB(b, k, j);
    return res;
}

Matrix *matrix_transpose(const Matrix *mat) {
    Matrix *tr = matrix_alloc(mat->cols, mat->rows);
    for(int i = 0; i < tr->rows; ++i)
        for(int j = 0; j < tr->cols; ++j)
            MSUB(tr, i, j) = MSUB(mat, j, i);
    return tr;
}

// Tests and conditions:

int matrix_are_equal(const Matrix *a, const Matrix *b) {
    if(a->rows != b->rows || a->cols != b->cols) 
        return 0;
    for(int i = 0; i < a->rows; ++i)
        for(int j = 0; j < a->cols; ++j)
            if(MSUB(a, i, j) != MSUB(b, i, j))
                return 0;
    return 1;
}

int matrix_is_square(const Matrix *mat) {
    return mat->rows == mat->cols;
}

int matrix_is_diagonally_dominant(const Matrix *mat) {
    for(int i = 0; i < mat->rows; ++i) {
        double s = 0;
        for(int j = 0; j < mat->cols; ++j) {
            if(i == j) continue;
            s += MSUB(mat, i, j);
        }
        if(s >= MSUB(mat, i, i)) 
            return 0;
    }
    return 1;
}

// Solving linear systems of equations:

static void swap_rows(Matrix *mat, int r1, int r2) {
    for(int j = 0; j < mat->cols; ++j) {
        double aux = MSUB(mat, r1, j);
        MSUB(mat, r1, j) = MSUB(mat, r2, j);
        MSUB(mat, r2, j) = aux;
    }
}

static void find_pivots(const Matrix *mat, int *pivots) {
    for(int i = 0; i < mat->rows; ++i)
        for(int j = 0; j < mat->cols; ++j) {
            if(MSUB(mat, i, j) != 0) {
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
        double c = MSUB(mat, i, p);
        for(int j = 0; j < mat->cols; ++j)
            MSUB(mat, i, j) /= c;
        // Step #3: gaussian elimination
        for(int k = 0; k < mat->rows; ++k) {
            if(k == i) continue;
            c = MSUB(mat, k, p);
            for(int j = 0; j < mat->cols; ++j)
                MSUB(mat, k, j) -= MSUB(mat, i, j) * c;
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
            MSUB(ab, i, j) = (j < a->cols) ?
                MSUB(a, i, j) : MSUB(b, i, j - a->cols);
    return ab;
}

// Textbook method of inversion
Matrix *matrix_inverse(const Matrix *mat) {
    Matrix *inv = matrix_alloc(mat->rows, mat->cols);
    Matrix *id = matrix_identity(mat->rows);
    Matrix *aug = matrix_augment(mat, id);
    // reduce(A|I) = I|inv(A)
    matrix_reduce(aug);
    for(int i = 0; i < mat->rows; ++i)
        for(int j = 0; j < mat->cols; ++j)
            MSUB(inv, i, j) = MSUB(aug, i, j + mat->cols);

    matrix_free(aug);
    matrix_free(id);
    return inv;
}

// Oh my god some code duplication ;-;

Matrix *matrix_solve(Matrix *a, Matrix *b, int *solution) {
    if(solution != NULL) *solution = 1;
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
        double c = MSUB(a, i, p);
        for(int j = 0; j < a->cols; ++j)
            MSUB(a, i, j) /= c;
        b->p[b->cols * i] /= c;
        // Step #3: gaussian elimination
        for(int k = 0; k < a->rows; ++k) {
            if(k == i) continue;
            c = MSUB(a, k, p);
            for(int j = 0; j < a->cols; ++j)
                MSUB(a, k, j) -= MSUB(a, i, j) * c;
            MSUB(b, k, 0) -= MSUB(b, i, 0) * c;
        }
    }
    // Step #4: finish things off
    Matrix *res = matrix_alloc(a->rows, 1);
    for(int i = 0; i < a->rows; ++i)
        if(pivots[i] >= a->rows && b->p[i] != 0) {
            // rank(A) < rank(A|B) => no solution
            if(solution != NULL) *solution = 0;
            return res;
        }
        else
            MSUB(res, i, 0) = MSUB(b, i, 0);
    return res;
}

Matrix *matrix_solve_copy(const Matrix *a, const Matrix *b, int *solution) {
    // Make copies first
    Matrix *_a = matrix_copy(a);
    Matrix *_b = matrix_copy(b);
    // Solve passing the copies in
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
            MSUB(res, i, 0) = MSUB(b, i, 0);
            for(int j = 0; j < res->rows; ++j) {
                if(i == j) continue;
                MSUB(res, i, 0) -= MSUB(a, i, j) * MSUB(res, j, 0);
            }
            MSUB(res, i, 0) /= MSUB(a, i, i);
        }
    }
    return res;
}

#undef MSUB
