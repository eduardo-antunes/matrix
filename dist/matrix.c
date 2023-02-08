/* I dedicate any and all copyright interest in this software to the
 * public domain. I make this dedication for the benefit of the public at
 * large and to the detriment of my heirs and successors. I intend this
 * dedication to be an overt act of relinquishment in perpetuity of all
 * present and future rights to this software under copyright law.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "matrix.h"

// Access elements of a matrix: (defined in matrix.h)

double matrix_get(const Matrix*, int, int);

double *matrix_get_ptr(const Matrix*, int, int);

void matrix_set(Matrix*, int, int, double);

// Initializers, free and print:

Matrix *matrix_alloc(int rows, int cols) {
    // Allocate the matrix structure
    Matrix *m = (Matrix *) malloc(sizeof(Matrix));
    if(m == NULL) return NULL;
    // Allocate the actual data of the matrix, which will be 
    // stored as a 1D-array in row-major order
    m->p = (double *) malloc(rows * cols * sizeof(Matrix));
    if(m->p == NULL) {
        free(m);
        return NULL;
    }
    // Fill in the dimensions of the matrix
    m->rows = rows;
    m->cols = cols;
    return m;
}

Matrix *matrix_calloc(int rows, int cols) {
    Matrix *m = matrix_alloc(rows, cols);
    if(m == NULL) return NULL;
    // Zero-initialize the data
    memset(m->p, 0, rows * cols * sizeof(double));
    return m;
}

Matrix *matrix_identity(int order) {
    Matrix *m = matrix_alloc(order, order);
    if(m == NULL) return NULL;
    for(int i = 0; i < order; ++i)
        for(int j = 0; j < order; ++j)
            // Data is initialized with the Kronecker delta
            m->p[m->cols * i + j] = (i == j) ? 1 : 0;
    return m;
}

Matrix *matrix_copy(const Matrix *src) {
    Matrix *c = matrix_alloc(src->rows, src->cols);
    if(c == NULL) return NULL;
    // Copy the data from the source matrix
    memcpy(c->p, src->p, src->rows * src->cols * sizeof(double));
    return c;
}

void matrix_print(const Matrix *m) {
    int last_col = m->cols - 1;
    for(int i = 0; i < m->rows; ++i) {
        for(int j = 0; j < last_col; ++j)
            printf("%g ", m->p[m->cols * i + j]);
        printf("%g\n", m->p[m->cols * i + last_col]);
    }
}

void matrix_free(Matrix *m) {
    free(m->p);
    free(m);
}

// Arithmetic operations:

void matrix_scale(Matrix *m, double k) {
    for(int i = 0; i < m->rows; ++i)
        for(int j = 0; j < m->cols; ++j)
            m->p[m->cols * i + j] *= k;
}

Matrix *matrix_scale_copy(const Matrix *m, double k) {
    Matrix *c = matrix_copy(m);
    for(int i = 0; i < c->rows; ++i)
        for(int j = 0; j < c->cols; ++j)
            c->p[c->cols * i + j] *= k;
    return c;
}

Matrix *matrix_add(const Matrix *a, const Matrix *b) {
    Matrix *res = matrix_alloc(a->rows, a->cols);
    for(int i = 0; i < res->rows; ++i)
        for(int j = 0; j < res->cols; ++j)
            res->p[res->cols * i + j] = a->p[a->cols * i + j] 
                + b->p[b->cols * i + j];
    return res;
}

Matrix *matrix_prod(const Matrix *a, const Matrix *b) {
    Matrix *res = matrix_calloc(a->rows, b->cols);
    for(int i = 0; i < a->rows; ++i)
        for(int j = 0; j < b->cols; ++j)
            for(int k = 0; k < a->cols; ++k)
                res->p[res->cols * i + j] += a->p[a->cols * i + k] 
                    * b->p[b->cols * k + j];
    return res;
}

Matrix *matrix_transpose(const Matrix *m) {
    Matrix *tr = matrix_alloc(m->cols, m->rows);
    for(int i = 0; i < tr->rows; ++i)
        for(int j = 0; j < tr->cols; ++j)
            tr->p[tr->cols * i + j] = m->p[m->cols * j + i];
    return tr;
}

// Tests and conditions:

bool matrix_are_equal(const Matrix *a, const Matrix *b) {
    if(a->rows != b->rows || a->cols != b->cols) 
        return false;
    for(int i = 0; i < a->rows; ++i)
        for(int j = 0; j < a->cols; ++j)
            if(a->p[a->cols * i + j] != b->p[b->cols * i + j])
                return false;
    return true;
}

bool matrix_is_square(const Matrix *mat) {
    return mat->rows == mat->cols;
}

bool matrix_is_diagonally_dominant(const Matrix *m) {
    for(int i = 0; i < m->rows; ++i) {
        double s = 0;
        for(int j = 0; j < m->cols; ++j) {
            if(i == j) continue;
            s += m->p[m->cols * i + j];
        }
        if(s >= m->p[m->cols * i + i]) 
            return false;
    }
    return true;
}

// Solving linear systems of equations:

static void swap_rows(Matrix *m, int r1, int r2) {
    for(int j = 0; j < m->cols; ++j) {
        double aux = m->p[m->cols * r1 + j];
        m->p[m->cols * r1 + j] = m->p[m->cols * r2 + j];
        m->p[m->cols * r2 + j] = aux;
    }
}

static void find_pivots(const Matrix *m, int *pivots) {
    for(int i = 0; i < m->rows; ++i)
        for(int j = 0; j < m->cols; ++j) {
            if(m->p[m->cols * i + j] != 0) {
                pivots[i] = j;
                break;
            }
            pivots[i] = m->cols;
        }
}

void matrix_reduce(Matrix *m) {
    int aux, pivots[m->rows];
    for(int i = 0; i < m->rows; ++i) {
        // Step #1: sort rows by insertion, with pivot indexes as keys
        find_pivots(m, pivots);
        for(int t = 1; t < m->rows; ++t) {
            aux = pivots[t];
            int u = t - 1;
            while(u >= 0 && aux < pivots[u]) {
                swap_rows(m, u, u + 1);
                pivots[u + 1] = pivots[u];
                --u;
            }
            pivots[u + 1] = aux;
        }
        // Step #2: divide the i-th row by the value of its pivot
        int p = pivots[i];
        if(p >= m->cols) break;
        double c = m->p[m->cols * i + p];
        for(int j = 0; j < m->cols; ++j)
            m->p[m->cols * i + j] /= c;
        // Step #3: gaussian elimination
        for(int k = 0; k < m->rows; ++k) {
            if(k == i) continue;
            c = m->p[m->cols * k + p];
            for(int j = 0; j < m->cols; ++j)
                m->p[m->cols * k + j] -= m->p[m->cols * i + j] * c;
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
            ab->p[ab->cols * i + j] = (j < a->cols) ?
                a->p[a->cols * i + j] : b->p[b->cols * i + j - a->cols];
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
            inv->p[inv->cols * i + j] = aug->p[aug->cols * i + j + mat->cols];

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
        double c = a->p[a->cols * i + p];
        for(int j = 0; j < a->cols; ++j)
            a->p[a->cols * i + j] /= c;
        b->p[b->cols * i] /= c;
        // Step #3: gaussian elimination
        for(int k = 0; k < a->rows; ++k) {
            if(k == i) continue;
            c = a->p[a->cols * k + p];
            for(int j = 0; j < a->cols; ++j)
                a->p[a->cols * k + j] -= a->p[a->cols * i + j] * c;
            b->p[b->cols * k] -= b->p[b->cols * i] * c;
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
            res->p[res->cols * i] = b->p[b->cols * i];
    return res;
}

Matrix *matrix_solve_copy(const Matrix *a, const Matrix *b, bool *solution) {
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
            res->p[res->cols * i] = b->p[b->cols * i];
            for(int j = 0; j < res->rows; ++j) {
                if(i == j) continue;
                res->p[res->cols * i] -= a->p[a->cols * i + j] 
                    * res->p[res->cols * j];
            }
            res->p[res->cols * i] /= a->p[a->cols * i + i];
        }
    }
    return res;
}

