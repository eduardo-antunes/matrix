/* I dedicate any and all copyright interest in this software to the
 * public domain. I make this dedication for the benefit of the public at
 * large and to the detriment of my heirs and successors. I intend this
 * dedication to be an overt act of relinquishment in perpetuity of all
 * present and future rights to this software under copyright law.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "matrix.h"

// Meta operations:

static void matrix_alloc(Matrix *mat, int rows, int cols) {
    mat->p = (double**) malloc(rows * sizeof(double*));
    for(int i = 0; i < rows; ++i)
        mat->p[i] = (double*) malloc(cols * sizeof(double));
    mat->rows = rows;
    mat->cols = cols;
}

void matrix_init(Matrix *mat, int rows, int cols) {
    matrix_alloc(mat, rows, cols);
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            mat->p[i][j] = 0;
}

void matrix_init_with(Matrix *mat, int rows, int cols, double p[rows][cols]) {
    matrix_alloc(mat, rows, cols);
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            mat->p[i][j] = p[i][j];
}

void matrix_init_id(Matrix *mat, int order) {
    matrix_alloc(mat, order, order);
    for(int i = 0; i < order; ++i)
        for(int j = 0; j < order; ++j)
            mat->p[i][j] = (i == j) ? 1 : 0;
}

Matrix matrix_copy(Matrix *mat) {
    Matrix copy;
    matrix_init(&copy, mat->rows, mat->cols);
    for(int i = 0; i < mat->rows; ++i)
        for(int j = 0; j < mat->cols; ++j)
            copy.p[i][j] = mat->p[i][j];
    return copy;
}

void matrix_free(Matrix *mat) {
    for(int i = 0; i < mat->rows; ++i)
        free(mat->p[i]);
    free(mat->p);
}

void matrix_print(const Matrix *mat) {
    for(int i = 0; i < mat->rows; ++i) {
        printf("%g", mat->p[i][0]);
        for(int j = 1; j < mat->cols; ++j)
            printf(" %g", mat->p[i][j]);
        printf("\n");
    }
}

// Basic operations:

bool matrix_square(Matrix *mat) {
    return mat->rows == mat->cols;
}

void matrix_scale(Matrix *mat, double k) {
    for(int i = 0; i < mat->rows; ++i)
        for(int j = 0; j < mat->cols; ++j)
            mat->p[i][j] *= k;
}

Matrix matrix_add(Matrix *a, Matrix *b) {
    Matrix res;
    matrix_init(&res, a->rows, a->cols);
    for(int i = 0; i < a->rows; ++i)
        for(int j = 0; j < a->cols; ++j)
            res.p[i][j] = a->p[i][j] + b->p[i][j];
    return res;
}

Matrix matrix_prod(Matrix *a, Matrix *b) {
    Matrix res;
    matrix_init(&res, a->rows, b->cols);
    for(int i = 0; i < a->rows; ++i)
        for(int j = 0; j < b->cols; ++j)
            for(int k = 0; k < a->cols; ++k)
                res.p[i][j] += a->p[i][k] * b->p[k][j];
    return res;
}

// Gaussian reduction operations:

void matrix_swap_rows(Matrix *mat, int r1, int r2) {
    double *aux = mat->p[r1];
    mat->p[r1] = mat->p[r2];
    mat->p[r2] = aux;
}

void matrix_scale_row(Matrix *mat, int r, double k) {
    for(int j = 0; j < mat->cols; ++j)
        mat->p[r][j] *= k;
}

void matrix_divide_row(Matrix *mat, int r, double k) {
    assert(k != 0);
    for(int j = 0; j < mat->cols; ++j)
        mat->p[r][j] /= k;
}

void matrix_add_rows(Matrix *mat, int rr, int rf, double k) {
    for(int j = 0; j < mat->cols; ++j)
        mat->p[rr][j] += k * mat->p[rf][j];
}

static void find_pivots(Matrix *mat, int *pivots) {
    for(int i = 0; i < mat->rows; ++i) {
        for(int j = 0; j < mat->cols; ++j) {
            if(mat->p[i][j] != 0) {
                pivots[i] = j;
                break;
            }
            pivots[i] = mat->cols;
        }
    }
}

static void sort_by_pivots(Matrix *mat, int *pivots) {
    int aux;
    for(int i = 1; i < mat->rows; ++i) {
        aux = pivots[i];
        int j = i - 1;
        while(j >= 0 && aux < pivots[j]) {
            matrix_swap_rows(mat, j, j + 1);
            pivots[j + 1] = pivots[j];
            --j;
        }
        pivots[j + 1] = aux;
    }
}

void matrix_reduce(Matrix *mat) {
    int j, pivots[mat->rows];
    find_pivots(mat, pivots);
    for(int i = 0; i < mat->rows; ++i) {
        // Step #1
        sort_by_pivots(mat, pivots);
        // Step #2
        j = pivots[i];
        if(j >= mat->cols) break;
        double t = mat->p[i][j];
        matrix_divide_row(mat, i, t);
        // Step #3
        for(int k = 0; k < mat->rows; ++k) {
            if(k == i) continue;
            t = -mat->p[k][j];
            matrix_add_rows(mat, k, i, t);
        }
        // Step #4
        find_pivots(mat, pivots);
    }
}

// System solving operations:

Matrix matrix_augment(Matrix *a, Matrix *b) {
    Matrix ab;
    matrix_init(&ab, a->rows, a->cols + b->cols);
    for(int i = 0; i < ab.rows; ++i) {
        for(int j = 0; j < ab.cols; ++j) {
            ab.p[i][j] = (j < a->cols) ?
                a->p[i][j] : b->p[i][j - a->cols];
        }
    }
    return ab;
}

Matrix matrix_solve(Matrix *a, Matrix *b, bool *has_solution) {
    *has_solution = true; // assume there is a solution
    Matrix sys = matrix_augment(a, b);
    matrix_reduce(&sys);

    Matrix x;
    matrix_init(&x, a->cols, 1);
    int h, last_col = sys.cols - 1;
    for(int i = 0; i < sys.rows; ++i) {
        h = 0;
        for(int j = 0; j < last_col; ++j)
            h += sys.p[i][j];
        if(h == 0) {
            // this indicates that rank(A) < rank(A|B)
            // thus, there is no solution to the system
            *has_solution = false;
            break;
        }
        x.p[i][0] = sys.p[i][last_col];
    }

    matrix_free(&sys);
    return x;
}

Matrix matrix_inverse(Matrix *mat) {
    assert(matrix_square(mat));
    Matrix id, inv;
    matrix_init(&inv, mat->rows, mat->cols);
    matrix_init_id(&id, mat->rows);

    Matrix aug = matrix_augment(mat, &id);
    // aug = mat|I, so by reducing it, we obtain I|inv
    matrix_reduce(&aug);
    for(int i = 0; i < mat->rows; ++i)
        for(int j = 0; j < mat->cols; ++j)
            inv.p[i][j] = aug.p[i][j + mat->cols];

    matrix_free(&aug);
    matrix_free(&id);
    return inv;
}
