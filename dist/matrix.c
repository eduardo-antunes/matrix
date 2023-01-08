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

void matrix_init_id(Matrix *mat, int order) {
    matrix_alloc(mat, order, order);
    for(int i = 0; i < order; ++i)
        for(int j = 0; j < order; ++j)
        mat->p[i][i] = (i == j) ? 1 : 0;
}

void matrix_init_with(Matrix *mat, int rows, int cols, double p[rows][cols]) {
    matrix_alloc(mat, rows, cols);
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            mat->p[i][j] = p[i][j];
}

void matrix_init_copy(Matrix *mat, const Matrix *src) {
    matrix_alloc(mat, src->rows, src->cols);
    for(int i = 0; i < src->rows; ++i)
        for(int j = 0; j < src->cols; ++j)
            mat->p[i][j] = src->p[i][j];
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
}

// Basic arithmetic operations:

void matrix_scale(Matrix *mat, double k) {
    for(int i = 0; i < mat->rows; ++i)
        for(int j = 0; j < mat->cols; ++j)
            mat->p[i][j] *= k;
}

void matrix_add(const Matrix *a, const Matrix *b, Matrix *res) {
    for(int i = 0; i < res->rows; ++i)
        for(int j = 0; j < res->cols; ++j)
            res->p[i][j] = a->p[i][j] + b->p[i][j];
}

void matrix_prod(const Matrix *a, const Matrix *b, Matrix *res) {
    for(int i = 0; i < a->rows; ++i)
        for(int j = 0; j < b->cols; ++j)
            for(int k = 0; k < a->cols; ++k)
                res->p[i][j] += a->p[i][k] * b->p[k][j];
}

// Solving linear systems of equations:

static void swap_rows(Matrix *mat, int r1, int r2) {
    double *aux = mat->p[r1];
    mat->p[r1] = mat->p[r2];
    mat->p[r2] = aux;
}

static void find_pivots(const Matrix *mat, int *pivots) {
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

// I'd like to improve this, but I'm stupid :)
bool matrix_solve(Matrix *a, Matrix *b, Matrix *res) {
    int aux, pivots[a->rows];
    find_pivots(a, pivots);
    for(int i = 0; i < a->rows; ++i) {
        // Step #1: permute rows
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
        // Step #2: "normalize" current row
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
        // Step #3.5: recover pivot indexes
        find_pivots(a, pivots);
    }
    // Step #4: finish things off
    for(int i = 0; i < a->rows; ++i)
        if(pivots[i] >= a->rows && b->p[i] != 0)
            return false;
        else
            res->p[i][0] = b->p[i][0];
    return true;
}
