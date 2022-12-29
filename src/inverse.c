/* I dedicate any and all copyright interest in this software to the
 * public domain. I make this dedication for the benefit of the public at
 * large and to the detriment of my heirs and successors. I intend this
 * dedication to be an overt act of relinquishment in perpetuity of all
 * present and future rights to this software under copyright law.
 */

#include <stdio.h>

#include "matrix.h"

int main(void) {
    int n;
    printf("Order of the matrix A: ");
    scanf("%d", &n);

    double p[n][n];
    printf("Matrix A (%dx%d):\n", n, n);
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            scanf("%lf", &p[i][j]);

    double q[n][1];
    printf("Matrix B (%dx1):\n", n);
    for(int i = 0; i < n; ++i)
        scanf("%lf", &q[i][0]);

    Matrix a, b, inv, prod;
    matrix_init_with(&a, n, n, p);
    inv = matrix_inverse(&a);
    printf("Matrix inv(A):\n");
    matrix_print(&inv);

    matrix_init_with(&b, n, 1, q);
    prod = matrix_prod(&inv, &b);
    printf("Solution of Ax = B:\n");
    matrix_print(&prod);

    matrix_free(&a);
    matrix_free(&b);
    matrix_free(&inv);
    matrix_free(&prod);
    return 0;
}
