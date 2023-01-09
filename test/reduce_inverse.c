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
    printf("Order of the matrix Q: ");
    scanf("%d", &n);

    double p[n][n];
    printf("Matrix Q:\n");
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            scanf("%lf", &p[i][j]);

    Matrix q, _q;
    matrix_init_with(&q, n, n, p);
    matrix_init_copy(&_q, &q);

    matrix_reduce(&_q);
    printf("Matrix Q in row reduced echelon form:\n");
    matrix_print(&_q);
    matrix_free(&_q);

    _q = matrix_inverse(&q);
    printf("Inverse of matrix Q:\n");
    matrix_print(&_q);

    matrix_free(&q);
    matrix_free(&_q);
    return 0;
}
