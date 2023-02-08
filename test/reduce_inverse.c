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

    double value;
    Matrix *q = matrix_alloc(n, n);
    printf("Matrix Q:\n");
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j) {
            scanf("%lf", &value);
            matrix_set(q, i, j, value);
        }

    Matrix *r = matrix_reduce_const(q);
    printf("Matrix Q in reduced echelon form:\n");
    matrix_print(r);
    matrix_free(r);

    r = matrix_inverse(q);
    printf("Inverse of matrix Q:\n");
    matrix_print(r);

    matrix_free(q);
    matrix_free(r);
    return 0;
}
