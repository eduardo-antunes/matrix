/* I dedicate any and all copyright interest in this software to the
 * public domain. I make this dedication for the benefit of the public at
 * large and to the detriment of my heirs and successors. I intend this
 * dedication to be an overt act of relinquishment in perpetuity of all
 * present and future rights to this software under copyright law.
 */

#include <stdio.h>
#include <stdbool.h>

#include "matrix.h"

int main(void) {
    double p[2][2], q[2][1];
    printf("Matrix A|B:\n");
    for(int i = 0; i < 2; ++i) {
        for(int j = 0; j < 2; ++j)
            scanf("%lf", &p[i][j]);
        scanf("%lf", &q[i][0]);
    }

    Matrix a, b, x;
    matrix_init_with(&a, 2, 2, p);
    matrix_init_with(&b, 2, 1, q);
    matrix_init(&x, 2, 1);

    bool solution = matrix_solve(&a, &b, &x);
    if(solution) {
        printf("Solution to Ax = B:\n");
        matrix_print(&x);
    } else {
        printf("No solution to Ax = B\n");
    }

    matrix_free(&a);
    matrix_free(&b);
    matrix_free(&x);
    return 0;
}
