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
    int equations, unknowns;
    printf("Equations: ");
    scanf("%d", &equations);
    printf("Unknowns: ");
    scanf("%d", &unknowns);

    double value;
    Matrix *a = matrix_alloc(equations, unknowns);
    Matrix *b = matrix_alloc(equations, 1);
    printf("Matrix A|B (%dx%d):\n", equations, unknowns + 1);
    for(int i = 0; i < equations; ++i) {
        for(int j = 0; j < unknowns; ++j) {
            scanf("%lf", &value);
            matrix_set(a, i, j, value);
        }
        scanf("%lf", &value);
        matrix_set(b, i, 0, value);
    }

    bool solution;
    Matrix *x = matrix_solve(a, b, &solution);
    if(solution) {
        printf("Solution to Ax = B:\n");
        matrix_print(x);
    } else {
        printf("No solution to Ax = B\n");
    }

    matrix_free(a);
    matrix_free(b);
    matrix_free(x);
    return 0;
}
