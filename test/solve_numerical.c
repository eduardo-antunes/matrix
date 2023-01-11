/* I dedicate any and all copyright interest in this software to the
 * public domain. I make this dedication for the benefit of the public at
 * large and to the detriment of my heirs and successors. I intend this
 * dedication to be an overt act of relinquishment in perpetuity of all
 * present and future rights to this software under copyright law.
 */

#include <stdio.h>

#include "matrix.h"

int main(void) {
    int equations;
    printf("Equations: ");
    scanf("%d", &equations);

    double value;
    Matrix *a = matrix_alloc(equations, equations);
    Matrix *b = matrix_alloc(equations, 1);
    printf("Matrix A|B (%dx%d):\n", equations, equations + 1);
    for(int i = 0; i < equations; ++i) {
        for(int j = 0; j < equations; ++j) {
            scanf("%lf", &value);
            matrix_set(a, i, j, value);
        }
        scanf("%lf", &value);
        matrix_set(b, i, 0, value);
    }

    if(!matrix_is_diagonally_dominant(a))
        printf("A is not diagonally dominant. Thus, convergence is not "
                "guaranteed, and the results might be imprecise\n");

    // You can pass any number less than or equal to 0 for the iterations
    // to use the default chosen by the library.
    Matrix *x = matrix_solve_numerical(a, b, -1);
    printf("Solution to Ax = B, obtained by Gauss-Seidel:\n");
    matrix_print(x);

    matrix_free(a);
    matrix_free(b);
    matrix_free(x);
    return 0;
}
