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

int main(void) {
    int equations, unknowns;
    printf("Equations: ");
    scanf("%d", &equations);
    printf("Unknowns: ");
    scanf("%d", &unknowns);

    double q[equations][1];
    double p[equations][unknowns];
    printf("Matrix A|B:\n");
    for(int i = 0; i < equations; ++i) {
        for(int j = 0; j < unknowns; ++j)
            scanf("%lf", &p[i][j]);
        scanf("%lf", &q[i][0]);
    }

    Matrix a, b;
    matrix_init_with(&a, equations, unknowns, p);
    matrix_init_with(&b, equations, 1, q);

    bool has_solution;
    Matrix x = matrix_solve(&a, &b, &has_solution);
    if(!has_solution)
        printf("There is no solution to Ax = B\n");
    else {
        printf("Solution to Ax = B:\n");
        matrix_print(&x);
    }

    matrix_free(&a);
    matrix_free(&b);
    matrix_free(&x);
    return 0;
}
