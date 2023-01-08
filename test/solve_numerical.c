#include <stdio.h>

#include "matrix.h"

int main(void) {
    int equations;
    printf("Equations: ");
    scanf("%d", &equations);

    double p[equations][equations], q[equations][1];
    printf("Matrix A|B (%dx%d):\n", equations, equations + 1);
    for(int i = 0; i < equations; ++i) {
        for(int j = 0; j < equations; ++j)
            scanf("%lf", &p[i][j]);
        scanf("%lf", &q[i][0]);
    }

    Matrix a, b;
    matrix_init_with(&a, equations, equations, p);
    matrix_init_with(&b, equations, 1, q);

    if(!matrix_is_diagonally_dominant(&a))
        printf("A is not diagonally dominant. Thus, convergence is not "
                "guaranteed, and the results might be imprecise\n");

    // You can pass any number less than or equal to 0 for the iterations
    // to use the default chosen by the library.
    Matrix x = matrix_solve_numerical(&a, &b, -1);
    printf("Solution to Ax = B, obtained by Gauss-Seidel:\n");
    matrix_print(&x);

    matrix_free(&a);
    matrix_free(&b);
    matrix_free(&x);
    return 0;
}
