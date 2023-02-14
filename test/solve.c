/*
 * Copyright 2022-2023 Eduardo Antunes dos Santos Vieira
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *  
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License. 
 */

#include <stdio.h>
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

    int solution;
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
