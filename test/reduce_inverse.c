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
