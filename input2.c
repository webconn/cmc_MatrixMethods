#include "input.h"

#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#define N 100
#define M 4

void input_form2(matrix_t *m, vector_t *f)
{
        input_form2_m(m);

        *f = vector_create(N);

        number_t x;
        fprintf(stderr, "[INPUT] Type X:\n");
        fscanf(stdin, NUMBER_READ_FORMAT, &x);

        for (int i=0; i<N; i++) {
                f->vector[i] = (number_t) N * exp(x / (i + 1)) * cos(x);
        }

        if (N <= 20) {
                fprintf(stderr, "[INPUT] Vector:\n");
                vector_print(stderr, *f, FORMAT_TEXT);
        } 
}

void input_form2_m(matrix_t *m)
{
        *m = matrix_create(N);

        number_t qm = 1.001 - 2 * M * 0.001;

        for (int i=0; i<N; i++) {
                for (int j=0; j<N; j++) {
                        if (i == j) {
                                m->matrix[i][j] = pow(qm - 1, i + j + 2);
                        } else {
                                m->matrix[i][j] = pow(qm, i + j + 2) + 0.1 * (j - i);
                        }
                }
        }
        
        if (N <= 20) {
                fprintf(stderr, "[INPUT] Matrix:\n");
                matrix_print(stderr, *m, FORMAT_TEXT);
        }
}
