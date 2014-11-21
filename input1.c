#include "input.h"

#include <stdio.h>
#include <stdlib.h>

#define M 8
#define N 20

void input_form1(matrix_t *m, vector_t *f)
{
        input_form1_m(m);

        *f = vector_create(N);
        for (int i=0; i<N; i++) {
                f->vector[i] = 200 + 50 * (i + 1);
        }
        
        if (N <= 10) {
                fprintf(stderr, "[INPUT] Vector:\n");
                vector_print(stderr, *f, FORMAT_TEXT);
        }
}

void input_form1_m(matrix_t *m)
{
        *m = matrix_create(N);

        for (int i=0; i<N; i++) {
                for (int j=0; j<N; j++) {
                        if (i == j) { 
                                m->matrix[i][j] = (number_t) N + (number_t) M*M + (number_t) (j + 1) / M + (number_t) (i + 1) / N;
                        } else {
                                m->matrix[i][j] = (number_t) (i + j + 2) / (M + N);
                        }
                }
        }

        if (N <= 10) {
                fprintf(stderr, "[INPUT] Matrix:\n");
                matrix_print(stderr, *m, FORMAT_TEXT);
        }

        
}
