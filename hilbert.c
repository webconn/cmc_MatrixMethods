#include "input.h"

void input_hilbert(matrix_t *m, vector_t *f)
{
        input_hilbert_m(m);
        int N = m->size;

        *f = vector_create(N);
        for (int i=0; i<N; i++) {
                f->vector[i] = i + 1;
        }
        
        if (N <= 10) {
                fprintf(stderr, "[INPUT] Vector:\n");
                vector_print(stderr, *f, FORMAT_TEXT);
        }
}

void input_hilbert_m(matrix_t *m)
{
        int N;
        fprintf(stderr, "[INPUT] Type N\n");
        scanf("%d", &N);
        
        *m = matrix_create(N);

        for (int i=0; i<N; i++) {
                for (int j=0; j<N; j++) {
                        m->matrix[i][j] = (number_t) 1 / ((number_t) i + j + 1);
                }
        }

        if (N <= 10) {
                fprintf(stderr, "[INPUT] Matrix:\n");
                matrix_print(stderr, *m, FORMAT_TEXT);
        }
}
