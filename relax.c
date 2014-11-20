#include "relax.h"

#include "matrix.h"
#include <math.h>

static void relax_iteration(matrix_t m, vector_t f, vector_t old, vector_t new, number_t omega)
{
        if (omega <= 0 || omega >= 2)
                return; /* such values are ineffective for solving */

        for (int i=0; i<old.size; i++) {
                number_t sum = f.vector[i];
                
                for (int j=0; j < i; j++) {
                        sum -= m.matrix[i][j] * new.vector[j];
                }

                for (int j=i+1; j<old.size; j++) {
                        sum -= m.matrix[i][j] * old.vector[j];
                }
                
                if (!NOT_ZERO(m.matrix[i][i])) {
                        fprintf(stderr, "ERROR: Zero on diagonal\n");
                        return;
                }

                sum /= m.matrix[i][i];

                new.vector[i] = sum;
        }

        for (int i=0; i<new.size; i++) {
                new.vector[i] = new.vector[i] * omega + old.vector[i] * (1 - omega);
        }
}

void relax_solve(matrix_t m, vector_t *f, number_t omega, number_t eps)
{
        vector_t x1, x2, diff, tmp_vect;
        x1 = vector_create(m.size);
        x2 = vector_create(m.size);
        diff = vector_create(m.size);
        
        number_t norm = 0;
        int iters = 0;
        do {
                relax_iteration(m, *f, x1, x2, omega);
                matrix_mulMatVector(m, x2, x1); /* now in x1 we have A*x(k+1) to compare with resulting f vector */
                vector_sub(x1, *f, diff);
                norm = vector_norm(diff);

                tmp_vect = x1;
                x1 = x2;
                x2 = tmp_vect;
        } while (iters++ < 50 && fabs(norm) >= eps);

        for (int i=0; i<f->size; i++) {
                f->vector[i] = x1.vector[i];
        }

        vector_free(x1);
        vector_free(x2);
        vector_free(diff);
}
