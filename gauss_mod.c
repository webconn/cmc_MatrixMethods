#include "gauss_mod.h"
#include <math.h>

static int *sequence;

static void gv_swap(int i, int j, void *arg)
{
        int t = sequence[i];
        sequence[i] = sequence[j];
        sequence[j] = t;
}

static void gv_div(int i, number_t d, void *arg)
{
        vector_t *f = (vector_t *) arg;

        f->vector[i] /= d;
}

static void gv_hit(int src, int dst, number_t mul, void *arg)
{
        vector_t *f = (vector_t *) arg;

        f->vector[dst] += mul * f->vector[src];
}

static gauss_mod_handlers_t vect_ops = {
        .swap = gv_swap,
        .div = gv_div,
        .hit = gv_hit
};

number_t gauss_mod_solve(matrix_t m, vector_t *f)
{
        if (f == NULL) {
                fprintf(stderr, "ERROR in gauss_solve(): No vector\n");
                return 0;
        }

        sequence = (int *) malloc (m.size * sizeof (int));
        for (int i=0; i<m.size; i++) {
                sequence[i] = i;
        }

        number_t d = gauss_mod_full(m, &vect_ops, f);

        vector_t old = vector_copy(*f);
        for (int i=0; i<m.size; i++) {
                f->vector[i] = old.vector[sequence[i]];
        }

        free(sequence);
        vector_free(old);
        
        return d;
}

static void gm_swap(int i, int j, void *arg)
{
        int t = sequence[i];
        sequence[i] = sequence[j];
        sequence[j] = t;
}

static void gm_div(int i, number_t div, void *arg)
{
        matrix_t *m = (matrix_t *) arg;

        matrix_divRow(*m, i, div);
}

static void gm_hit(int src, int dst, number_t mul, void *arg)
{
        matrix_t *m = (matrix_t *) arg;

        matrix_hitRow(*m, src, dst, mul);
}

static gauss_mod_handlers_t matr_ops = {
        .swap = gm_swap,
        .div = gm_div,
        .hit = gm_hit
};

matrix_t gauss_mod_invert(matrix_t m)
{
        sequence = (int *) malloc (m.size * sizeof (int));

        matrix_t result = matrix_create(m.size);
        for (int i=0; i<m.size; i++) {
                result.matrix[i][i] = 1;
                sequence[i] = i;
        }
        gauss_mod_full(m, &matr_ops, &result);

        int *old_rows = (int *) malloc (m.size * sizeof (int));
        
        for (int i=0; i<m.size; i++) {
                old_rows[i] = result.rows[i];
        }

        for (int i=0; i<m.size; i++) {
                result.rows[i] = old_rows[sequence[i]];
        }

        free(old_rows);
        free(sequence);
        return result;
}

number_t gauss_mod_full(matrix_t m, gauss_mod_handlers_t *h, void *arg)
{
        number_t det = gauss_mod_direct(m, h, arg);
        gauss_mod_reverse(m, h, arg);
        return det;
}

number_t gauss_mod_direct(matrix_t m, gauss_mod_handlers_t *h, void *arg)
{
        number_t det = 1;

        for (int i = 0; i < m.size; i++) {
                int base_elem = i;

                number_t max_elem = fabs(matrix_elem(m, i, i));
                int max_elem_ind = i;
                for (int j = i + 1; j < m.size; j++) {
                        if (fabs(matrix_elem(m, i, j)) > max_elem) {
                                max_elem = fabs(matrix_elem(m, i, j));
                                max_elem_ind = j;
                        }
                }

                if (!NOT_ZERO(max_elem)) {
                        fprintf(stderr, "ERROR: Matrix is singular\n");
                        return 0;
                }

                matrix_swapCols(m, i, max_elem_ind);
                h->swap(i, max_elem_ind, arg); /* vector_exchangeElems(f, i, base); */

                if (i != max_elem_ind)
                        det = -det;

                number_t divider = matrix_elem(m, i, i);
                det *= divider;
                h->div(i, divider, arg); /* f.vector[i] /= divider; */
                matrix_divRow(m, i, divider);

                for (int j = i + 1; j < m.size; j++) {
                        if (NOT_ZERO(matrix_elem(m, j, i))) {
                                number_t mul = matrix_elem(m, j, i);
                                matrix_hitRow(m, i, j, -mul);
                                h->hit(i, j, -mul, arg); /* f.vector[j] -= mul * f.vector[i]; */
                        }
                }
        }

        return det;
}

void gauss_mod_reverse(matrix_t m, gauss_mod_handlers_t *h, void *arg)
{
        for (int i = m.size - 1; i >= 0; i--) {
                for (int j = i - 1; j >= 0; j--) {
                        if (NOT_ZERO(matrix_elem(m, j, i))) {
                                h->hit(i, j, -matrix_elem(m, j, i), arg); /* f.vector[j] -= matrix_elem(m, j, i) * f.vector[i]; */
                                matrix_elem(m, j, i) = 0;
                        }
                }
        }
}
