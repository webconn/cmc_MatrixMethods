#include "gauss.h"
#include "matrix.h"

/** Секция векторных операций для метода Гаусса */
static void gv_swap(int i, int j, void *arg)
{
        vector_t *f = (vector_t *) arg;
        
        number_t t = f->vector[i];
        f->vector[i] = f->vector[j];
        f->vector[j] = t;
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

static gauss_handlers_t vect_ops = {
        .swap = gv_swap,
        .div = gv_div,
        .hit = gv_hit
};

number_t gauss_solve(matrix_t m, vector_t *f)
{
        if (f == NULL) {
                fprintf(stderr, "ERROR in gauss_solve(): No vector\n");
                return 0;
        }

        number_t d = gauss_full(m, &vect_ops, f);
        
        return d;
}

/** Секция матричных операций для метода Гаусса (обратная матрица) */
static void gm_swap(int i, int j, void *arg)
{
        matrix_t *m = (matrix_t *) arg;

        matrix_exchangeRows(*m, i, j);
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

static gauss_handlers_t matr_ops = {
        .swap = gm_swap,
        .div = gm_div,
        .hit = gm_hit
};

matrix_t gauss_invert(matrix_t m)
{
        matrix_t result = matrix_create(m.size);
        for (int i=0; i<m.size; i++) {
                result.matrix[i][i] = 1;
        }
        gauss_full(m, &matr_ops, &result);
        return result;
}

/** Секция общих операций метода Гаусса */

number_t gauss_full(matrix_t m, gauss_handlers_t *h, void *arg)
{
        number_t det = gauss_direct(m, h, arg);
        gauss_reverse(m, h, arg);
        return det;
}

number_t gauss_direct(matrix_t m, gauss_handlers_t *h, void *arg)
{
        number_t det = 1;

        /* Пройдём всю матрицу построчно, на каждом этапе выбирая первую строку, i-й элемент которой
         * не равен 0 */
        for (int i = 0; i < m.size; i++) {
                int base = i;

                while (!NOT_ZERO(matrix_elem(m, base, i)) && base < m.size)
                        base++;

                if (base == m.size) {
                        fprintf(stderr, "ERROR: Matrix is singular\n");
                        return 0;
                }

                matrix_exchangeRows(m, i, base);
                h->swap(i, base, arg); /* vector_exchangeElems(f, i, base); */

                if (i != base)
                        det = -det;

                /* Теперь работаем с i-й строкой матрицы - здесь диагональный элемент ненулевой */
                /* Делим значения строки на значение первого элемента. На это же значение умножаем
                 * значение будущего определителя матрицы */
                number_t divider = matrix_elem(m, i, i);
                det *= divider;
                h->div(i, divider, arg); /* f.vector[i] /= divider; */
                matrix_divRow(m, i, divider);

                /* Вычитаем строку из всех тех оставшихся, в которых i-й элемент ненулевой
                 * (с домножением на соответствующий коэффициент */
                for (int j = i + 1; j < m.size; j++) {
                        if (NOT_ZERO(matrix_elem(m, j, i))) {
                                number_t mul = matrix_elem(m, j, i);
                                matrix_hitRow(m, i, j, -mul);
                                h->hit(i, j, -mul, arg); /* f.vector[j] -= mul * f.vector[i]; */
                        }
                }

                /*
                matrix_print(stderr, m);
                fputc('\n', stderr);
                vector_print(stderr, f);
                fputc('\n', stderr);
                fputc('\n', stderr);
                */
                
        }

        return det;
}

void gauss_reverse(matrix_t m, gauss_handlers_t *h, void *arg)
{
        /* Обратный ход метода Гаусса */
        for (int i = m.size - 1; i >= 0; i--) {
                for (int j = i - 1; j >= 0; j--) {
                        if (NOT_ZERO(matrix_elem(m, j, i))) {
                                h->hit(i, j, -matrix_elem(m, j, i), arg); /* f.vector[j] -= matrix_elem(m, j, i) * f.vector[i]; */
                                matrix_elem(m, j, i) = 0;
                        }
                }
        }
}
