#include "gauss_mod.h"
#include <math.h>

/* Массив, в котором будет содержаться перестановка переменных при перестановке
 * строк в модифицированном методе Гаусса */
static int *sequence;

/** Секция векторных операций (для решения СЛАУ) */

/* Перестановка переменных в векторе (будет произведена в самом конце) */
static void gv_swap(int i, int j, void *arg)
{
        int t = sequence[i];
        sequence[i] = sequence[j];
        sequence[j] = t;
}

/* Деление элемента вектора на число d */
static void gv_div(int i, number_t d, void *arg)
{
        if (!NOT_ZERO(d))
                return;

        vector_t *f = (vector_t *) arg;

        f->vector[i] /= d;
}


/* Линейная комбинация: к элементу вектора dst прибавляется элемент
 * вектора src, умноженный на коэффициент mul */
static void gv_hit(int src, int dst, number_t mul, void *arg)
{
        vector_t *f = (vector_t *) arg;

        f->vector[dst] += mul * f->vector[src];
}

/* Таблица операций для решения СЛАУ модифицированным методом Гаусса */
static gauss_mod_handlers_t vect_ops = {
        .swap = gv_swap,
        .div = gv_div,
        .hit = gv_hit
};

/* Решение СЛАУ модицифицированным методом Гаусса */
number_t gauss_mod_solve(matrix_t m, vector_t *f)
{
        if (f == NULL) {
                fprintf(stderr, "ERROR in gauss_solve(): No vector\n");
                return 0;
        }

        /* Создаём массив - перестановку переменных */
        sequence = (int *) malloc (m.size * sizeof (int));
        for (int i=0; i<m.size; i++) {
                sequence[i] = i;
        }

        /* Проводим полный цикл модифицированного метода Гаусса */
        number_t d = gauss_mod_full(m, &vect_ops, f);

        /* Копируем полученный вектор значений; в нём переменные идут в 
         * первоначальном порядке. Затем переставляем значения согласно
         * новому порядку переменных */
        vector_t old = vector_copy(*f);
        for (int i=0; i<m.size; i++) {
                f->vector[sequence[i]] = old.vector[i];
        }

        /* Освобождаем память */
        free(sequence);
        vector_free(old);
        
        return d;
}

/** Секция матричных операций (для вычисления обратной матрицы) */

/* Перестановка строк матрицы (откладывается до конца вычислений) */
static void gm_swap(int i, int j, void *arg)
{
        int t = sequence[i];
        sequence[i] = sequence[j];
        sequence[j] = t;
}

/* Деление строки матрицы на число div */
static void gm_div(int i, number_t div, void *arg)
{
        if (!NOT_ZERO(div))
                return;

        matrix_t *m = (matrix_t *) arg;

        matrix_divRow(*m, i, div);
}

/* Линейная комбинация строк матрицы: к строке dst прибавляется строка
 * src, умноженная на число mul */
static void gm_hit(int src, int dst, number_t mul, void *arg)
{
        matrix_t *m = (matrix_t *) arg;

        matrix_hitRow(*m, src, dst, mul);
}

/* Таблица операций для вычисления обратной матрицы модифицированным
 * методом Гаусса */
static gauss_mod_handlers_t matr_ops = {
        .swap = gm_swap,
        .div = gm_div,
        .hit = gm_hit
};

/* Вычисление обратной матрицы модицифированным методом Гаусса */
matrix_t gauss_mod_invert(matrix_t m)
{
        /* Создаём единичную матрицу для подготовки результата, а также
         * выделяем память для хранения перестановки столбцов исходной матрицы.
         * В итоге эта перестановка будет применяться к строкам обратной */
        sequence = (int *) malloc (m.size * sizeof (int));

        matrix_t result = matrix_create(m.size);
        for (int i=0; i<m.size; i++) {
                result.matrix[i][i] = 1;
                sequence[i] = i;
        }

        /* Выполняем полный цикл модиф. метода Гаусса: прямой и обратный ход */
        gauss_mod_full(m, &matr_ops, &result);

        /* Сохраняем предыдущий порядок следования строк обратной матрицы */
        int *old_rows = (int *) malloc (m.size * sizeof (int));
        for (int i=0; i<m.size; i++) {
                old_rows[i] = result.rows[i];
        }

        /* Переставляем строки обратной матрицы согласно перестановке столбцов
         * исходной матрицы */
        for (int i=0; i<m.size; i++) {
                result.rows[i] = old_rows[sequence[i]];
        }

        /* Освобождаем память */
        free(old_rows);
        free(sequence);

        return result;
}

/** Секция общих операций модифицированного метода Гаусса */

/* Полный цикл модицифированнного метода Гаусса: прямой и обратный ход */
number_t gauss_mod_full(matrix_t m, gauss_mod_handlers_t *h, void *arg)
{
        /* Прямой ход; по пути считаем определитель матрицы */
        number_t det = gauss_mod_direct(m, h, arg);

        /* Обратный ход */
        gauss_mod_reverse(m, h, arg);
        return det;
}

/* Прямой ход модифицированного метода Гаусса */
number_t gauss_mod_direct(matrix_t m, gauss_mod_handlers_t *h, void *arg)
{
        number_t det = 1;

        /* Проходим матрицу построчно; на каждом этапе в строке выбираем столбец с
         * наибольшим по модулю элементом и перемещаем этот столбец в текущее
         * положение. Таким образом, на диагонали окажутся наибольшие по модулю
         * элементы - главные элементы матрицы. */
        for (int i = 0; i < m.size; i++) {
                int base_elem = i;

                /* Выбираем максимальный по модулю элемент в строке */
                number_t max_elem = fabs(matrix_elem(m, i, i));
                int max_elem_ind = i;
                for (int j = i + 1; j < m.size; j++) {
                        if (fabs(matrix_elem(m, i, j)) > max_elem) {
                                max_elem = fabs(matrix_elem(m, i, j));
                                max_elem_ind = j;
                        }
                }

                /* Если наибольший по модулю элемент нулевой - матрица вырождена */
                if (!NOT_ZERO(max_elem)) {
                        fprintf(stderr, "ERROR: Matrix is singular\n");
                        return 0;
                }

                /* Меняем текущий столбец на столбец с главным элементом */
                matrix_swapCols(m, i, max_elem_ind);
                h->swap(i, max_elem_ind, arg);

                /* Если произошла смена столбцов - меняем знак определителя */
                if (i != max_elem_ind)
                        det = -det;

                /* Делим значения строки на значение главного элемента. На это же
                 * значение умножаем значение будущего определителя матрицы */
                number_t divider = matrix_elem(m, i, i);
                det *= divider;
                h->div(i, divider, arg);
                matrix_divRow(m, i, divider);

                /* Вычитаем строку из тех оставшихся, в которых в данном столбце 
                 * остались ненулевые элементы */
                for (int j = i + 1; j < m.size; j++) {
                        if (NOT_ZERO(matrix_elem(m, j, i))) {
                                number_t mul = matrix_elem(m, j, i);
                                matrix_hitRow(m, i, j, -mul);
                                h->hit(i, j, -mul, arg);
                        }
                }
        }

        return det;
}

/* Обратный ход модифицированного метода Гаусса; ничем не отличается от такого
 * для оригинального метода */
void gauss_mod_reverse(matrix_t m, gauss_mod_handlers_t *h, void *arg)
{
        for (int i = m.size - 1; i >= 0; i--) {
                for (int j = i - 1; j >= 0; j--) {
                        if (NOT_ZERO(matrix_elem(m, j, i))) {
                                h->hit(i, j, -matrix_elem(m, j, i), arg);
                                matrix_elem(m, j, i) = 0;
                        }
                }
        }
}
