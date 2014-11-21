#include "relax.h"

#include "matrix.h"
#include <math.h>

/* Итерация метода верхней релаксации */
static void relax_iteration(matrix_t m, vector_t f, vector_t old, 
                                vector_t new, number_t omega)
{
        /* Ниже реализовано вычисление значения по формуле */
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

        /* Подводим параметр w */
        for (int i=0; i<new.size; i++) {
                new.vector[i] = new.vector[i] * omega + old.vector[i] * (1 - omega);
        }
}

/* Решение СЛАУ методом верхней релаксации с заданной точностью и параметром */
void relax_solve(matrix_t m, vector_t *f, number_t omega, number_t eps)
{
        /* Готовим память для векторов итераций и для подсчёта невязки */
        vector_t x1, x2, diff, tmp_vect;
        x1 = vector_create(m.size);
        x2 = vector_create(m.size);
        diff = vector_create(m.size);
        
        number_t norm = 0;
        int iters = 0;
        do {
                /* Проводим итерацию метода */
                relax_iteration(m, *f, x1, x2, omega);

                /* Вычисление невязки: считаем A*x2 и помещаем в x1 */ 
                matrix_mulMatVector(m, x2, x1);
                /* Вычитаем f */
                vector_sub(x1, *f, diff);
                /* Считаем норму невязки */
                norm = vector_norm(diff);

                /* Меняем местами старый и новый вектора для следующей итерации */
                tmp_vect = x1;
                x1 = x2;
                x2 = tmp_vect;

                /* Выводим сообщение о текущей итерации (анализ сходимости) */
                fprintf(stderr, "[RELAX] Iteration %d, residual " NUMBER_WRITE_FORMAT "\n", iters + 1, norm);

        } while (iters++ < 50 && norm >= eps);

        /* Помещаем результат в f, откуда его заберут снаружи */
        for (int i=0; i<f->size; i++) {
                f->vector[i] = x1.vector[i];
        }

        /* Освобождаем память */
        vector_free(x1);
        vector_free(x2);
        vector_free(diff);
}
