#ifndef INCLUDE_GAUSS_H
#define INCLUDE_GAUSS_H

#include "matrix.h"

typedef struct {
        void (*swap)(int, int, void*);
        void (*div)(int, number_t, void*);
        void (*hit)(int, int, number_t, void*);
} gauss_handlers_t;

number_t gauss_full(matrix_t m, gauss_handlers_t *h, void *arg);

number_t gauss_solve(matrix_t m, vector_t *f);

matrix_t gauss_invert(matrix_t m);

number_t gauss_direct(matrix_t m, gauss_handlers_t *h, void *arg);
void gauss_reverse(matrix_t m, gauss_handlers_t *h, void *arg);

#endif
