#ifndef INCLUDE_RELAX_H
#define INCLUDE_RELAX_H

#include "matrix.h"

void relax_solve(matrix_t m, vector_t *f, number_t omega, number_t eps);

#endif
