#include "gauss_mod.h"

number_t gauss_mod_full(matrix_t m, gauss_mod_handlers_t *h, void *arg)
{
        return 0;
}

number_t gauss_mod_solve(matrix_t m, vector_t *f)
{
        return 0;
}

matrix_t gauss_mod_invert(matrix_t m)
{
        return matrix_create(m.size);
}

number_t gauss_mod_direct(matrix_t m, gauss_mod_handlers_t *h, void *arg)
{
        return 0;
}

void gauss_mod_reverse(matrix_t m, gauss_mod_handlers_t *h, void *arg)
{

}
