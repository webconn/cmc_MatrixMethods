#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include "gauss.h"

int main(int argc, char *argv[])
{
        /* 1. read matrix */
        printf("Enter N, then matrix NxN:\n");
        matrix_t m = matrix_read(stdin);

        matrix_t inv = gauss_invert(m);

        printf("Your matrix is:\n");
        matrix_print(stdout, m);

        printf("Result:\n");
        matrix_print(stdout, inv);

        matrix_free(m);
        matrix_free(inv);

        return 0;
}
