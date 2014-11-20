#ifndef INCLUDE_MATRIX_H
#define INCLUDE_MATRIX_H

#include <stdio.h>
#include <stdlib.h>

typedef double number_t;
#define NUMBER_READ_FORMAT "%lf"
#define NUMBER_WRITE_FORMAT "%.5g"
#define ZERO_EPS (0.000001)
#define NOT_ZERO(n) (((n) >= ZERO_EPS) || ((n) <= -ZERO_EPS))

#define matrix_elem(m, i, j) ((m).matrix[(m).rows[(i)]][(m).cols[(j)]])

typedef struct {
        unsigned int size;
        unsigned int *rows;
        unsigned int *cols;
        number_t **matrix;
} matrix_t;

typedef struct {
        unsigned int size;
        number_t *vector;
} vector_t;

typedef enum m_formats {
        FORMAT_TEXT = 0,
        FORMAT_LATEX = 1,
        FORMAT_END = 2
} format_t;

matrix_t matrix_create(int n);
matrix_t matrix_read(FILE *stream);
matrix_t matrix_readN(FILE *stream, int n);
void matrix_print(FILE *stream, matrix_t m, format_t f);
matrix_t matrix_copy(matrix_t source);
void matrix_exchangeRows(matrix_t m, int a, int b);
void matrix_swapCols(matrix_t m, int a, int b);
void matrix_divRow(matrix_t m, int j, number_t div);
void matrix_hitRow(matrix_t m, int src, int dst, number_t mul);
void matrix_mulMatVector(matrix_t m, vector_t f, vector_t result);

vector_t vector_create(int n);
vector_t vector_read(FILE *stream);
vector_t vector_readN(FILE *stream, int n);
vector_t vector_copy(vector_t source);
void vector_print(FILE *stream, vector_t v, format_t f);
void vector_exchangeElems(vector_t v, int a, int b);
void vector_sub(vector_t a, vector_t b, vector_t result);
number_t vector_norm(vector_t v);

void matrix_free(matrix_t matrix);
void vector_free(vector_t vector);

#endif
