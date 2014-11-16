#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

matrix_t matrix_read(FILE *stream)
{
        /* 1. read matrix size */
        int N;
        fscanf(stream, "%d", &N);
        return matrix_readN(stream, N);
}

matrix_t matrix_create(int n)
{
        return matrix_readN(NULL, n);
}

matrix_t matrix_readN(FILE *stream, int N)
{
        /* 2. allocate memory for sequences */
        unsigned int *rows, *cols;
        rows = (unsigned int *) malloc(N * sizeof (unsigned int));
        cols = (unsigned int *) malloc(N * sizeof (unsigned int));

        /* 2. allocate memory for matrix */
        number_t **matrix = (number_t **) malloc(N * sizeof (number_t *));

        for (unsigned int i=0; i<N; i++) {
                matrix[i] = (number_t *) malloc(N * sizeof (number_t));
        }

        /* 3. read matrix data */
        for (unsigned int i=0; i<N; i++) {
                for (unsigned int j=0; j<N; j++) {
                        if (stream != NULL && fscanf(stream, NUMBER_READ_FORMAT, &matrix[i][j]) == 0) {
                                fprintf(stderr, "ERROR: Wrong input stream (unexpected EOF)\n");
                                exit(1);
                        } else if (stream == NULL) {
                                matrix[i][j] = 0;
                        }
                }
                rows[i] = i;
                cols[i] = i;
        }

        matrix_t ret = {
                .matrix = matrix,
                .size = N,
                .cols = cols,
                .rows = rows
        };

        return ret;
}

void matrix_print(FILE *stream, matrix_t m, format_t f)
{
        for (int i=0; i<m.size; i++) {
                for (int j=0; j<m.size; j++) {
                        fprintf(stream, NUMBER_WRITE_FORMAT " ", m.matrix[m.rows[i]][m.cols[j]]);
                }
                fputc('\n', stream);
        }
}

matrix_t matrix_copy(matrix_t source)
{
        number_t **m = (number_t **) malloc(source.size * sizeof (number_t *));
        unsigned int *rows = (unsigned int *) malloc(source.size * sizeof (unsigned int));
        unsigned int *cols = (unsigned int *) malloc(source.size * sizeof (unsigned int));

        for (int i=0; i<source.size; i++) {
                m[i] = (number_t *) malloc(source.size * sizeof (number_t));

                for (int j=0; j<source.size; j++) {
                        m[i][j] = source.matrix[i][j];
                }

                rows[i] = source.rows[i];
                cols[i] = source.cols[i];
        }

        matrix_t ret = {
                .matrix = m,
                .size = source.size,
                .rows = rows,
                .cols = cols
        };

        return ret;
}

void matrix_divRow(matrix_t m, int j, number_t div)
{
        if (!NOT_ZERO(div))
                return;

        for (int i = 0; i < m.size; i++) {
                m.matrix[m.rows[j]][i] /= div;
        }
}

void matrix_hitRow(matrix_t m, int src, int dst, number_t mul)
{
        for (int i = 0; i < m.size; i++) {
                m.matrix[m.rows[dst]][i] += mul * m.matrix[m.rows[src]][i];
        }
}

vector_t vector_create(int N)
{
        return vector_readN(NULL, N);
}

vector_t vector_read(FILE *stream)
{
        int n;
        fscanf(stream, "%d", &n);
        return vector_readN(stream, n);
}

vector_t vector_readN(FILE *stream, int n)
{
        number_t *vect = (number_t *) malloc(n * sizeof (number_t));

        for (int i = 0; i < n; i++) {
                if (stream) {
                        fscanf(stream, NUMBER_READ_FORMAT, &vect[i]);
                } else {
                        vect[i] = 0;
                }
        }

        vector_t v = {
                .size = n,
                .vector = vect
        };

        return v;
}

void vector_print(FILE *stream, vector_t v, format_t f)
{
        for (int i = 0; i < v.size; i++) {
                fprintf(stream, NUMBER_WRITE_FORMAT " ", v.vector[i]);
        }
        fputc('\n', stream);
}

vector_t vector_copy(vector_t source)
{
        vector_t ret = vector_create(source.size);

        for (int i=0; i<ret.size; i++) {
                ret.vector[i] = source.vector[i];
        }

        return ret;
}

void vector_exchangeElems(vector_t v, int a, int b)
{
        number_t tmp = v.vector[a];
        v.vector[a] = v.vector[b];
        v.vector[b] = tmp;
}

void matrix_exchangeRows(matrix_t m, int a, int b)
{
        unsigned int tmp = m.rows[a];
        m.rows[a] = m.rows[b];
        m.rows[b] = tmp;
}

void matrix_swapCols(matrix_t m, int a, int b)
{
        unsigned int tmp = m.cols[a];
        m.cols[a] = m.cols[b];
        m.cols[b] = tmp;
}

void vector_free(vector_t v)
{
        free(v.vector);
}

void matrix_free(matrix_t m)
{
        for (unsigned int i=0; i<m.size; i++)
                free(m.matrix[i]);
        free(m.matrix);
        free(m.rows);
        free(m.cols);
}
