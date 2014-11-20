#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#define _GNU_SOURCE
#include <getopt.h>

#include "matrix.h"

#include "gauss.h"
#include "gauss_mod.h"
#include "relax.h"

static struct option longopts[] = {
        { .name = "help", .has_arg = 0, .flag = NULL, .val = 'h' },
        { .name = "method", .has_arg = 1, .flag = NULL, .val = 'm' },
        { .name = "operation", .has_arg = 1, .flag = NULL, .val = 'o' },
        { .name = "format", .has_arg = 1, .flag = NULL, .val = 'f' }
};

static enum m_methods {
        METHOD_GAUSS = 0,
        METHOD_GAUSS_MOD = 1,
        METHOD_RELAXATION = 2,
        METHOD_END
} method;

static char *method_names[] = {
        "gauss", "gauss_mod", "relax"
};

static char *ops_names[] = {
        "det", "solve", "invert"
};

static enum m_ops {
        OPERATION_DET = 0,
        OPERATION_SOLVE = 1,
        OPERATION_INVERT = 2,
        OPERATION_END
} operation;

static char *format_names[] = {
        "text", "latex"
};

static format_t format;

void op_solve(enum m_methods met, format_t format)
{
        if (met < 0 || met >= METHOD_END) {
                fprintf(stderr, "[ERROR] Unknown method: %d\n", met);
                return;
        }

        fprintf(stderr, "[INPUT] Type N, than matrix (NxN)\n");
        matrix_t m = matrix_read(stdin);
        
        fprintf(stderr, "[INPUT] Type vector f (length %d)\n", m.size);
        vector_t f = vector_readN(stdin, m.size);

        number_t omega = 1;
        number_t eps = 0.0001;
        switch (met) {
                case METHOD_GAUSS:
                        gauss_solve(m, &f);
                        break;
                case METHOD_GAUSS_MOD:
                        gauss_mod_solve(m, &f);
                        break;
                case METHOD_RELAXATION:
                        fprintf(stderr, "[INPUT] Type 'omega' for over relaxation method\n");
                        scanf(NUMBER_READ_FORMAT, &omega);
                        fprintf(stderr, "[INPUT] Type precision coefficient (eps)\n");
                        scanf(NUMBER_READ_FORMAT, &eps);
                        relax_solve(m, &f, omega, eps);
                        break;
        }

        fprintf(stderr, "[OUTPUT] Result\n");
        vector_print(stdout, f, format);

        matrix_free(m);
        vector_free(f);
}

void op_det(enum m_methods met, format_t format)
{
        if (met != METHOD_GAUSS && met != METHOD_GAUSS_MOD) {
                fprintf(stderr, "[INPUT] Determinant calculation methods: gauss and gauss_mod\n");
                return;
        }

        fprintf(stderr, "[INPUT] Type N, than matrix (NxN)\n");
        matrix_t m = matrix_read(stdin);
        vector_t f = vector_create(m.size);

        number_t det = 0;

        switch (met) {
                case METHOD_GAUSS:
                        det = gauss_solve(m, &f);
                        break;
                case METHOD_GAUSS_MOD:
                        det = gauss_mod_solve(m, &f);
                        break;
        }

        fprintf(stderr, "[OUTPUT] Result\n");
        printf(NUMBER_WRITE_FORMAT "\n", det);

        matrix_free(m);
        vector_free(f);
}

void op_invert(enum m_methods met, format_t format)
{
        if (met != METHOD_GAUSS && met != METHOD_GAUSS_MOD) {
                fprintf(stderr, "[INPUT] Invertion methods: gauss and gauss_mod\n");
                return;
        }

        fprintf(stderr, "[INPUT] Type N, than matrix (NxN)\n");
        matrix_t m = matrix_read(stdin);
        matrix_t inv;

        switch (met) {
                case METHOD_GAUSS:
                        inv = gauss_invert(m);
                        break;
                case METHOD_GAUSS_MOD:
                        inv = gauss_mod_invert(m);
                        break;
        }

        fprintf(stderr, "[OUTPUT] Result\n");
        matrix_print(stdout, inv, format);

        matrix_free(m);
        matrix_free(inv);
}

char *argv0;

void print_help()
{
        fprintf(stderr, "Usage: %s -o <operation> -m <method> [-f <format>]\n\n", argv0);
        fprintf(stderr, " -o, --operation=<operation>\tOperation name: det, solve, invert\n");
        fprintf(stderr, " -m, --method=<method>\t\tMethod: gauss, gauss_mod, relax (only for solve operation)\n");
        fprintf(stderr, " -f, --format=<format>\t\tOutput format: text, latex\n");
        fprintf(stderr, " -h, --help\t\t\tPrint this message\n\n");
}


int main(int argc, char *argv[])
{
        int c;
        int flag_gotmet = 0, flag_gotop = 0, flag_gotformat = 0;
        argv0 = argv[0];
        
        while ((c = getopt_long(argc, argv, "hm:o:f:", longopts, NULL)) > 0) {
                switch (c) {
                        case 'h':
                                print_help();
                                exit(0);
                        case 'm':
                                flag_gotmet = 1;
                                method = 0;

                                while (method != METHOD_END && strcmp(optarg, method_names[method]))
                                        method++;

                                if (method == METHOD_END) {
                                        fprintf(stderr, "ERROR: Unknown method: %s\n\n", optarg);
                                        print_help();
                                        exit(1);
                                }
                                break;
                        case 'o':
                                flag_gotop = 1;
                                operation = 0;

                                while (operation != OPERATION_END && strcmp(optarg, ops_names[operation]))
                                        operation++;

                                if (operation == OPERATION_END) {
                                        fprintf(stderr, "ERROR: Unknown operation: %s\n\n", optarg);
                                        print_help();
                                        exit(1);
                                }
                                break;
                        case 'f':
                                flag_gotformat = 1;
                                format = 0;

                                while (format != FORMAT_END && strcmp(optarg, format_names[format]))
                                        format++;

                                if (format == FORMAT_END) {
                                        fprintf(stderr, "ERROR: Unknown format: %s\n\n", optarg);
                                        print_help();
                                        exit(1);
                                }
                                break;
                }
        }

        if (!flag_gotmet || !flag_gotop) {
                print_help();
                exit(0);
        }

        if (!flag_gotformat) {
                format = FORMAT_TEXT;
        }

        /* select operation and method */
        switch (operation) {
                case OPERATION_SOLVE:
                        op_solve(method, format);
                        break;
                case OPERATION_DET:
                        op_det(method, format);
                        break;
                case OPERATION_INVERT:
                        op_invert(method, format);
                        break;
        }

        return 0;
}
