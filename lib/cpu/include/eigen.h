#ifndef __EIGEN_H__
#define __EIGEN_H__

// #include <iostream>
#include "link.h"
#include <cmath>
// #include "matrix.h"
// #include <complex.h>
// #include <sys/time.h>

typedef struct scomplex_s {
    FLOAT re;
    FLOAT im;
} scomplex_t;

void dirac_mult(scomplex_t* res, const scomplex_t* src, int place, matrix* data_conf, int x_size, int y_size, int z_size, int t_size, FLOAT mass, FLOAT mu_q);
void dirac_mult_conj(scomplex_t* res, const scomplex_t* src, int place, matrix* data_conf, int x_size, int y_size, int z_size, int t_size);
int complex_place(link1<matrix>& link);
// void matrix_mult_complex(matrix A, const scomplex_t* a, scomplex_t* a1);
void matrix_mult_complex1(matrix A, const scomplex_t* a, scomplex_t* a1, int i, FLOAT border_sign);
FLOAT test_module(const scomplex_t* vec, int size);// for testing eigenvectors
void test_eigenvector(const scomplex_t* eigenvector, scomplex_t eigenvalue, int size, matrix* data_conf, int x_size, int y_size, int z_size, int t_size, FLOAT mass, FLOAT mu_q, FLOAT tolerance);
FLOAT eta_sign(int mu, link1<matrix>& link);
FLOAT eta_sign_5(link1<matrix>& link);
// void scomplex_mult_add(scomplex_t* a1, const scomplex_t* b, const scomplex_t* c, scomplex_t* z1, scomplex_t* z2);

#endif