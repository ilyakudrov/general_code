#ifndef __EIGEN_H__
#define __EIGEN_H__

#ifdef DOUBLE
#define FLOAT double
#else
#define FLOAT float
#endif

#include "link.h"
#include <cmath>

typedef struct scomplex_number {
  FLOAT re;
  FLOAT im;
} complex;

std::vector<su2> make_matrix_staggered(const std::vector<su2> &conf,
                                       FLOAT mu_q);
std::vector<complex>
matrix_multiplication_staggered(std::vector<su2> &matrix,
                                const std::vector<complex> &vec_input);
void matrix_multiplication_su2(const complex *vec_input, complex *vec_output,
                               su2 &A);
FLOAT eta_sign(int mu, link1 &link);

#endif