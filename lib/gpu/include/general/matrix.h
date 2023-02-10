#pragma once

#include "cuda_host_device.h"

namespace matrix_gpu {

struct complex_t {
  double real;
  double imag;

  CUDA_DEVICE complex_t(const double real1, const double imag1);
  CUDA_DEVICE complex_t();

  CUDA_DEVICE double module();

  CUDA_DEVICE double norm2();

  CUDA_DEVICE double angle();

  CUDA_DEVICE complex_t conj();

  CUDA_DEVICE complex_t negative();

  CUDA_DEVICE complex_t mult_by_imag(double x);

  CUDA_DEVICE complex_t sqrt_complex();

  CUDA_DEVICE complex_t &operator+=(const complex_t &a);

  CUDA_DEVICE complex_t &operator-=(const complex_t &a);

  CUDA_DEVICE complex_t &operator*=(const complex_t &a);

  CUDA_DEVICE complex_t &operator*=(const double &a);

  CUDA_DEVICE complex_t &operator/=(const complex_t &a);

  CUDA_DEVICE complex_t &operator/=(const double &a);
};

CUDA_DEVICE complex_t operator+(const complex_t &a, const complex_t &b);

CUDA_DEVICE complex_t operator-(const complex_t &a, const complex_t &b);

CUDA_DEVICE complex_t operator*(const complex_t &a, const complex_t &b);

CUDA_DEVICE complex_t operator*(const double &a, const complex_t &b);

CUDA_DEVICE complex_t operator*(const complex_t &a, const double &b);

CUDA_DEVICE complex_t operator^(const complex_t &a, const complex_t &b);

CUDA_DEVICE complex_t operator%(const complex_t &a, const complex_t &b);

CUDA_DEVICE complex_t operator&(const complex_t &a, const complex_t &b);

CUDA_DEVICE complex_t operator/(const complex_t &a, const double &b);

CUDA_DEVICE complex_t operator/(const complex_t &a, const complex_t &b);

// su2 matrix in sigma matrices representation (a0 + ai * sigma[i])
class su2 {
public:
  double a0, a1, a2, a3;
  CUDA_DEVICE su2(double b0, double b1, double b2, double b3);
  CUDA_DEVICE su2();

  // calculate trace of the matrix
  CUDA_DEVICE double tr();

  CUDA_DEVICE double multiply_tr(const su2 &B);

  CUDA_DEVICE double multiply_tr_adjoint(const su2 &B);

  // calculate inverse of the matrix
  CUDA_DEVICE su2 inverse();

  // calculate conjugate of the matrix
  CUDA_DEVICE su2 conj() const;

  // gets projection onto su2 group
  CUDA_DEVICE su2 proj();

  // calculates module of vector in sigma matrices representation
  CUDA_DEVICE double module();

  // left and right multiplication of the matrix on sigma[3] matrix
  CUDA_DEVICE su2 sigma3_mult() const;

  CUDA_DEVICE su2 &operator+=(const su2 &A);
};

CUDA_DEVICE su2 operator+(const su2 &A, const su2 &B);
CUDA_DEVICE su2 operator-(const su2 &A, const su2 &B);
CUDA_DEVICE su2 operator*(const double &x, const su2 &A);
CUDA_DEVICE su2 operator*(const su2 &A, const double &x);

// matrix multiplication A * B
CUDA_DEVICE su2 operator*(const su2 &A, const su2 &B);

// matrix multiplication A * B.conj()
CUDA_DEVICE su2 operator^(const su2 &A, const su2 &B);

// matrix multiplication A.conj() * B
CUDA_DEVICE su2 operator%(const su2 &A, const su2 &B);

// abelian variable (module * exp(i * phi))
class abelian {
public:
  double r, phi;
  CUDA_DEVICE abelian();
  CUDA_DEVICE abelian(double r1, double phi1);

  // trace
  CUDA_DEVICE double tr();

  CUDA_DEVICE double multiply_tr(const abelian &B);

  CUDA_DEVICE double multiply_tr_adjoint(const abelian &B);

  // inverse
  CUDA_DEVICE abelian inverse();
  // conjugated
  CUDA_DEVICE abelian conj() const;
  CUDA_DEVICE abelian proj();
  CUDA_DEVICE double module();
};

CUDA_DEVICE abelian operator+(const abelian &A, const abelian &B);
CUDA_DEVICE abelian operator-(const abelian &A, const abelian &B);
CUDA_DEVICE abelian operator*(const double &x, const abelian &A);
CUDA_DEVICE abelian operator*(const abelian &A, const double &x);

CUDA_DEVICE abelian operator*(const abelian &A, const abelian &B);
CUDA_DEVICE abelian operator^(const abelian &A, const abelian &B);
CUDA_DEVICE abelian operator%(const abelian &A, const abelian &B);

// su3 matrix in 3x3 complex matrix representation
class su3 {
public:
  complex_t matrix[3][3];
  CUDA_DEVICE su3(complex_t B[3][3]);
  CUDA_DEVICE su3();

  // calculate trace of the matrix
  CUDA_DEVICE double tr();

  CUDA_DEVICE complex_t tr_complex();

  CUDA_DEVICE double multiply_tr(const su3 &B);

  CUDA_DEVICE double multiply_tr_adjoint(const su3 &B);

  // calculate inverse of the matrix
  CUDA_DEVICE su3 inverse();

  // calculate conjugate of the matrix
  CUDA_DEVICE su3 conj() const;

  CUDA_DEVICE su3 mult_by_imag(double x);

  // gets projection onto su3 group
  CUDA_DEVICE su3 proj();

  // muliply from left and right by lambda3
  CUDA_DEVICE su3 lambda3_mult();

  // muliply from left and right by lambda8
  CUDA_DEVICE su3 lambda8_mult();

  // calculates module of vector in sigma matrices representation
  CUDA_DEVICE double module();

  CUDA_DEVICE complex_t determinant();

  CUDA_DEVICE complex_t unitarity_check();
};

CUDA_DEVICE su3 operator+(const su3 &A, const su3 &B);
CUDA_DEVICE su3 operator-(const su3 &A, const su3 &B);
CUDA_DEVICE su3 operator*(const double &x, const su3 &A);
CUDA_DEVICE su3 operator*(const su3 &A, const double &x);

CUDA_DEVICE su3 operator*(const complex_t &x, const su3 &A);
CUDA_DEVICE su3 operator*(const su3 &A, const complex_t &x);

// matrix multiplication A * B
CUDA_DEVICE su3 operator*(const su3 &A, const su3 &B);

// matrix multiplication A * B.conj()
CUDA_DEVICE su3 operator^(const su3 &A, const su3 &B);

// matrix multiplication A.conj() * B
CUDA_DEVICE su3 operator%(const su3 &A, const su3 &B);

// diagonal su3 matrix
class su3_abelian {
public:
  complex_t matrix[3];
  CUDA_DEVICE su3_abelian(complex_t B[3]);
  CUDA_DEVICE su3_abelian();

  // calculate trace of the matrix
  CUDA_DEVICE double tr();

  CUDA_DEVICE complex_t tr_complex();

  CUDA_DEVICE double multiply_tr(const su3_abelian &B);

  CUDA_DEVICE double multiply_tr_adjoint(const su3_abelian &B);

  // calculate inverse of the matrix
  CUDA_DEVICE su3_abelian inverse();

  // calculate conjugate of the matrix
  CUDA_DEVICE su3_abelian conj() const;

  CUDA_DEVICE su3_abelian mult_by_imag(double x);

  // gets projection onto su3 group
  CUDA_DEVICE su3_abelian proj();

  // muliply from left and right by lambda3
  CUDA_DEVICE su3_abelian lambda3_mult();

  // muliply from left and right by lambda8
  CUDA_DEVICE su3_abelian lambda8_mult();

  // calculates module of vector in sigma matrices representation
  CUDA_DEVICE double module();

  CUDA_DEVICE complex_t determinant();

  CUDA_DEVICE complex_t unitarity_check();
};

CUDA_DEVICE su3_abelian operator+(const su3_abelian &A, const su3_abelian &B);
CUDA_DEVICE su3_abelian operator-(const su3_abelian &A, const su3_abelian &B);
CUDA_DEVICE su3_abelian operator*(const double &x, const su3_abelian &A);
CUDA_DEVICE su3_abelian operator*(const su3_abelian &A, const double &x);

CUDA_DEVICE su3_abelian operator*(const complex_t &x, const su3_abelian &A);
CUDA_DEVICE su3_abelian operator*(const su3_abelian &A, const complex_t &x);

// matrix multiplication A * B
CUDA_DEVICE su3_abelian operator*(const su3_abelian &A, const su3_abelian &B);

// matrix multiplication A * B.conj()
CUDA_DEVICE su3_abelian operator^(const su3_abelian &A, const su3_abelian &B);

// matrix multiplication A.conj() * B
CUDA_DEVICE su3_abelian operator%(const su3_abelian &A, const su3_abelian &B);

// 3D vector realisation for spin model.
class spin {
public:
  double a1, a2, a3;

  CUDA_DEVICE spin(double a1, double a2, double a3);
  CUDA_DEVICE spin(su2 U);
  CUDA_DEVICE spin();

  CUDA_DEVICE double norm();
  // set norm equal to unity
  CUDA_DEVICE void normalize();
  // reflect spin-vector through the V-vector axis
  CUDA_DEVICE double reflect(spin &V);
  // reflect without difference calculation
  CUDA_DEVICE void reflect_fast(spin &V);
  // rorate spin-vector to the same direction as V-vector
  CUDA_DEVICE double parallel(spin &V);
  // parallel without difference calculation
  CUDA_DEVICE void parallel_fast(spin &V);
  // check identity matrix
  CUDA_DEVICE bool IsUnit();
  // contribution from neighbour site
  CUDA_DEVICE spin contribution(const su2 &A) const;
  // contribution from neighbour site in negative direction (conjugated matrix)
  CUDA_DEVICE spin contribution_conj(const su2 &A) const;
  CUDA_DEVICE void contribution1(const su2 &A, const spin &b);
  CUDA_DEVICE void contribution1_conj(const su2 &A, const spin &b);
  // calculation of gauge matrix G from spin vector.
  CUDA_DEVICE su2 GetGaugeMatrix();
};

CUDA_DEVICE spin operator+(const spin &A, const spin &B);
CUDA_DEVICE spin operator-(const spin &A, const spin &B);
CUDA_DEVICE spin operator*(const double &x, const spin &A);
CUDA_DEVICE spin operator*(const spin &A, const double &x);
CUDA_DEVICE double operator*(const spin &A, const spin &B);

CUDA_DEVICE complex_t::complex_t(const double real1, const double imag1) {
  real = real1;
  imag = imag1;
}

CUDA_DEVICE complex_t::complex_t() {
  real = 0;
  imag = 0;
}

CUDA_DEVICE double complex_t::module() {
  return sqrt(real * real + imag * imag);
}

CUDA_DEVICE double complex_t::norm2() { return real * real + imag * imag; }

CUDA_DEVICE double complex_t::angle() { return atan2(imag, real); }

CUDA_DEVICE complex_t complex_t::conj() { return complex_t(real, -imag); }

CUDA_DEVICE complex_t complex_t::negative() { return complex_t(-real, -imag); }

CUDA_DEVICE complex_t complex_t::mult_by_imag(double x) {
  return complex_t(-imag * x, real * x);
}

CUDA_DEVICE complex_t complex_t::sqrt_complex() {
  double module = sqrt(sqrt(real * real + imag * imag));

  double phi = atan2(imag, real) / 2;

  return complex_t(module * cos(phi), module * sin(phi));
}

CUDA_DEVICE complex_t operator+(const complex_t &a, const complex_t &b) {
  return complex_t(a.real + b.real, a.imag + b.imag);
}

CUDA_DEVICE complex_t operator-(const complex_t &a, const complex_t &b) {
  return complex_t(a.real - b.real, a.imag - b.imag);
}

CUDA_DEVICE complex_t operator*(const complex_t &a, const complex_t &b) {
  return complex_t(a.real * b.real - a.imag * b.imag,
                   a.real * b.imag + a.imag * b.real);
}

CUDA_DEVICE complex_t operator*(const double &a, const complex_t &b) {
  return complex_t(a * b.real, a * b.imag);
}

CUDA_DEVICE complex_t operator*(const complex_t &a, const double &b) {
  return complex_t(a.real * b, a.imag * b);
}

CUDA_DEVICE complex_t operator^(const complex_t &a, const complex_t &b) {
  return complex_t(a.real * b.real + a.imag * b.imag,
                   a.imag * b.real - a.real * b.imag);
}

CUDA_DEVICE complex_t operator%(const complex_t &a, const complex_t &b) {
  return complex_t(a.real * b.real + a.imag * b.imag,
                   a.real * b.imag - a.imag * b.real);
}

CUDA_DEVICE complex_t operator&(const complex_t &a, const complex_t &b) {
  return complex_t(a.real * b.real - a.imag * b.imag,
                   -a.imag * b.real - a.real * b.imag);
}

CUDA_DEVICE complex_t operator/(const complex_t &a, const double &b) {
  return complex_t(a.real / b, a.imag / b);
}

CUDA_DEVICE complex_t operator/(const complex_t &a, const complex_t &b) {
  double norm_factor = 1. / (b.real * b.real + b.imag * b.imag);
  return complex_t((a.real * b.real + a.imag * b.imag) * norm_factor,
                   (a.imag * b.real - a.real * b.imag) * norm_factor);
}

CUDA_DEVICE complex_t &complex_t::operator+=(const complex_t &a) {
  real += a.real;
  imag += a.imag;

  return *this;
}

CUDA_DEVICE complex_t &complex_t::operator-=(const complex_t &a) {
  real -= a.real;
  imag -= a.imag;

  return *this;
}

CUDA_DEVICE complex_t &complex_t::operator*=(const complex_t &a) {
  double tmp1 = real * a.real - imag * a.imag;
  imag = real * a.imag + imag * a.real;
  real = tmp1;

  return *this;
}

CUDA_DEVICE complex_t &complex_t::operator*=(const double &a) {
  real *= a;
  imag *= a;

  return *this;
}

CUDA_DEVICE complex_t &complex_t::operator/=(const complex_t &a) {
  double norm_factor = 1. / (a.real * a.real + a.imag * a.imag);
  double tmp = (real * a.real + imag * a.imag) * norm_factor;

  imag = (imag * a.real - real * a.imag) * norm_factor;
  real = tmp;

  return *this;
}

CUDA_DEVICE complex_t &complex_t::operator/=(const double &a) {
  real /= a;
  imag /= a;

  return *this;
}

// su2 methods
CUDA_DEVICE su2::su2() {
  a0 = 1.;
  a1 = 0;
  a2 = 0;
  a3 = 0;
}

CUDA_DEVICE su2::su2(double b0, double b1, double b2, double b3) {
  a0 = b0;
  a1 = b1;
  a2 = b2;
  a3 = b3;
}

CUDA_DEVICE double su2::tr() { return a0; }

CUDA_DEVICE double su2::multiply_tr(const su2 &B) {
  return a0 * B.a0 + a1 * B.a1 + a2 * B.a2 + a3 * B.a3;
}

CUDA_DEVICE double su2::multiply_tr_adjoint(const su2 &B) {
  double trace = a0 * B.a0 + a1 * B.a1 + a2 * B.a2 + a3 * B.a3;
  return trace * trace - 1;
}

CUDA_DEVICE su2 su2::inverse() {
  double rho = a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3;
  return su2(a0 / rho, -a1 / rho, -a2 / rho, -a3 / rho);
}
CUDA_DEVICE double su2::module() {
  return a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3;
}
CUDA_DEVICE su2 su2::conj() const { return su2(a0, -a1, -a2, -a3); }
CUDA_DEVICE su2 su2::proj() {
  double rho = sqrt(a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3);
  return su2(a0 / rho, a1 / rho, a2 / rho, a3 / rho);
}
CUDA_DEVICE su2 su2::sigma3_mult() const { return su2(a0, -a1, -a2, a3); }

CUDA_DEVICE su2 operator+(const su2 &A, const su2 &B) {
  return su2(A.a0 + B.a0, A.a1 + B.a1, A.a2 + B.a2, A.a3 + B.a3);
}
CUDA_DEVICE su2 operator-(const su2 &A, const su2 &B) {
  return su2(A.a0 - B.a0, A.a1 - B.a1, A.a2 - B.a2, A.a3 - B.a3);
}
CUDA_DEVICE su2 operator*(const double &x, const su2 &A) {
  return su2(A.a0 * x, A.a1 * x, A.a2 * x, A.a3 * x);
}
CUDA_DEVICE su2 operator*(const su2 &A, const double &x) {
  return su2(A.a0 * x, A.a1 * x, A.a2 * x, A.a3 * x);
}

CUDA_DEVICE su2 operator*(const su2 &A, const su2 &B) {
  return su2(A.a0 * B.a0 - A.a1 * B.a1 - A.a2 * B.a2 - A.a3 * B.a3,
             A.a0 * B.a1 + B.a0 * A.a1 + A.a3 * B.a2 - A.a2 * B.a3,
             A.a0 * B.a2 + B.a0 * A.a2 + A.a1 * B.a3 - A.a3 * B.a1,
             A.a0 * B.a3 + B.a0 * A.a3 + A.a2 * B.a1 - A.a1 * B.a2);
}
CUDA_DEVICE su2 operator^(const su2 &A, const su2 &B) {
  return su2(A.a0 * B.a0 + A.a1 * B.a1 + A.a2 * B.a2 + A.a3 * B.a3,
             -A.a0 * B.a1 + B.a0 * A.a1 - A.a3 * B.a2 + A.a2 * B.a3,
             -A.a0 * B.a2 + B.a0 * A.a2 - A.a1 * B.a3 + A.a3 * B.a1,
             -A.a0 * B.a3 + B.a0 * A.a3 - A.a2 * B.a1 + A.a1 * B.a2);
}
CUDA_DEVICE su2 operator%(const su2 &A, const su2 &B) {
  return su2(A.a0 * B.a0 + A.a1 * B.a1 + A.a2 * B.a2 + A.a3 * B.a3,
             A.a0 * B.a1 - B.a0 * A.a1 - A.a3 * B.a2 + A.a2 * B.a3,
             A.a0 * B.a2 - B.a0 * A.a2 - A.a1 * B.a3 + A.a3 * B.a1,
             A.a0 * B.a3 - B.a0 * A.a3 - A.a2 * B.a1 + A.a1 * B.a2);
}

CUDA_DEVICE su2 &su2::operator+=(const su2 &A) {
  a0 += A.a0;
  a1 += A.a1;
  a2 += A.a2;
  a3 += A.a3;

  return *this;
}

// abelian methods
CUDA_DEVICE abelian::abelian() {
  r = 1;
  phi = 0;
}

CUDA_DEVICE abelian::abelian(double r1, double phi1) {
  r = r1;
  phi = phi1;
}

CUDA_DEVICE double abelian::tr() { return r * cos(phi); }

CUDA_DEVICE double abelian::multiply_tr(const abelian &B) {
  return r * B.r * cos(phi - B.phi);
}

CUDA_DEVICE double abelian::multiply_tr_adjoint(const abelian &B) {
  double trace = r * B.r * cos(2 * (phi - B.phi));
  return trace * trace;
}

CUDA_DEVICE abelian abelian::inverse() { return abelian(1 / r, -phi); }
CUDA_DEVICE double abelian::module() { return r; }
CUDA_DEVICE abelian abelian::conj() const { return abelian(r, -phi); }
CUDA_DEVICE abelian abelian::proj() { return abelian(1, phi); }

CUDA_DEVICE abelian operator+(const abelian &A, const abelian &B) {
  return abelian(sqrt((A.r * sin(A.phi) + B.r * sin(B.phi)) *
                          (A.r * sin(A.phi) + B.r * sin(B.phi)) +
                      (A.r * cos(A.phi) + B.r * cos(B.phi)) *
                          (A.r * cos(A.phi) + B.r * cos(B.phi))),
                 atan2(A.r * sin(A.phi) + B.r * sin(B.phi),
                       A.r * cos(A.phi) + B.r * cos(B.phi)));
}
CUDA_DEVICE abelian operator-(const abelian &A, const abelian &B) {
  return abelian(sqrt((A.r * sin(A.phi) - B.r * sin(B.phi)) *
                          (A.r * sin(A.phi) - B.r * sin(B.phi)) +
                      (A.r * cos(A.phi) - B.r * cos(B.phi)) *
                          (A.r * cos(A.phi) - B.r * cos(B.phi))),
                 atan2(A.r * sin(A.phi) - B.r * sin(B.phi),
                       A.r * cos(A.phi) - B.r * cos(B.phi)));
}
CUDA_DEVICE abelian operator*(const double &x, const abelian &A) {
  return abelian(A.r * x, A.phi);
}
CUDA_DEVICE abelian operator*(const abelian &A, const double &x) {
  return abelian(A.r * x, A.phi);
}

CUDA_DEVICE abelian operator*(const abelian &A, const abelian &B) {
  return abelian(A.r * B.r, A.phi + B.phi);
}
CUDA_DEVICE abelian operator^(const abelian &A, const abelian &B) {
  return abelian(A.r * B.r, A.phi - B.phi);
}
CUDA_DEVICE abelian operator%(const abelian &A, const abelian &B) {
  return abelian(A.r * B.r, B.phi - A.phi);
}

// su3 methods
CUDA_DEVICE su3::su3() {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i != j)
        matrix[i][j] = complex_t(0, 0);
      else
        matrix[i][j] = complex_t(1., 0);
    }
  }
}

CUDA_DEVICE su3::su3(complex_t B[3][3]) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      matrix[i][j] = B[i][j];
    }
  }
}

CUDA_DEVICE double su3::tr() {
  return (matrix[0][0].real + matrix[1][1].real + matrix[2][2].real) / 3;
}

CUDA_DEVICE complex_t su3::tr_complex() {
  return (matrix[0][0] + matrix[1][1] + matrix[2][2]) / 3;
}

CUDA_DEVICE double su3::multiply_tr(const su3 &B) {
  double trace = 0;
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++) {
      trace += matrix[i][k].real * B.matrix[i][k].real +
               matrix[i][k].imag * B.matrix[i][k].imag;
    }
  }
  return trace / 3;
}

CUDA_DEVICE double multiply_tr(const su3 &A, const su3 &B) {
  double trace = 0;
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++) {
      trace += A.matrix[i][k].real * B.matrix[k][i].real -
               A.matrix[i][k].imag * B.matrix[k][i].imag;
    }
  }
  return trace / 3;
}

CUDA_DEVICE double su3::multiply_tr_adjoint(const su3 &B) {
  complex_t trace = complex_t(0, 0);
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++) {
      trace += matrix[i][k] ^ B.matrix[i][k];
    }
  }

  return trace.norm2() - 1;
}

CUDA_DEVICE su3 su3::inverse() {
  su3 B;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      B.matrix[i][j].real = matrix[j][i].real;
      B.matrix[i][j].imag = -matrix[j][i].imag;
    }
  }
  return B;
}

CUDA_DEVICE complex_t determinant_part(complex_t B[3][3], int i, int j, int k) {
  return B[0][i] * (B[1][j] * B[2][k] - B[2][j] * B[1][k]);
}

CUDA_DEVICE double su3::module() {
  complex_t determinant;
  determinant = determinant_part(matrix, 0, 1, 2);
  determinant = determinant + determinant_part(matrix, 1, 0, 2);
  determinant = determinant + determinant_part(matrix, 2, 0, 1);
  return determinant.real;
}

CUDA_DEVICE complex_t su3::determinant() {
  complex_t determinant;
  determinant = determinant_part(matrix, 0, 1, 2);
  determinant = determinant - determinant_part(matrix, 1, 0, 2);
  determinant = determinant + determinant_part(matrix, 2, 0, 1);
  return determinant;
}

CUDA_DEVICE complex_t su3::unitarity_check() {
  su3 A = this->conj();
  A = A * *this;

  complex_t diagonal = complex_t(0, 0), non_diagonal = complex_t(0, 0);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j)
        diagonal = diagonal + A.matrix[i][j];
      else
        non_diagonal = non_diagonal + A.matrix[i][j];
    }
  }
  return non_diagonal;
}

CUDA_DEVICE su3 su3::conj() const {
  su3 B;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      B.matrix[i][j].real = matrix[j][i].real;
      B.matrix[i][j].imag = -matrix[j][i].imag;
    }
  }
  return B;
}

CUDA_DEVICE su3 su3::mult_by_imag(double x) {
  su3 C;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C.matrix[i][j] = matrix[i][j].mult_by_imag(x);
    }
  }
  return C;
}

/*su3 su3::proj() {
  su3 A = *this;
  A = (1. / sqrt((A ^ A).tr() / 3)) * A;
  complex_t x;
  for (int i = 0; i < 4; i++) {
    A = 1.5 * A - 0.5 * ((A ^ A) * A);
    // A = A * (complex_t(1, 0) - (1. / 3) * (A.determinant() - complex_t(1,
    // 0)));
    x = complex_t(1, -(A).determinant().imag / 3);
    // A = A * x;
  }

  return A;
}*/

CUDA_DEVICE void project_su2(complex_t v[][2]) {
  complex_t u[2][2];

  // v*v^dagger
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      u[i][j] = complex_t(0, 0);
      for (int k = 0; k < 2; k++) {
        u[i][j] += (v[i][k] ^ v[j][k]);
      }
    }
  }

  double a = u[0][0].real + u[1][1].real;
  double b =
      sqrt((u[0][0].real - u[1][1].real) * (u[0][0].real - u[1][1].real) +
           4 * (u[0][1].real * u[1][0].real - u[0][1].imag * u[1][0].imag));

  double y1 = (a + b) / 2;
  double y2 = (a - b) / 2;

  complex_t w[2][2];

  w[0][0] = complex_t((y1 - u[1][1].real), 0) / u[1][0];
  w[0][1] = complex_t((y2 - u[1][1].real), 0) / u[1][0];

  double w_norm = 1 / sqrt(w[0][0].norm2() + 1);
  w[0][0] *= w_norm;
  w[1][0] = complex_t(w_norm, 0);
  w_norm = 1 / sqrt(w[0][1].norm2() + 1);
  w[0][1] *= w_norm;
  w[1][1] = complex_t(w_norm, 0);

  // v^dagger * v
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      u[i][j] = complex_t(0, 0);
      for (int k = 0; k < 2; k++) {
        u[i][j] += (v[k][i] % v[k][j]);
      }
    }
  }

  a = u[0][0].real + u[1][1].real;
  b = sqrt((u[0][0].real - u[1][1].real) * (u[0][0].real - u[1][1].real) +
           4 * (u[0][1].real * u[1][0].real - u[0][1].imag * u[1][0].imag));

  y1 = (a + b) / 2;
  y2 = (a - b) / 2;

  complex_t r[2][2];

  r[0][0] = complex_t((y1 - u[1][1].real), 0) / u[1][0];
  r[0][1] = complex_t((y2 - u[1][1].real), 0) / u[1][0];

  double r_norm = 1 / sqrt(r[0][0].norm2() + 1);
  r[0][0] *= r_norm;
  r[1][0] = complex_t(r_norm, 0);
  r_norm = 1 / sqrt(r[0][1].norm2() + 1);
  r[0][1] *= r_norm;
  r[1][1] = complex_t(r_norm, 0);

  complex_t dw = w[0][0] * w[1][1] - w[1][0] * w[0][1];
  complex_t dr = r[0][0] * r[1][1] - r[1][0] * r[0][1];

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      w[i][j] /= dw.sqrt_complex();
      r[i][j] /= dr.sqrt_complex();
    }
  }

  dw = w[0][0] * w[1][1] - w[1][0] * w[0][1];
  dr = r[0][0] * r[1][1] - r[1][0] * r[0][1];

  complex_t ax1[2][2];

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      ax1[i][j] = complex_t(0, 0);
      for (int k = 0; k < 2; k++) {
        ax1[i][j] += (w[k][i] % v[k][j]);
      }
    }
  }

  complex_t ax2[2][2];

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      ax2[i][j] = complex_t(0, 0);
      for (int k = 0; k < 2; k++) {
        ax2[i][j] += ax1[i][k] * r[k][j];
      }
    }
  }

  complex_t tmp = ax2[0][0] + ax2[1][1].conj();

  ax2[0][0] = tmp / tmp.module();
  ax2[1][1] = ax2[0][0].conj();
  ax2[0][1] = complex_t(0, 0);
  ax2[1][0] = complex_t(0, 0);

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      ax1[i][j] = complex_t(0, 0);
      for (int k = 0; k < 2; k++) {
        ax1[i][j] += w[i][k] * ax2[k][j];
      }
    }
  }

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      v[i][j] = complex_t(0, 0);
      for (int k = 0; k < 2; k++) {
        v[i][j] += ax1[i][k] ^ r[j][k];
      }
    }
  }
}

CUDA_DEVICE su3 su3::proj() {

  su3 A;

  complex_t tmp1(0, 0);
  double norm_factor = 0;

  // 1st orthonomalized vector
  for (int i = 0; i < 3; i++) {
    norm_factor += matrix[i][0].norm2();
  }
  norm_factor = 1 / sqrt(norm_factor);

  for (int i = 0; i < 3; i++) {
    A.matrix[i][0] = matrix[i][0] * norm_factor;
  }

  // 2nd orthonomalized vector
  tmp1 = complex_t(0, 0);
  for (int i = 0; i < 3; i++) {
    tmp1 -= matrix[i][1] ^ A.matrix[i][0];
  }

  for (int i = 0; i < 3; i++) {
    A.matrix[i][1] = matrix[i][1] + tmp1 * A.matrix[i][0];
  }

  norm_factor = 0;
  for (int i = 0; i < 3; i++) {
    norm_factor += A.matrix[i][1].norm2();
  }
  norm_factor = 1 / sqrt(norm_factor);

  for (int i = 0; i < 3; i++) {
    A.matrix[i][1] = A.matrix[i][1] * norm_factor;
  }

  // 3rd orthonomalized vector
  tmp1 = complex_t(0, 0);
  complex_t tmp2(0, 0);
  for (int i = 0; i < 3; i++) {
    tmp1 -= matrix[i][2] ^ A.matrix[i][0];
    tmp2 -= matrix[i][2] ^ A.matrix[i][1];
  }

  for (int i = 0; i < 3; i++) {
    A.matrix[i][2] =
        matrix[i][2] + tmp1 * A.matrix[i][0] + tmp2 * A.matrix[i][1];
  }

  norm_factor = 0;
  for (int i = 0; i < 3; i++) {
    norm_factor += A.matrix[i][2].norm2();
  }
  norm_factor = 1. / sqrt(norm_factor);

  for (int i = 0; i < 3; i++) {
    A.matrix[i][2] = A.matrix[i][2] * norm_factor;
  }

  // to make det = 1
  complex_t determinant = A.determinant();

  double ph = atan2(determinant.imag, determinant.real) / 3;

  tmp1 = complex_t(cos(ph), -sin(ph));

  A = A * tmp1;

  su3 v = A % *this;

  complex_t v1[2][2];

  for (int iter = 0; iter < 3; iter++) {

    // first subgroup

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        v1[i][j] = v.matrix[i][j];
      }
    }

    project_su2(v1);

    su3 s3;

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        s3.matrix[i][j] = v1[i][j];
      }
    }

    v = s3 % v;
    A = A * s3;

    // second subgroup

    v1[0][0] = v.matrix[0][0];
    v1[1][1] = v.matrix[2][2];
    v1[1][0] = v.matrix[2][0];
    v1[0][1] = v.matrix[0][2];

    project_su2(v1);

    s3 = su3();

    s3.matrix[0][0] = v1[0][0];
    s3.matrix[2][2] = v1[1][1];
    s3.matrix[2][0] = v1[1][0];
    s3.matrix[0][2] = v1[0][1];

    v = s3 % v;
    A = A * s3;

    // third subgroup

    v1[0][0] = v.matrix[1][1];
    v1[1][1] = v.matrix[2][2];
    v1[1][0] = v.matrix[2][1];
    v1[0][1] = v.matrix[1][2];

    project_su2(v1);

    s3 = su3();

    s3.matrix[1][1] = v1[0][0];
    s3.matrix[2][2] = v1[1][1];
    s3.matrix[2][1] = v1[1][0];
    s3.matrix[1][2] = v1[0][1];

    v = s3 % v;
    A = A * s3;
  }

  return A;
}

CUDA_DEVICE su3 su3::lambda3_mult() {

  su3 A;

  A.matrix[0][0] = matrix[0][0];
  A.matrix[0][1] = matrix[0][1].negative();
  A.matrix[1][0] = matrix[1][0].negative();
  A.matrix[1][1] = matrix[1][1];

  A.matrix[2][2] = complex_t(0, 0);

  return su3(A);
}

CUDA_DEVICE su3 su3::lambda8_mult() {

  su3 A(matrix);

  A.matrix[0][0] = A.matrix[0][0] / 3;
  A.matrix[0][1] = A.matrix[0][1] / 3;
  A.matrix[1][0] = A.matrix[1][0] / 3;
  A.matrix[1][1] = A.matrix[1][1] / 3;

  A.matrix[0][2] = A.matrix[0][2] * (-2. / 3);
  A.matrix[1][2] = A.matrix[1][2] * (-2. / 3);
  A.matrix[2][0] = A.matrix[2][0] * (-2. / 3);
  A.matrix[2][1] = A.matrix[2][1] * (-2. / 3);

  A.matrix[2][2] = A.matrix[2][2] * (4. / 3);

  return su3(A);
}

CUDA_DEVICE su3 operator+(const su3 &A, const su3 &B) {
  su3 C;
  complex_t a;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C.matrix[i][j] = A.matrix[i][j] + B.matrix[i][j];
    }
  }
  return C;
}
CUDA_DEVICE su3 operator-(const su3 &A, const su3 &B) {
  su3 C;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C.matrix[i][j] = A.matrix[i][j] - B.matrix[i][j];
    }
  }
  return C;
}
CUDA_DEVICE su3 operator*(const double &x, const su3 &A) {
  su3 C;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C.matrix[i][j] = x * A.matrix[i][j];
    }
  }
  return C;
}
CUDA_DEVICE su3 operator*(const su3 &A, const double &x) {
  su3 C;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C.matrix[i][j] = x * A.matrix[i][j];
    }
  }
  return C;
}
CUDA_DEVICE su3 operator*(const complex_t &x, const su3 &A) {
  su3 C;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C.matrix[i][j] = x * A.matrix[i][j];
    }
  }
  return C;
}
CUDA_DEVICE su3 operator*(const su3 &A, const complex_t &x) {
  su3 C;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C.matrix[i][j] = x * A.matrix[i][j];
    }
  }
  return C;
}

CUDA_DEVICE su3 operator*(const su3 &A, const su3 &B) {
  su3 C;
  complex_t a;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      a = complex_t(0, 0);
      for (int k = 0; k < 3; k++) {
        a = a + A.matrix[i][k] * B.matrix[k][j];
      }
      C.matrix[i][j] = a;
    }
  }
  return C;
}

CUDA_DEVICE su3 operator^(const su3 &A, const su3 &B) {
  su3 C;
  complex_t a;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      a = complex_t(0, 0);
      for (int k = 0; k < 3; k++) {
        a = a + (A.matrix[i][k] ^ B.matrix[j][k]);
      }
      C.matrix[i][j] = a;
    }
  }
  return C;
}

CUDA_DEVICE su3 operator%(const su3 &A, const su3 &B) {
  su3 C;
  complex_t a;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      a = complex_t(0, 0);
      for (int k = 0; k < 3; k++) {
        a = a + (A.matrix[k][i] % B.matrix[k][j]);
      }
      C.matrix[i][j] = a;
    }
  }
  return C;
}

// su3_abelian methods
CUDA_DEVICE su3_abelian::su3_abelian(complex_t B[3]) {
  for (int i = 0; i < 3; i++) {
    matrix[i] = B[i];
  }
}

CUDA_DEVICE su3_abelian::su3_abelian() {
  for (int i = 0; i < 3; i++) {
    matrix[i] = complex_t(1, 0);
  }
}

CUDA_DEVICE double su3_abelian::tr() {
  return (matrix[0].real + matrix[1].real + matrix[2].real) / 3;
}

CUDA_DEVICE complex_t su3_abelian::tr_complex() {
  return (matrix[0] + matrix[1] + matrix[2]) / 3;
}

CUDA_DEVICE double su3_abelian::multiply_tr(const su3_abelian &B) {
  double trace = 0;
  for (int i = 0; i < 3; i++) {
    trace +=
        matrix[i].real * B.matrix[i].real + matrix[i].imag * B.matrix[i].imag;
  }
  return trace / 3;
}

CUDA_DEVICE double su3_abelian::multiply_tr_adjoint(const su3_abelian &B) {
  double trace = 0;
  complex_t tmp = complex_t(0, 0);
  for (int i = 0; i < 3; i++) {
    tmp = matrix[i] ^ B.matrix[i];
    trace += tmp.real * tmp.real - tmp.imag * tmp.imag;
  }
  return trace;
}

CUDA_DEVICE su3_abelian su3_abelian::inverse() {
  su3_abelian B;
  for (int i = 0; i < 3; i++) {
    B.matrix[i].real = matrix[i].real;
    B.matrix[i].imag = -matrix[i].imag;
  }
  return B;
}

CUDA_DEVICE double su3_abelian::module() {
  complex_t tmp = matrix[0] * matrix[1];
  return tmp.real * matrix[2].real - tmp.imag * matrix[2].imag;
}

CUDA_DEVICE complex_t su3_abelian::determinant() {
  return matrix[0] * matrix[1] * matrix[2];
}

CUDA_DEVICE su3_abelian su3_abelian::conj() const {
  su3_abelian B;
  for (int i = 0; i < 3; i++) {
    B.matrix[i].real = matrix[i].real;
    B.matrix[i].imag = -matrix[i].imag;
  }
  return B;
}

CUDA_DEVICE su3_abelian su3_abelian::mult_by_imag(double x) {
  su3_abelian C;
  for (int i = 0; i < 3; i++) {
    C.matrix[i] = matrix[i].mult_by_imag(x);
  }
  return C;
}

CUDA_DEVICE su3_abelian su3_abelian::proj() {
  su3_abelian A;
  for (int i = 0; i < 3; i++) {
    A.matrix[i] = matrix[i] / matrix[i].module();
  }

  return A;
}

CUDA_DEVICE su3_abelian operator+(const su3_abelian &A, const su3_abelian &B) {
  su3_abelian C;
  complex_t a;
  for (int i = 0; i < 3; i++) {
    C.matrix[i] = A.matrix[i] + B.matrix[i];
  }
  return C;
}
CUDA_DEVICE su3_abelian operator-(const su3_abelian &A, const su3_abelian &B) {
  su3_abelian C;
  for (int i = 0; i < 3; i++) {
    C.matrix[i] = A.matrix[i] - B.matrix[i];
  }
  return C;
}
CUDA_DEVICE su3_abelian operator*(const double &x, const su3_abelian &A) {
  su3_abelian C;
  for (int i = 0; i < 3; i++) {
    C.matrix[i] = x * A.matrix[i];
  }
  return C;
}
CUDA_DEVICE su3_abelian operator*(const su3_abelian &A, const double &x) {
  su3_abelian C;
  for (int i = 0; i < 3; i++) {
    C.matrix[i] = x * A.matrix[i];
  }
  return C;
}
CUDA_DEVICE su3_abelian operator*(const complex_t &x, const su3_abelian &A) {
  su3_abelian C;
  for (int i = 0; i < 3; i++) {
    C.matrix[i] = x * A.matrix[i];
  }
  return C;
}
CUDA_DEVICE su3_abelian operator*(const su3_abelian &A, const complex_t &x) {
  su3_abelian C;
  for (int i = 0; i < 3; i++) {
    C.matrix[i] = x * A.matrix[i];
  }
  return C;
}

CUDA_DEVICE su3_abelian operator*(const su3_abelian &A, const su3_abelian &B) {
  su3_abelian C;

  C.matrix[0] = A.matrix[0] * B.matrix[0];
  C.matrix[1] = A.matrix[1] * B.matrix[1];
  C.matrix[2] = A.matrix[2] * B.matrix[2];

  return C;
}

CUDA_DEVICE su3_abelian operator^(const su3_abelian &A, const su3_abelian &B) {
  su3_abelian C;

  C.matrix[0] = A.matrix[0] ^ B.matrix[0];
  C.matrix[1] = A.matrix[1] ^ B.matrix[1];
  C.matrix[2] = A.matrix[2] ^ B.matrix[2];

  return C;
}

CUDA_DEVICE su3_abelian operator%(const su3_abelian &A, const su3_abelian &B) {
  su3_abelian C;

  C.matrix[0] = A.matrix[0] % B.matrix[0];
  C.matrix[1] = A.matrix[1] % B.matrix[1];
  C.matrix[2] = A.matrix[2] % B.matrix[2];

  return C;
}

// spin methods
CUDA_DEVICE spin::spin() {
  a1 = (double)0.0;
  a2 = (double)0.0;
  a3 = (double)1.0;
}
CUDA_DEVICE spin::spin(su2 U) {
  a1 = 2 * (U.a0 * U.a2 + U.a3 * U.a1);
  a2 = 2 * (U.a3 * U.a2 + U.a0 * U.a1);
  a3 = 1. - 2 * (U.a2 * U.a2 + U.a1 * U.a1);
}
CUDA_DEVICE spin::spin(double a1, double a2, double a3)
    : a1(a1), a2(a2), a3(a3) {}

CUDA_DEVICE double spin::norm() { return sqrt(a1 * a1 + a2 * a2 + a3 * a3); }

CUDA_DEVICE void spin::normalize() {
  double norm = sqrt(a1 * a1 + a2 * a2 + a3 * a3);
  a1 = a1 / norm;
  a2 = a2 / norm;
  a3 = a3 / norm;
}

// Function reflect spin vector through arbitrary spin vector
// it's more convenient to just change the spin variable
CUDA_DEVICE double spin::reflect(spin &V) {
  double tmp1 = (V.a1 * V.a1 - V.a2 * V.a2 - V.a3 * V.a3) * a1 +
                2 * V.a1 * (V.a2 * a2 + V.a3 * a3);

  double tmp2 = (-V.a1 * V.a1 + V.a2 * V.a2 - V.a3 * V.a3) * a2 +
                2 * V.a2 * (V.a1 * a1 + V.a3 * a3);

  double tmp3 = (-V.a1 * V.a1 - V.a2 * V.a2 + V.a3 * V.a3) * a3 +
                2 * V.a3 * (V.a1 * a1 + V.a2 * a2);

  double norm2 = V.norm() * V.norm();

  double b1 = tmp1 / norm2;
  double b2 = tmp2 / norm2;
  double b3 = tmp3 / norm2;

  double c = b1 * a1 + b2 * a2 + b3 * a3;

  a1 = tmp1 / norm2;
  a2 = tmp2 / norm2;
  a3 = tmp3 / norm2;

  return 1 - c;
}

CUDA_DEVICE void spin::reflect_fast(spin &V) {
  double tmp1 = (V.a1 * V.a1 - V.a2 * V.a2 - V.a3 * V.a3) * a1 +
                2 * V.a1 * (V.a2 * a2 + V.a3 * a3);

  double tmp2 = (-V.a1 * V.a1 + V.a2 * V.a2 - V.a3 * V.a3) * a2 +
                2 * V.a2 * (V.a1 * a1 + V.a3 * a3);

  double tmp3 = (-V.a1 * V.a1 - V.a2 * V.a2 + V.a3 * V.a3) * a3 +
                2 * V.a3 * (V.a1 * a1 + V.a2 * a2);

  double norm2 = V.norm() * V.norm();

  a1 = tmp1 / norm2;
  a2 = tmp2 / norm2;
  a3 = tmp3 / norm2;
}

CUDA_DEVICE double spin::parallel(spin &V) {
  double vNorm = V.norm();

  double b1 = V.a1 / vNorm;
  double b2 = V.a2 / vNorm;
  double b3 = V.a3 / vNorm;

  double c = b1 * a1 + b2 * a2 + b3 * a3;

  a1 = b1;
  a2 = b2;
  a3 = b3;

  return 1 - c;
}

CUDA_DEVICE void spin::parallel_fast(spin &V) {
  double vNorm = V.norm();

  a1 = V.a1 / vNorm;
  a2 = V.a2 / vNorm;
  a3 = V.a3 / vNorm;
}

CUDA_DEVICE bool spin::IsUnit() {
  return (this->norm() - 1.) > 1e-6 ? false : true;
}

// vector which spin variable is multiplied on
// it's contribution from single neighbour site in one direction
CUDA_DEVICE spin spin::contribution(const su2 &A) const {
  double q1 = A.a1 * A.a1;
  double q2 = A.a2 * A.a2;
  double q3 = A.a3 * A.a3;
  double q12 = 2. * A.a1 * A.a2;
  double q13 = 2. * A.a1 * A.a3;
  double q23 = 2. * A.a2 * A.a3;
  double q01 = 2. * A.a0 * A.a1;
  double q02 = 2. * A.a0 * A.a2;
  double q03 = 2. * A.a0 * A.a3;
  double A5 = A.a0 * A.a0 - q1 - q2 - q3;

  // matrix A(i, j) = tr(sigma(i) * U * sigma(j) * U.conj) is multiplied by
  // spin b variable from the right (A * b)
  return spin(a1 * (A5 + 2 * q1) + a2 * (q12 + q03) + a3 * (q13 - q02),
              a1 * (q12 - q03) + a2 * (A5 + 2 * q2) + a3 * (q23 + q01),
              a1 * (q13 + q02) + a2 * (q23 - q01) + a3 * (A5 + 2 * q3));
}

// supposed to be used for spin variable on current lattice cite
CUDA_DEVICE spin spin::contribution_conj(const su2 &A) const {
  double q1 = A.a1 * A.a1;
  double q2 = A.a2 * A.a2;
  double q3 = A.a3 * A.a3;
  double q12 = 2. * A.a1 * A.a2;
  double q13 = 2. * A.a1 * A.a3;
  double q23 = 2. * A.a2 * A.a3;
  double q01 = 2. * A.a0 * A.a1;
  double q02 = 2. * A.a0 * A.a2;
  double q03 = 2. * A.a0 * A.a3;
  double A5 = A.a0 * A.a0 - q1 - q2 - q3;

  // matrix A(i, j) = tr(sigma(i) * U * sigma(j) * U.conj) is multiplied by
  //  spin b variable from the left (b * A)
  return spin(a1 * (A5 + 2 * q1) + a2 * (q12 - q03) + a3 * (q13 + q02),
              a1 * (q12 + q03) + a2 * (A5 + 2 * q2) + a3 * (q23 - q01),
              a1 * (q13 - q02) + a2 * (q23 + q01) + a3 * (A5 + 2 * q3));
}

CUDA_DEVICE void spin::contribution1(const su2 &A, const spin &b) {
  double q1 = A.a1 * A.a1;
  double q2 = A.a2 * A.a2;
  double q3 = A.a3 * A.a3;
  double q12 = A.a1 * A.a2;
  double q13 = A.a1 * A.a3;
  double q23 = A.a2 * A.a3;
  double q01 = A.a0 * A.a1;
  double q02 = A.a0 * A.a2;
  double q03 = A.a0 * A.a3;
  double A5 = A.a0 * A.a0 - q1 - q2 - q3;

  a1 += b.a1 * (A5 + 2 * q1) + 2 * b.a2 * (q12 + q03) + 2 * b.a3 * (q13 - q02);
  a2 += 2 * b.a1 * (q12 - q03) + b.a2 * (A5 + 2 * q2) + 2 * b.a3 * (q23 + q01);
  a3 += 2 * b.a1 * (q13 + q02) + 2 * b.a2 * (q23 - q01) + b.a3 * (A5 + 2 * q3);
}

CUDA_DEVICE void spin::contribution1_conj(const su2 &A, const spin &b) {
  double q1 = A.a1 * A.a1;
  double q2 = A.a2 * A.a2;
  double q3 = A.a3 * A.a3;
  double q12 = A.a1 * A.a2;
  double q13 = A.a1 * A.a3;
  double q23 = A.a2 * A.a3;
  double q01 = A.a0 * A.a1;
  double q02 = A.a0 * A.a2;
  double q03 = A.a0 * A.a3;
  double A5 = A.a0 * A.a0 - q1 - q2 - q3;

  a1 += b.a1 * (A5 + 2 * q1) + 2 * b.a2 * (q12 - q03) + 2 * b.a3 * (q13 + q02);
  a2 += 2 * b.a1 * (q12 + q03) + b.a2 * (A5 + 2 * q2) + 2 * b.a3 * (q23 - q01);
  a3 += 2 * b.a1 * (q13 - q02) + 2 * b.a2 * (q23 + q01) + b.a3 * (A5 + 2 * q3);
}

// Calculation have been done for specific element of gauge matrix G with g_3 =
// 0 Here:
//      | g_0 + I g_1    g_2 + I g_3 |
// G =  |                            |
//      |-g_2 + I g_3    g_0 - I g_1 |
//
// spinVector = (a, b, c) = (a1, a2, a3)
//
// g_0 = - a/sqrt(2-2*c)
// g_1 = - b/sqrt(2-2*c)
// g_2 = - sqrt(2-2*c) / 2
// g_3 = 0
CUDA_DEVICE su2 spin::GetGaugeMatrix() {
  double sq = sqrt(2. - 2 * a3);
  return su2(-a1 / sq, 0., -sq / 2., -a2 / sq);
}
CUDA_DEVICE spin operator*(const double &x, const spin &A) {
  return spin(A.a1 * x, A.a2 * x, A.a3 * x);
}
CUDA_DEVICE spin operator*(const spin &A, const double &x) {
  return spin(A.a1 * x, A.a2 * x, A.a3 * x);
}
CUDA_DEVICE spin operator+(const spin &A, const spin &B) {
  return spin(A.a1 + B.a1, A.a2 + B.a2, A.a3 + B.a3);
}
CUDA_DEVICE spin operator-(const spin &A, const spin &B) {
  return spin(A.a1 - B.a1, A.a2 - B.a2, A.a3 - B.a3);
}
CUDA_DEVICE double operator*(const spin &A, const spin &B) {
  return A.a1 * B.a1 + A.a2 * B.a2 + A.a3 * B.a3;
}

} // namespace matrix_gpu