#pragma once

#include <complex>
#include <iostream>

struct complex_t {
  double real;
  double imag;

  complex_t(const double real1, const double imag1);
  complex_t();
};

complex_t operator+(const complex_t &a, const complex_t &b);

complex_t operator-(const complex_t &a, const complex_t &b);

complex_t operator*(const complex_t &a, const complex_t &b);

complex_t operator*(const double &a, const complex_t &b);

complex_t operator*(const complex_t &a, const double &b);

complex_t operator^(const complex_t &a, const complex_t &b);

complex_t operator%(const complex_t &a, const complex_t &b);

complex_t operator/(const complex_t &a, const double &b);

std::ostream &operator<<(std::ostream &os, const complex_t &a);

// su2 matrix in sigma matrices representation (a0 + ai * sigma[i])
class su2 {
public:
  double a0, a1, a2, a3;
  su2(double b0, double b1, double b2, double b3);
  su2();

  // calculate trace of the matrix
  double tr();

  double multiply_tr(const su2 *B);

  // calculate inverse of the matrix
  su2 inverse();

  // calculate conjugate of the matrix
  su2 conj() const;

  // gets projection onto su2 group
  su2 proj();

  // calculates module of vector in sigma matrices representation
  double module();

  // left and right multiplication of the matrix on sigma[3] matrix
  su2 sigma3_mult() const;
};

su2 operator+(const su2 &A, const su2 &B);
su2 operator-(const su2 &A, const su2 &B);
su2 operator*(const double &x, const su2 &A);
su2 operator*(const su2 &A, const double &x);

// matrix multiplication A * B
su2 operator*(const su2 &A, const su2 &B);

// matrix multiplication A * B
su2 operator*(const su2 &A, const su2 *B);

// matrix multiplication A * B.conj()
su2 operator^(const su2 &A, const su2 *B);

// matrix multiplication A.conj() * B
su2 operator%(const su2 &A, const su2 *B);

std::ostream &operator<<(std::ostream &os, const su2 &A);

// abelian variable (module * exp(i * phi))
class abelian {
public:
  double r, phi;
  abelian();
  abelian(double r1, double phi1);

  // trace
  double tr();

  double multiply_tr(const abelian *B);

  // inverse
  abelian inverse();
  // conjugated
  abelian conj() const;
  abelian proj();
  double module();
};

abelian operator+(const abelian &A, const abelian &B);
abelian operator-(const abelian &A, const abelian &B);
abelian operator*(const double &x, const abelian &A);
abelian operator*(const abelian &A, const double &x);

abelian operator*(const abelian &A, const abelian &B);
abelian operator*(const abelian &A, const abelian *B);
abelian operator^(const abelian &A, const abelian *B);
abelian operator%(const abelian &A, const abelian *B);

std::ostream &operator<<(std::ostream &os, const abelian &A);

// su3 matrix in 3x3 complex matrix representation
class su3_full {
public:
  std::complex<double> matrix[3][3];
  su3_full(std::complex<double> B[3][3]);
  su3_full();

  // calculate trace of the matrix
  double tr();

  double multiply_tr(const su3_full *B);

  // calculate inverse of the matrix
  su3_full inverse();

  // calculate conjugate of the matrix
  su3_full conj() const;

  // gets projection onto su3 group
  su3_full proj();

  // calculates module of vector in sigma matrices representation
  double module();

  std::complex<double> determinant();

  std::complex<double> unitarity_check();
};

su3_full operator+(const su3_full &A, const su3_full &B);
su3_full operator-(const su3_full &A, const su3_full &B);
su3_full operator*(const double &x, const su3_full &A);
su3_full operator*(const su3_full &A, const double &x);

// matrix multiplication A * B
su3_full operator*(const su3_full &A, const su3_full &B);

// matrix multiplication A * B
su3_full operator*(const su3_full &A, const su3_full *B);

// matrix multiplication A * B.conj()
su3_full operator^(const su3_full &A, const su3_full *B);

// matrix multiplication A.conj() * B
su3_full operator%(const su3_full &A, const su3_full *B);

std::ostream &operator<<(std::ostream &os, const su3_full &A);

// su3 matrix in 3x3 complex matrix representation
class su3 {
public:
  complex_t matrix[3][3];
  su3(complex_t B[3][3]);
  su3();

  // calculate trace of the matrix
  double tr();

  double multiply_tr(const su3 *B);

  // calculate inverse of the matrix
  su3 inverse();

  // calculate conjugate of the matrix
  su3 conj() const;

  // gets projection onto su3 group
  su3 proj();

  // calculates module of vector in sigma matrices representation
  double module();

  complex_t determinant();

  complex_t unitarity_check();
};

su3 operator+(const su3 &A, const su3 &B);
su3 operator-(const su3 &A, const su3 &B);
su3 operator*(const double &x, const su3 &A);
su3 operator*(const su3 &A, const double &x);

// matrix multiplication A * B
su3 operator*(const su3 &A, const su3 &B);

// matrix multiplication A * B
su3 operator*(const su3 &A, const su3 *B);

// matrix multiplication A * B.conj()
su3 operator^(const su3 &A, const su3 *B);

// matrix multiplication A.conj() * B
su3 operator%(const su3 &A, const su3 *B);

std::ostream &operator<<(std::ostream &os, const su3 &A);