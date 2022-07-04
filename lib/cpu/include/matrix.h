#pragma once

#include <complex>
#include <iostream>

struct complex_t {
  double real;
  double imag;

  complex_t(const double real1, const double imag1);
  complex_t();

  double module();

  double angle();

  complex_t mult_by_imag(double x);
};

complex_t operator+(const complex_t &a, const complex_t &b);

complex_t operator-(const complex_t &a, const complex_t &b);

complex_t operator*(const complex_t &a, const complex_t &b);

complex_t operator*(const double &a, const complex_t &b);

complex_t operator*(const complex_t &a, const double &b);

complex_t operator^(const complex_t &a, const complex_t &b);

complex_t operator%(const complex_t &a, const complex_t &b);

complex_t operator&(const complex_t &a, const complex_t &b);

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

  double multiply_tr(const su2 &B);

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

// matrix multiplication A * B.conj()
su2 operator^(const su2 &A, const su2 &B);

// matrix multiplication A.conj() * B
su2 operator%(const su2 &A, const su2 &B);

std::ostream &operator<<(std::ostream &os, const su2 &A);

// abelian variable (module * exp(i * phi))
class abelian {
public:
  double r, phi;
  abelian();
  abelian(double r1, double phi1);

  // trace
  double tr();

  double multiply_tr(const abelian &B);

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
abelian operator^(const abelian &A, const abelian &B);
abelian operator%(const abelian &A, const abelian &B);

std::ostream &operator<<(std::ostream &os, const abelian &A);

// su3 matrix in 3x3 complex matrix representation
class su3 {
public:
  complex_t matrix[3][3];
  su3(complex_t B[3][3]);
  su3();

  // calculate trace of the matrix
  double tr();

  double multiply_tr(const su3 &B);

  // calculate inverse of the matrix
  su3 inverse();

  // calculate conjugate of the matrix
  su3 conj() const;

  su3 mult_by_imag(double x);

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

su3 operator*(const complex_t &x, const su3 &A);
su3 operator*(const su3 &A, const complex_t &x);

// matrix multiplication A * B
su3 operator*(const su3 &A, const su3 &B);

// matrix multiplication A * B.conj()
su3 operator^(const su3 &A, const su3 &B);

// matrix multiplication A.conj() * B
su3 operator%(const su3 &A, const su3 &B);

std::ostream &operator<<(std::ostream &os, const su3 &A);

// 3D vector realisation for spin model.
class spin {
public:
  double a1, a2, a3;

  spin(double a1, double a2, double a3);
  spin(su2 U);
  spin();

  double norm();
  // set norm equal to unity
  void normalize();
  // reflect spin-vector through the V-vector axis
  double reflect(spin &V);
  // reflect without difference calculation
  void reflect_fast(spin &V);
  // rorate spin-vector to the same direction as V-vector
  double parallel(spin &V);
  // parallel without difference calculation
  void parallel_fast(spin &V);
  // check identity matrix
  bool IsUnit();
  // contribution from neighbour site
  spin contribution(const su2 &A) const;
  // contribution from neighbour site in negative direction (conjugated matrix)
  spin contribution_conj(const su2 &A) const;
  void contribution1(const su2 &A, const spin &b);
  void contribution1_conj(const su2 &A, const spin &b);
  // calculation of gauge matrix G from spin vector.
  su2 GetGaugeMatrix();
};

spin operator+(const spin &A, const spin &B);
spin operator-(const spin &A, const spin &B);
spin operator*(const double &x, const spin &A);
spin operator*(const spin &A, const double &x);
double operator*(const spin &A, const spin &B);

std::ostream &operator<<(std::ostream &os, const spin &A);