#pragma once

#ifdef DOUBLE
#define FLOAT double
#else
#define FLOAT float
#endif

#include <iostream>

// su2 matrix in sigma matrices representation (a0 + ai * sigma[i])
class su2 {
public:
  FLOAT a0, a1, a2, a3;
  su2(FLOAT b0, FLOAT b1, FLOAT b2, FLOAT b3);
  su2();

  // calculate trace of the matrix
  FLOAT tr();

  // calculate inverse of the matrix
  su2 inverse();

  // calculate conjugate of the matrix
  su2 conj() const;

  // gets projection onto su2 group
  su2 proj();

  // calculates module of vector in sigma matrices representation
  FLOAT module();

  // left and right multiplication of the matrix on sigma[3] matrix
  su2 sigma3_mult() const;
};

su2 operator+(const su2 &A, const su2 &B);
su2 operator-(const su2 &A, const su2 &B);
su2 operator*(const FLOAT &x, const su2 &A);
su2 operator*(const su2 &A, const FLOAT &x);

// matrix multiplication A * B
su2 operator*(const su2 &A, const su2 &B);

// matrix multiplication A * B
su2 operator*(const su2 &A, const su2 *B);

// matrix multiplication A * B.conj()
su2 operator^(const su2 &A, const su2 *B);

std::ostream &operator<<(std::ostream &os, const su2 &A);

// abelian variable (module * exp(i * phi))
class abelian {
public:
  FLOAT r, phi;
  abelian();
  abelian(FLOAT r1, FLOAT phi1);

  // trace
  FLOAT tr();

  // inverse
  abelian inverse();
  // conjugated
  abelian conj() const;
  abelian proj();
  FLOAT module();
};

abelian operator+(const abelian &A, const abelian &B);
abelian operator-(const abelian &A, const abelian &B);
abelian operator*(const FLOAT &x, const abelian &A);
abelian operator*(const abelian &A, const FLOAT &x);

abelian operator*(const abelian &A, const abelian &B);
abelian operator*(const abelian &A, const abelian *B);
abelian operator^(const abelian &A, const abelian *B);

std::ostream &operator<<(std::ostream &os, const abelian &A);