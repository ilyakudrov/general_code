#pragma once

#ifdef DOUBLE
#define FLOAT double
#else
#define FLOAT float
#endif

#include <iostream>

using namespace std;

class su2 {
public:
  FLOAT a0, a1, a2, a3;
  su2(FLOAT b0, FLOAT b1, FLOAT b2, FLOAT b3);
  su2();

  FLOAT tr();
  su2 inverse();
  su2 conj();
  su2 proj();
  FLOAT module();
};

su2 operator+(const su2 &A, const su2 &B);
su2 operator-(const su2 &A, const su2 &B);
su2 operator*(const FLOAT &x, const su2 &A);
su2 operator*(const su2 &A, const FLOAT &x);
su2 operator*(const su2 &A, const su2 &B);

ostream &operator<<(ostream &os, const su2 &A);

class abelian {
public:
  FLOAT r, phi;
  abelian();
  abelian(FLOAT r1, FLOAT phi1);

  FLOAT tr();
  abelian inverse();
  abelian conj();
  abelian proj();
  FLOAT module();
};

abelian operator+(const abelian &A, const abelian &B);
abelian operator-(const abelian &A, const abelian &B);
abelian operator*(const FLOAT &x, const abelian &A);
abelian operator*(const abelian &A, const FLOAT &x);
abelian operator*(const abelian &A, const abelian &B);

ostream &operator<<(ostream &os, const abelian &A);
