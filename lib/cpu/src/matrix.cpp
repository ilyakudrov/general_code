#include "../include/matrix.h"
#include <cmath>

// su2 methods
su2::su2() {
  a0 = 1;
  a1 = 0;
  a2 = 0;
  a3 = 0;
}

su2::su2(FLOAT b0, FLOAT b1, FLOAT b2, FLOAT b3) {
  a0 = b0;
  a1 = b1;
  a2 = b2;
  a3 = b3;
}

FLOAT su2::tr() { return 2 * a0; }

su2 su2::inverse() {
  double rho = a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3;
  return su2(a0 / rho, -a1 / rho, -a2 / rho, -a3 / rho);
}
FLOAT su2::module() { return a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3; }
su2 su2::conj() const { return su2(a0, -a1, -a2, -a3); }
su2 su2::proj() {
  FLOAT rho = a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3;
  return su2(a0 / powf(rho, 0.5), a1 / powf(rho, 0.5), a2 / powf(rho, 0.5),
             a3 / powf(rho, 0.5));
}
su2 su2::sigma3_mult() const { return su2(a0, -a1, -a2, a3); }

su2 operator+(const su2 &A, const su2 &B) {
  return su2(A.a0 + B.a0, A.a1 + B.a1, A.a2 + B.a2, A.a3 + B.a3);
};
su2 operator-(const su2 &A, const su2 &B) {
  return su2(A.a0 - B.a0, A.a1 - B.a1, A.a2 - B.a2, A.a3 - B.a3);
};
su2 operator*(const FLOAT &x, const su2 &A) {
  return su2(A.a0 * x, A.a1 * x, A.a2 * x, A.a3 * x);
};
su2 operator*(const su2 &A, const FLOAT &x) {
  return su2(A.a0 * x, A.a1 * x, A.a2 * x, A.a3 * x);
};
su2 operator*(const su2 &A, const su2 &B) {
  return su2(A.a0 * B.a0 - A.a1 * B.a1 - A.a2 * B.a2 - A.a3 * B.a3,
             A.a0 * B.a1 + B.a0 * A.a1 + A.a3 * B.a2 - A.a2 * B.a3,
             A.a0 * B.a2 + B.a0 * A.a2 + A.a1 * B.a3 - A.a3 * B.a1,
             A.a0 * B.a3 + B.a0 * A.a3 + A.a2 * B.a1 - A.a1 * B.a2);
};

ostream &operator<<(ostream &os, const su2 &A) {
  os << "a0 = " << A.a0 << " "
     << "a1 = " << A.a1 << " "
     << "a2 = " << A.a2 << " "
     << "a3 = " << A.a3;
  return os;
}

// abelian methods
abelian::abelian() {
  r = 1;
  phi = 0;
}

abelian::abelian(FLOAT r1, FLOAT phi1) {
  r = r1;
  phi = phi1;
}

FLOAT abelian::tr() { return cos(phi); }

abelian abelian::inverse() { return abelian(1 / r, -phi); }
FLOAT abelian::module() { return r; }
abelian abelian::conj() const { return abelian(r, -phi); }
abelian abelian::proj() { return abelian(1, phi); }

abelian operator+(const abelian &A, const abelian &B) {
  return abelian(sqrt((A.r * sin(A.phi) + B.r * sin(B.phi)) *
                          (A.r * sin(A.phi) + B.r * sin(B.phi)) +
                      (A.r * cos(A.phi) + B.r * cos(B.phi)) *
                          (A.r * cos(A.phi) + B.r * cos(B.phi))),
                 atan2(A.r * sin(A.phi) + B.r * sin(B.phi),
                       A.r * cos(A.phi) + B.r * cos(B.phi)));
};
abelian operator-(const abelian &A, const abelian &B) {
  return abelian(sqrt((A.r * sin(A.phi) - B.r * sin(B.phi)) *
                          (A.r * sin(A.phi) - B.r * sin(B.phi)) +
                      (A.r * cos(A.phi) - B.r * cos(B.phi)) *
                          (A.r * cos(A.phi) - B.r * cos(B.phi))),
                 atan2(A.r * sin(A.phi) - B.r * sin(B.phi),
                       A.r * cos(A.phi) - B.r * cos(B.phi)));
};
abelian operator*(const FLOAT &x, const abelian &A) {
  return abelian(A.r * x, A.phi);
};
abelian operator*(const abelian &A, const FLOAT &x) {
  return abelian(A.r * x, A.phi);
};
abelian operator*(const abelian &A, const abelian &B) {
  return abelian(A.r * B.r, A.phi + B.phi);
};

ostream &operator<<(ostream &os, const abelian &A) {
  os << "r = " << A.r << " "
     << "phi = " << A.phi << " " << endl;
  return os;
}