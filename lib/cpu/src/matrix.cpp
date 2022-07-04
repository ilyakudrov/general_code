#include "../include/matrix.h"
#include <cmath>

complex_t::complex_t(const double real1, const double imag1) {
  real = real1;
  imag = imag1;
}

complex_t::complex_t() {
  real = 0;
  imag = 0;
}

double complex_t::module() { return sqrt(real * real + imag * imag); }

double complex_t::angle() { return atan2(imag, real); }

complex_t complex_t::mult_by_imag(double x) {
  return complex_t(-imag * x, real * x);
}

complex_t operator+(const complex_t &a, const complex_t &b) {
  return complex_t(a.real + b.real, a.imag + b.imag);
}

complex_t operator-(const complex_t &a, const complex_t &b) {
  return complex_t(a.real - b.real, a.imag - b.imag);
}

complex_t operator*(const complex_t &a, const complex_t &b) {
  return complex_t(a.real * b.real - a.imag * b.imag,
                   a.real * b.imag + a.imag * b.real);
}

complex_t operator*(const double &a, const complex_t &b) {
  return complex_t(a * b.real, a * b.imag);
}

complex_t operator*(const complex_t &a, const double &b) {
  return complex_t(a.real * b, a.imag * b);
}

complex_t operator^(const complex_t &a, const complex_t &b) {
  return complex_t(a.real * b.real + a.imag * b.imag,
                   a.imag * b.real - a.real * b.imag);
}

complex_t operator%(const complex_t &a, const complex_t &b) {
  return complex_t(a.real * b.real + a.imag * b.imag,
                   a.real * b.imag - a.imag * b.real);
}

complex_t operator&(const complex_t &a, const complex_t &b) {
  return complex_t(a.real * b.real - a.imag * b.imag,
                   -a.imag * b.real - a.real * b.imag);
}

complex_t operator/(const complex_t &a, const double &b) {
  return complex_t(a.real / b, a.imag / b);
}

std::ostream &operator<<(std::ostream &os, const complex_t &a) {
  os << "(" << a.real << ", " << a.imag << ")";
  return os;
}

// su2 methods
su2::su2() {
  a0 = 1.;
  a1 = 0;
  a2 = 0;
  a3 = 0;
}

su2::su2(double b0, double b1, double b2, double b3) {
  a0 = b0;
  a1 = b1;
  a2 = b2;
  a3 = b3;
}

double su2::tr() { return a0; }

double su2::multiply_tr(const su2 &B) {
  return a0 * B.a0 + a1 * B.a1 + a2 * B.a2 + a3 * B.a3;
}

su2 su2::inverse() {
  double rho = a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3;
  return su2(a0 / rho, -a1 / rho, -a2 / rho, -a3 / rho);
}
double su2::module() { return a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3; }
su2 su2::conj() const { return su2(a0, -a1, -a2, -a3); }
su2 su2::proj() {
  double rho = sqrt(a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3);
  return su2(a0 / rho, a1 / rho, a2 / rho, a3 / rho);
}
su2 su2::sigma3_mult() const { return su2(a0, -a1, -a2, a3); }

su2 operator+(const su2 &A, const su2 &B) {
  return su2(A.a0 + B.a0, A.a1 + B.a1, A.a2 + B.a2, A.a3 + B.a3);
}
su2 operator-(const su2 &A, const su2 &B) {
  return su2(A.a0 - B.a0, A.a1 - B.a1, A.a2 - B.a2, A.a3 - B.a3);
}
su2 operator*(const double &x, const su2 &A) {
  return su2(A.a0 * x, A.a1 * x, A.a2 * x, A.a3 * x);
}
su2 operator*(const su2 &A, const double &x) {
  return su2(A.a0 * x, A.a1 * x, A.a2 * x, A.a3 * x);
}

su2 operator*(const su2 &A, const su2 &B) {
  return su2(A.a0 * B.a0 - A.a1 * B.a1 - A.a2 * B.a2 - A.a3 * B.a3,
             A.a0 * B.a1 + B.a0 * A.a1 + A.a3 * B.a2 - A.a2 * B.a3,
             A.a0 * B.a2 + B.a0 * A.a2 + A.a1 * B.a3 - A.a3 * B.a1,
             A.a0 * B.a3 + B.a0 * A.a3 + A.a2 * B.a1 - A.a1 * B.a2);
}
su2 operator^(const su2 &A, const su2 &B) {
  return su2(A.a0 * B.a0 + A.a1 * B.a1 + A.a2 * B.a2 + A.a3 * B.a3,
             -A.a0 * B.a1 + B.a0 * A.a1 - A.a3 * B.a2 + A.a2 * B.a3,
             -A.a0 * B.a2 + B.a0 * A.a2 - A.a1 * B.a3 + A.a3 * B.a1,
             -A.a0 * B.a3 + B.a0 * A.a3 - A.a2 * B.a1 + A.a1 * B.a2);
}
su2 operator%(const su2 &A, const su2 &B) {
  return su2(A.a0 * B.a0 + A.a1 * B.a1 + A.a2 * B.a2 + A.a3 * B.a3,
             A.a0 * B.a1 - B.a0 * A.a1 - A.a3 * B.a2 + A.a2 * B.a3,
             A.a0 * B.a2 - B.a0 * A.a2 - A.a1 * B.a3 + A.a3 * B.a1,
             A.a0 * B.a3 - B.a0 * A.a3 - A.a2 * B.a1 + A.a1 * B.a2);
}

std::ostream &operator<<(std::ostream &os, const su2 &A) {
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

abelian::abelian(double r1, double phi1) {
  r = r1;
  phi = phi1;
}

double abelian::tr() { return r * cos(phi); }

double abelian::multiply_tr(const abelian &B) {
  return r * B.r * cos(phi - B.phi);
}

abelian abelian::inverse() { return abelian(1 / r, -phi); }
double abelian::module() { return r; }
abelian abelian::conj() const { return abelian(r, -phi); }
abelian abelian::proj() { return abelian(1, phi); }

abelian operator+(const abelian &A, const abelian &B) {
  return abelian(sqrt((A.r * sin(A.phi) + B.r * sin(B.phi)) *
                          (A.r * sin(A.phi) + B.r * sin(B.phi)) +
                      (A.r * cos(A.phi) + B.r * cos(B.phi)) *
                          (A.r * cos(A.phi) + B.r * cos(B.phi))),
                 atan2(A.r * sin(A.phi) + B.r * sin(B.phi),
                       A.r * cos(A.phi) + B.r * cos(B.phi)));
}
abelian operator-(const abelian &A, const abelian &B) {
  return abelian(sqrt((A.r * sin(A.phi) - B.r * sin(B.phi)) *
                          (A.r * sin(A.phi) - B.r * sin(B.phi)) +
                      (A.r * cos(A.phi) - B.r * cos(B.phi)) *
                          (A.r * cos(A.phi) - B.r * cos(B.phi))),
                 atan2(A.r * sin(A.phi) - B.r * sin(B.phi),
                       A.r * cos(A.phi) - B.r * cos(B.phi)));
}
abelian operator*(const double &x, const abelian &A) {
  return abelian(A.r * x, A.phi);
}
abelian operator*(const abelian &A, const double &x) {
  return abelian(A.r * x, A.phi);
}

abelian operator*(const abelian &A, const abelian &B) {
  return abelian(A.r * B.r, A.phi + B.phi);
}
abelian operator^(const abelian &A, const abelian &B) {
  return abelian(A.r * B.r, A.phi - B.phi);
}
abelian operator%(const abelian &A, const abelian &B) {
  return abelian(A.r * B.r, B.phi - A.phi);
}

std::ostream &operator<<(std::ostream &os, const abelian &A) {
  os << "r = " << A.r << " "
     << "phi = " << A.phi << " " << std::endl;
  return os;
}

// su3 methods
su3::su3() {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i != j)
        matrix[i][j] = complex_t(0, 0);
      else
        matrix[i][j] = complex_t(1., 0);
    }
  }
}

su3::su3(complex_t B[3][3]) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      matrix[i][j] = B[i][j];
    }
  }
}

double su3::tr() {
  return (matrix[0][0].real + matrix[1][1].real + matrix[2][2].real) / 3;
}

double su3::multiply_tr(const su3 &B) {
  double trace = 0;
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++) {
      trace += matrix[i][k].real * B.matrix[i][k].real +
               matrix[i][k].imag * B.matrix[i][k].imag;
    }
  }
  return trace / 3;
}

su3 su3::inverse() {
  su3 B;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      B.matrix[i][j].real = matrix[j][i].real;
      B.matrix[i][j].imag = -matrix[j][i].imag;
    }
  }
  return B;
}

complex_t determinant_part(complex_t B[3][3], int i, int j, int k) {
  return B[0][i] * (B[1][j] * B[2][k] - B[2][j] * B[1][k]);
}

double su3::module() {
  complex_t determinant;
  determinant = determinant_part(matrix, 0, 1, 2);
  determinant = determinant + determinant_part(matrix, 1, 0, 2);
  determinant = determinant + determinant_part(matrix, 2, 0, 1);
  return determinant.real;
}

complex_t su3::determinant() {
  complex_t determinant;
  determinant = determinant_part(matrix, 0, 1, 2);
  determinant = determinant - determinant_part(matrix, 1, 0, 2);
  determinant = determinant - determinant_part(matrix, 2, 0, 1);
  return determinant;
}

complex_t su3::unitarity_check() {
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

su3 su3::conj() const {
  su3 B;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      B.matrix[i][j].real = matrix[j][i].real;
      B.matrix[i][j].imag = -matrix[j][i].imag;
    }
  }
  return B;
}

su3 su3::mult_by_imag(double x) {
  su3 C;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C.matrix[i][j] = matrix[i][j].mult_by_imag(x);
    }
  }
  return C;
}

su3 su3::proj() {
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
}

su3 operator+(const su3 &A, const su3 &B) {
  su3 C;
  complex_t a;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C.matrix[i][j] = A.matrix[i][j] + B.matrix[i][j];
    }
  }
  return C;
}
su3 operator-(const su3 &A, const su3 &B) {
  su3 C;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C.matrix[i][j] = A.matrix[i][j] - B.matrix[i][j];
    }
  }
  return C;
}
su3 operator*(const double &x, const su3 &A) {
  su3 C;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C.matrix[i][j] = x * A.matrix[i][j];
    }
  }
  return C;
}
su3 operator*(const su3 &A, const double &x) {
  su3 C;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C.matrix[i][j] = x * A.matrix[i][j];
    }
  }
  return C;
}
su3 operator*(const complex_t &x, const su3 &A) {
  su3 C;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C.matrix[i][j] = x * A.matrix[i][j];
    }
  }
  return C;
}
su3 operator*(const su3 &A, const complex_t &x) {
  su3 C;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      C.matrix[i][j] = x * A.matrix[i][j];
    }
  }
  return C;
}

su3 operator*(const su3 &A, const su3 &B) {
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

su3 operator^(const su3 &A, const su3 &B) {
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

su3 operator%(const su3 &A, const su3 &B) {
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

std::ostream &operator<<(std::ostream &os, const su3 &A) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      os << "(" << A.matrix[i][j].real << ", " << A.matrix[i][j].imag << ") ";
    }
    os << std::endl;
  }
  return os;
}

spin::spin() {
  a1 = (double)0.0;
  a2 = (double)0.0;
  a3 = (double)1.0;
}
spin::spin(su2 U) {
  a1 = 2 * (U.a0 * U.a2 + U.a3 * U.a1);
  a2 = 2 * (U.a3 * U.a2 + U.a0 * U.a1);
  a3 = 1. - 2 * (U.a2 * U.a2 + U.a1 * U.a1);
}
spin::spin(double a1, double a2, double a3) : a1(a1), a2(a2), a3(a3) {}

double spin::norm() { return sqrt(a1 * a1 + a2 * a2 + a3 * a3); }

void spin::normalize() {
  double norm = sqrt(a1 * a1 + a2 * a2 + a3 * a3);
  a1 = a1 / norm;
  a2 = a2 / norm;
  a3 = a3 / norm;
}

// Function reflect spin vector through arbitrary spin vector
// it's more convenient to just change the spin variable
double spin::reflect(spin &V) {
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

void spin::reflect_fast(spin &V) {
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

double spin::parallel(spin &V) {
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

void spin::parallel_fast(spin &V) {
  double vNorm = V.norm();

  a1 = V.a1 / vNorm;
  a2 = V.a2 / vNorm;
  a3 = V.a3 / vNorm;
}

bool spin::IsUnit() { return (this->norm() - 1.) > 1e-6 ? false : true; }

// vector which spin variable is multiplied on
// it's contribution from single neighbour site in one direction
spin spin::contribution(const su2 &A) const {
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
spin spin::contribution_conj(const su2 &A) const {
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

void spin::contribution1(const su2 &A, const spin &b) {
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

void spin::contribution1_conj(const su2 &A, const spin &b) {
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
su2 spin::GetGaugeMatrix() {
  double sq = sqrt(2. - 2 * a3);
  return su2(-a1 / sq, 0., -sq / 2., -a2 / sq);
}
spin operator*(const double &x, const spin &A) {
  return spin(A.a1 * x, A.a2 * x, A.a3 * x);
}
spin operator*(const spin &A, const double &x) {
  return spin(A.a1 * x, A.a2 * x, A.a3 * x);
}
spin operator+(const spin &A, const spin &B) {
  return spin(A.a1 + B.a1, A.a2 + B.a2, A.a3 + B.a3);
}
spin operator-(const spin &A, const spin &B) {
  return spin(A.a1 - B.a1, A.a2 - B.a2, A.a3 - B.a3);
}
double operator*(const spin &A, const spin &B) {
  return A.a1 * B.a1 + A.a2 * B.a2 + A.a3 * B.a3;
}

std::ostream &operator<<(std::ostream &os, const spin &A) {
  os << "a1 = " << A.a1 << " "
     << "a2 = " << A.a2 << " "
     << "a3 = " << A.a3 << " ";
  return os;
}
