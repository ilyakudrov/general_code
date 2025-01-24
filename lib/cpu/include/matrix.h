#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <iostream>

// su2 matrix in sigma matrices representation (a0 + ai * sigma[i])
class su2 {
public:
  inline static int data_size = 4;
  double a0, a1, a2, a3;
  su2() {
    a0 = 1.;
    a1 = 0;
    a2 = 0;
    a3 = 0;
  }

  su2(double b0, double b1, double b2, double b3) {
    a0 = b0;
    a1 = b1;
    a2 = b2;
    a3 = b3;
  }

  su2(double a) {
    a0 = a;
    a1 = a;
    a2 = a;
    a3 = a;
  }

  su2(const std::vector<double> &data) {
    a0 = data[0];
    a1 = data[1];
    a2 = data[2];
    a3 = data[3];
  }

  std::vector<double> get_data() {
    std::vector<double> data(4);
    data[0] = a0;
    data[1] = a1;
    data[2] = a2;
    data[3] = a3;
    return data;
  }

  double tr() { return a0; }

  double multiply_conj_tr(const su2 &B) const {
    return a0 * B.a0 + a1 * B.a1 + a2 * B.a2 + a3 * B.a3;
  }

  double multiply_tr(const su2 &B) const {
    return a0 * B.a0 - a1 * B.a1 - a2 * B.a2 - a3 * B.a3;
  }

  double multiply_conj_tr_adjoint(const su2 &B) const {
    double trace = a0 * B.a0 + a1 * B.a1 + a2 * B.a2 + a3 * B.a3;
    return 4 * trace * trace - 1;
  }

  su2 inverse() const {
    double rho = a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3;
    return su2(a0 / rho, -a1 / rho, -a2 / rho, -a3 / rho);
  }
  double module() const { return a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3; }
  su2 conj() const { return su2(a0, -a1, -a2, -a3); }
  su2 proj() const {
    double rho = sqrt(a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3);
    return su2(a0 / rho, a1 / rho, a2 / rho, a3 / rho);
  }
  su2 sigma3_mult() const { return su2(a0, -a1, -a2, a3); }

  friend su2 operator+(const su2 &A, const su2 &B) {
    return su2(A.a0 + B.a0, A.a1 + B.a1, A.a2 + B.a2, A.a3 + B.a3);
  }
  friend su2 operator-(const su2 &A, const su2 &B) {
    return su2(A.a0 - B.a0, A.a1 - B.a1, A.a2 - B.a2, A.a3 - B.a3);
  }
  friend su2 operator*(const double &x, const su2 &A) {
    return su2(A.a0 * x, A.a1 * x, A.a2 * x, A.a3 * x);
  }
  friend su2 operator*(const su2 &A, const double &x) {
    return su2(A.a0 * x, A.a1 * x, A.a2 * x, A.a3 * x);
  }

  friend su2 operator*(const su2 &A, const su2 &B) {
    return su2(A.a0 * B.a0 - A.a1 * B.a1 - A.a2 * B.a2 - A.a3 * B.a3,
               A.a0 * B.a1 + B.a0 * A.a1 + A.a3 * B.a2 - A.a2 * B.a3,
               A.a0 * B.a2 + B.a0 * A.a2 + A.a1 * B.a3 - A.a3 * B.a1,
               A.a0 * B.a3 + B.a0 * A.a3 + A.a2 * B.a1 - A.a1 * B.a2);
  }
  friend su2 operator^(const su2 &A, const su2 &B) {
    return su2(A.a0 * B.a0 + A.a1 * B.a1 + A.a2 * B.a2 + A.a3 * B.a3,
               -A.a0 * B.a1 + B.a0 * A.a1 - A.a3 * B.a2 + A.a2 * B.a3,
               -A.a0 * B.a2 + B.a0 * A.a2 - A.a1 * B.a3 + A.a3 * B.a1,
               -A.a0 * B.a3 + B.a0 * A.a3 - A.a2 * B.a1 + A.a1 * B.a2);
  }
  friend su2 operator%(const su2 &A, const su2 &B) {
    return su2(A.a0 * B.a0 + A.a1 * B.a1 + A.a2 * B.a2 + A.a3 * B.a3,
               A.a0 * B.a1 - B.a0 * A.a1 - A.a3 * B.a2 + A.a2 * B.a3,
               A.a0 * B.a2 - B.a0 * A.a2 - A.a1 * B.a3 + A.a3 * B.a1,
               A.a0 * B.a3 - B.a0 * A.a3 - A.a2 * B.a1 + A.a1 * B.a2);
  }

  friend std::ostream &operator<<(std::ostream &os, const su2 &A) {
    os << "a0 = " << A.a0 << " " << "a1 = " << A.a1 << " " << "a2 = " << A.a2
       << " " << "a3 = " << A.a3;
    return os;
  }
};

class abelian {
public:
  inline static int data_size = 1;
  double r, phi;
  abelian() {
    r = 1;
    phi = 0;
  }

  abelian(double r1, double phi1) {
    r = r1;
    phi = phi1;
  }

  abelian(double a) {
    r = a;
    phi = 0;
  }

  abelian(const std::vector<double> &data) {
    r = 1;
    phi = data[0];
  }

  std::vector<double> get_data() {
    std::vector<double> data(1);
    data[0] = phi;
    return data;
  }

  double tr() const { return r * cos(phi); }

  double multiply_conj_tr(const abelian &B) const {
    return r * B.r * cos(phi - B.phi);
  }

  double multiply_tr(const abelian &B) const {
    return r * B.r * cos(phi + B.phi);
  }

  double multiply_conj_tr_adjoint(const abelian &B) const {
    double trace = r * B.r * cos(2 * (phi - B.phi));
    return trace;
  }

  abelian inverse() const { return abelian(1 / r, -phi); }
  double module() const { return r; }
  abelian conj() const { return abelian(r, -phi); }
  abelian proj() const { return abelian(1, phi); }

  friend abelian operator+(const abelian &A, const abelian &B) {
    return abelian(sqrt((A.r * sin(A.phi) + B.r * sin(B.phi)) *
                            (A.r * sin(A.phi) + B.r * sin(B.phi)) +
                        (A.r * cos(A.phi) + B.r * cos(B.phi)) *
                            (A.r * cos(A.phi) + B.r * cos(B.phi))),
                   atan2(A.r * sin(A.phi) + B.r * sin(B.phi),
                         A.r * cos(A.phi) + B.r * cos(B.phi)));
  }
  friend abelian operator-(const abelian &A, const abelian &B) {
    return abelian(sqrt((A.r * sin(A.phi) - B.r * sin(B.phi)) *
                            (A.r * sin(A.phi) - B.r * sin(B.phi)) +
                        (A.r * cos(A.phi) - B.r * cos(B.phi)) *
                            (A.r * cos(A.phi) - B.r * cos(B.phi))),
                   atan2(A.r * sin(A.phi) - B.r * sin(B.phi),
                         A.r * cos(A.phi) - B.r * cos(B.phi)));
  }
  friend abelian operator*(const double &x, const abelian &A) {
    return abelian(A.r * x, A.phi);
  }
  friend abelian operator*(const abelian &A, const double &x) {
    return abelian(A.r * x, A.phi);
  }

  friend abelian operator*(const abelian &A, const abelian &B) {
    return abelian(A.r * B.r, A.phi + B.phi);
  }
  friend abelian operator^(const abelian &A, const abelian &B) {
    return abelian(A.r * B.r, A.phi - B.phi);
  }
  friend abelian operator%(const abelian &A, const abelian &B) {
    return abelian(A.r * B.r, B.phi - A.phi);
  }

  friend std::ostream &operator<<(std::ostream &os, const abelian &A) {
    os << "r = " << A.r << " " << "phi = " << A.phi << " " << std::endl;
    return os;
  }
};

class su3 {
public:
  inline static int data_size = 18;
  Eigen::Matrix3cd matrix;
  su3() {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (i != j)
          matrix(i, j) = std::complex<double>(0, 0);
        else
          matrix(i, j) = std::complex<double>(1., 0);
      }
    }
  }

  su3(const Eigen::Matrix3cd &B) { matrix = B; }

  su3(double a) {
    matrix << std::complex<double>(a, 0), std::complex<double>(a, 0),
        std::complex<double>(a, 0), std::complex<double>(a, 0),
        std::complex<double>(a, 0), std::complex<double>(a, 0),
        std::complex<double>(a, 0), std::complex<double>(a, 0),
        std::complex<double>(a, 0);
  }

  su3(const std::vector<double> &data) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        matrix(i, j) = std::complex<double>(data[(j + i * 3) * 2],
                                            data[(j + i * 3) * 2 + 1]);
      }
    }
  }

  std::vector<double> get_data() {
    std::vector<double> data(18);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        data[(j + i * 3) * 2] = matrix(i, j).real();
        data[(j + i * 3) * 2 + 1] = matrix(i, j).imag();
      }
    }
    return data;
  }

  double tr() const { return matrix.trace().real() / 3.; }

  std::complex<double> tr_complex() const { return matrix.trace() / 3.; }

  double multiply_conj_tr(const su3 &B) const {
    return (matrix * B.matrix.adjoint()).trace().real() / 3;
  }

  double multiply_tr(const su3 &B) const {
    return (matrix * B.matrix).trace().real() / 3;
  }

  double multiply_conj_tr_adjoint(const su3 &B) const {
    return std::norm((matrix * B.matrix.adjoint()).trace()) - 1;
  }

  su3 conj() const { return su3(matrix.adjoint()); }

  su3 mult_by_imag(const double x) const {
    return su3(matrix * std::complex<double>(0, x));
  }

  su3 proj() {
    Eigen::Matrix3cd A;
    A = (3. / 2) * matrix - 0.5 * matrix * matrix.adjoint() * matrix;
    A = (std::complex<double>(1, 0) -
         (1. / 3) * (A.determinant() - std::complex<double>(1, 0))) *
        A;
    for (int i = 0; i < 5; i++) {
      A = (3. / 2) * A - 0.5 * A * A.adjoint() * A;
      A = (std::complex<double>(1, 0) -
           (1. / 3) * (A.determinant() - std::complex<double>(1, 0))) *
          A;
    }
    return su3(A);
  }

  su3 lambda3_mult() {

    Eigen::Matrix3cd A = Eigen::Matrix3cd::Identity(3, 3);

    A(0, 0) = matrix(0, 0);
    A(0, 1) = -matrix(0, 1);
    A(1, 0) = matrix(1, 0);
    A(1, 1) = -matrix(1, 1);
    A(2, 0) = A(2, 0);
    A(2, 1) = -A(2, 1);

    A(2, 2) = std::complex<double>(0, 0);

    return su3(A);
  }

  su3 lambda8_mult() {

    Eigen::Matrix3cd A = matrix;

    A(0, 0) = A(0, 0) / 3.;
    A(0, 1) = A(0, 1) / 3.;
    A(1, 0) = A(1, 0) / 3.;
    A(1, 1) = A(1, 1) / 3.;

    A(0, 2) = A(0, 2) * (-2. / 3);
    A(1, 2) = A(1, 2) * (-2. / 3);
    A(2, 0) = A(2, 0) * (-2. / 3);
    A(2, 1) = A(2, 1) * (-2. / 3);

    A(2, 2) = A(2, 2) * (4. / 3);

    return su3(A);
  }

  friend su3 operator+(const su3 &A, const su3 &B) {
    return su3(A.matrix + B.matrix);
  }
  friend su3 operator-(const su3 &A, const su3 &B) {
    return su3(A.matrix - B.matrix);
  }
  friend su3 operator*(const double &x, const su3 &A) {
    return su3(x * A.matrix);
  }
  friend su3 operator*(const su3 &A, const double &x) {
    return su3(A.matrix * x);
  }
  friend su3 operator*(const std::complex<double> &x, const su3 &A) {
    return su3(x * A.matrix);
  }
  friend su3 operator*(const su3 &A, const std::complex<double> &x) {
    return su3(A.matrix * x);
  }

  friend su3 operator*(const su3 &A, const su3 &B) {
    return su3(A.matrix * B.matrix);
  }

  friend su3 operator^(const su3 &A, const su3 &B) {
    return su3(A.matrix * B.matrix.adjoint().eval());
  }

  friend su3 operator%(const su3 &A, const su3 &B) {
    return su3(A.matrix.adjoint().eval() * B.matrix);
  }

  friend std::ostream &operator<<(std::ostream &os, const su3 &A) {
    os << A.matrix << std::endl;
    return os;
  }
};

class su3_abelian {
public:
  inline static int data_size = 3;
  std::complex<double> matrix[3];
  su3_abelian(std::complex<double> B[3]) {
    for (int i = 0; i < 3; i++) {
      matrix[i] = B[i];
    }
  }
  su3_abelian(std::vector<std::complex<double>> &B) {
    for (int i = 0; i < 3; i++) {
      matrix[i] = B[i];
    }
  }

  su3_abelian() {
    for (int i = 0; i < 3; i++) {
      matrix[i] = std::complex<double>(1, 0);
    }
  }

  su3_abelian(double a) {
    matrix[0] = std::complex<double>(a, 0);
    matrix[1] = std::complex<double>(a, 0);
    matrix[2] = std::complex<double>(a, 0);
  }

  su3_abelian(const std::vector<double> &data) {
    for (int i = 0; i < 3; i++) {
      matrix[i] = std::complex<double>(cos(data[i]), sin(data[i]));
    }
  }

  std::vector<double> get_data() {
    std::vector<double> data(3);
    for (int i = 0; i < 3; i++) {
      data[i] = atan2(matrix[i].imag(), matrix[i].real());
    }
    return data;
  }

  double tr() const {
    return (matrix[0].real() + matrix[1].real() + matrix[2].real()) / 3;
  }

  std::complex<double> tr_complex() const {
    return (matrix[0] + matrix[1] + matrix[2]) / 3.;
  }

  double multiply_conj_tr(const su3_abelian &B) const {
    double trace = 0;
    for (int i = 0; i < 3; i++) {
      trace += matrix[i].real() * B.matrix[i].real() +
               matrix[i].imag() * B.matrix[i].imag();
    }
    return trace / 3;
  }

  double multiply_tr(const su3_abelian &B) const {
    double trace = 0;
    for (int i = 0; i < 3; i++) {
      trace += matrix[i].real() * B.matrix[i].real() -
               matrix[i].imag() * B.matrix[i].imag();
    }
    return trace / 3;
  }

  double multiply_conj_tr_adjoint(const su3_abelian &B) const {
    std::complex<double> tmp = std::complex<double>(0, 0);
    for (int i = 0; i < 3; i++) {
      tmp += matrix[i] * std::conj(B.matrix[i]);
    }
    return std::norm(tmp) - 1;
  }

  su3_abelian inverse() const {
    su3_abelian B;
    for (int i = 0; i < 3; i++) {
      B.matrix[i] = std::conj(matrix[i]);
    }
    return B;
  }

  double module() const {
    std::complex<double> tmp = matrix[0] * matrix[1];
    return tmp.real() * matrix[2].real() - tmp.imag() * matrix[2].imag();
  }

  std::complex<double> determinant() const {
    return matrix[0] * matrix[1] * matrix[2];
  }

  su3_abelian conj() const {
    su3_abelian B;
    for (int i = 0; i < 3; i++) {
      B.matrix[i] = std::conj(matrix[i]);
    }
    return B;
  }

  su3_abelian mult_by_imag(double x) const {
    su3_abelian C;
    for (int i = 0; i < 3; i++) {
      C.matrix[i] = matrix[i] * std::complex<double>(0, x);
    }
    return C;
  }

  su3_abelian proj() const {
    su3_abelian A;
    std::vector<double> angles(3);
    double sum = 0;
    for (int i = 0; i < 3; i++) {
      angles[i] = atan2(matrix[i].imag(), matrix[i].real());
      sum += angles[i];
    }
    while (sum >= M_PI) {
      sum -= 2 * M_PI;
    }
    while (sum < -M_PI) {
      sum += 2 * M_PI;
    }
    for (int i = 0; i < 3; i++) {
      A.matrix[i] = std::complex<double>(cos(angles[i] - sum / 3),
                                         sin(angles[i] - sum / 3));
    }
    return A;
  }

  friend su3_abelian operator+(const su3_abelian &A, const su3_abelian &B) {
    su3_abelian C;
    for (int i = 0; i < 3; i++) {
      C.matrix[i] = A.matrix[i] + B.matrix[i];
    }
    return C;
  }
  friend su3_abelian operator-(const su3_abelian &A, const su3_abelian &B) {
    su3_abelian C;
    for (int i = 0; i < 3; i++) {
      C.matrix[i] = A.matrix[i] - B.matrix[i];
    }
    return C;
  }
  friend su3_abelian operator*(const double &x, const su3_abelian &A) {
    su3_abelian C;
    for (int i = 0; i < 3; i++) {
      C.matrix[i] = x * A.matrix[i];
    }
    return C;
  }
  friend su3_abelian operator*(const su3_abelian &A, const double &x) {
    su3_abelian C;
    for (int i = 0; i < 3; i++) {
      C.matrix[i] = x * A.matrix[i];
    }
    return C;
  }
  friend su3_abelian operator*(const std::complex<double> &x,
                               const su3_abelian &A) {
    su3_abelian C;
    for (int i = 0; i < 3; i++) {
      C.matrix[i] = x * A.matrix[i];
    }
    return C;
  }
  friend su3_abelian operator*(const su3_abelian &A,
                               const std::complex<double> &x) {
    su3_abelian C;
    for (int i = 0; i < 3; i++) {
      C.matrix[i] = x * A.matrix[i];
    }
    return C;
  }

  friend su3_abelian operator*(const su3_abelian &A, const su3_abelian &B) {
    su3_abelian C;
    C.matrix[0] = A.matrix[0] * B.matrix[0];
    C.matrix[1] = A.matrix[1] * B.matrix[1];
    C.matrix[2] = A.matrix[2] * B.matrix[2];
    return C;
  }

  friend su3_abelian operator^(const su3_abelian &A, const su3_abelian &B) {
    su3_abelian C;
    C.matrix[0] = A.matrix[0] * std::conj(B.matrix[0]);
    C.matrix[1] = A.matrix[1] * std::conj(B.matrix[1]);
    C.matrix[2] = A.matrix[2] * std::conj(B.matrix[2]);
    return C;
  }

  friend su3_abelian operator%(const su3_abelian &A, const su3_abelian &B) {
    su3_abelian C;
    C.matrix[0] = std::conj(A.matrix[0]) * B.matrix[0];
    C.matrix[1] = std::conj(A.matrix[1]) * B.matrix[1];
    C.matrix[2] = std::conj(A.matrix[2]) * B.matrix[2];
    return C;
  }

  friend std::ostream &operator<<(std::ostream &os, const su3_abelian &A) {
    for (int i = 0; i < 3; i++) {
      os << "(" << A.matrix[i].real() << ", " << A.matrix[i].imag() << ") ";
    }
    os << std::endl;
    return os;
  }
};

class su3_angles {
public:
  inline static int data_size = 3;
  double matrix[3];
  su3_angles(double B[3]) {
    for (int i = 0; i < 3; i++) {
      matrix[i] = B[i];
    }
  }

  su3_angles() {
    for (int i = 0; i < 3; i++) {
      matrix[i] = 0;
    }
  }
  su3_angles(std::vector<double> B) {
    for (int i = 0; i < 3; i++) {
      matrix[i] = B[i];
    }
  }
  su3_angles(double a) {
    matrix[0] = 0;
    matrix[1] = 0;
    matrix[2] = 0;
  }

  std::vector<double> get_data() {
    std::vector<double> data(3);
    for (int i = 0; i < 3; i++) {
      data[i] = matrix[i];
    }
    return data;
  }

  double tr() const {
    return (cos(matrix[0]) + cos(matrix[1]) + cos(matrix[2])) / 3;
  }

  std::complex<double> tr_complex() const {
    return std::complex<double>(
        (cos(matrix[0]) + cos(matrix[1]) + cos(matrix[2])) / 3,
        (sin(matrix[0]) + sin(matrix[1]) + sin(matrix[2])) / 3);
  }

  double multiply_conj_tr(const su3_angles &B) const {
    double trace = 0;
    for (int i = 0; i < 3; i++) {
      trace += cos(matrix[i] - B.matrix[i]);
    }
    return trace / 3;
  }

  double multiply_tr(const su3_angles &B) const {
    double trace = 0;
    for (int i = 0; i < 3; i++) {
      trace += cos(matrix[i] + B.matrix[i]);
    }
    return trace / 3;
  }

  double multiply_conj_tr_adjoint(const su3_angles &B) const {
    std::complex<double> tmp(0, 0);
    for (int i = 0; i < 3; i++) {
      tmp += std::complex<double>(cos(matrix[i] - B.matrix[i]),
                                  sin(matrix[i] - B.matrix[i]));
    }
    return std::norm(tmp) - 1;
  }

  su3_angles inverse() const {
    su3_angles B;
    for (int i = 0; i < 3; i++) {
      B.matrix[i] = -matrix[i];
    }
    return B;
  }

  double module() const { return 1; }

  std::complex<double> determinant() const {
    return std::complex<double>(1, 0);
  }

  su3_angles conj() const {
    su3_angles B;
    for (int i = 0; i < 3; i++) {
      B.matrix[i] = -matrix[i];
    }
    return B;
  }

  su3_angles mult_by_imag(double x) const {
    su3_angles C;
    for (int i = 0; i < 3; i++) {
      C.matrix[i] = matrix[i] + M_PI;
    }
    return C;
  }

  su3_angles proj() const {
    su3_angles A;
    double sum = 0;
    for (int i = 0; i < 3; i++) {
      sum += matrix[i];
    }
    while (sum >= M_PI) {
      sum -= 2 * M_PI;
    }
    while (sum < -M_PI) {
      sum += 2 * M_PI;
    }
    for (int i = 0; i < 3; i++) {
      A.matrix[i] = matrix[i] - sum / 3;
    }
    return A;
  }

  friend su3_angles operator+(const su3_angles &A, const su3_angles &B) {
    su3_angles C;
    for (int i = 0; i < 3; i++) {
      C.matrix[i] = atan2(sin(A.matrix[i]) + sin(B.matrix[i]),
                          cos(A.matrix[i]) + cos(B.matrix[i]));
    }
    return C;
  }
  friend su3_angles operator-(const su3_angles &A, const su3_angles &B) {
    su3_angles C;
    for (int i = 0; i < 3; i++) {
      C.matrix[i] = atan2(sin(A.matrix[i]) - sin(B.matrix[i]),
                          cos(A.matrix[i]) - cos(B.matrix[i]));
    }
    return C;
  }
  friend su3_angles operator*(const double &x, const su3_angles &A) {
    su3_angles C;
    for (int i = 0; i < 3; i++) {
      C.matrix[i] = A.matrix[i];
    }
    return C;
  }
  friend su3_angles operator*(const su3_angles &A, const double &x) {
    su3_angles C;
    for (int i = 0; i < 3; i++) {
      C.matrix[i] = A.matrix[i];
    }
    return C;
  }
  friend su3_angles operator*(const std::complex<double> &x,
                              const su3_angles &A) {
    su3_angles C;
    for (int i = 0; i < 3; i++) {
      C.matrix[i] = atan2(x.imag(), x.real()) + A.matrix[i];
    }
    return C;
  }
  friend su3_angles operator*(const su3_angles &A,
                              const std::complex<double> &x) {
    su3_angles C;
    for (int i = 0; i < 3; i++) {
      C.matrix[i] = atan2(x.imag(), x.real()) + A.matrix[i];
    }
    return C;
  }
  friend su3_angles operator*(const su3_angles &A, const su3_angles &B) {
    su3_angles C;
    for (int i = 0; i < 3; i++) {
      C.matrix[i] = A.matrix[i] + B.matrix[i];
    }
    return C;
  }
  friend su3_angles operator^(const su3_angles &A, const su3_angles &B) {
    su3_angles C;
    for (int i = 0; i < 3; i++) {
      C.matrix[i] = A.matrix[i] - B.matrix[i];
    }
    return C;
  }
  friend su3_angles operator%(const su3_angles &A, const su3_angles &B) {
    su3_angles C;
    for (int i = 0; i < 3; i++) {
      C.matrix[i] = B.matrix[i] - A.matrix[i];
    }
    return C;
  }
  friend std::ostream &operator<<(std::ostream &os, const su3_angles &A) {
    for (int i = 0; i < 3; i++) {
      os << A.matrix[i] << " ";
    }
    os << std::endl;
    return os;
  }
};

// 3D vector realisation for spin model.
class spin {
public:
  double a1, a2, a3;

  spin() {
    a1 = (double)0.0;
    a2 = (double)0.0;
    a3 = (double)1.0;
  }
  spin(su2 U) {
    a1 = 2 * (U.a0 * U.a2 + U.a3 * U.a1);
    a2 = 2 * (U.a3 * U.a2 + U.a0 * U.a1);
    a3 = 1. - 2 * (U.a2 * U.a2 + U.a1 * U.a1);
  }
  spin(double a1, double a2, double a3) : a1(a1), a2(a2), a3(a3) {}

  double norm() const { return sqrt(a1 * a1 + a2 * a2 + a3 * a3); }

  void normalize() {
    double norm = sqrt(a1 * a1 + a2 * a2 + a3 * a3);
    a1 = a1 / norm;
    a2 = a2 / norm;
    a3 = a3 / norm;
  }

  // Function reflect spin vector through arbitrary spin vector
  // it's more convenient to just change the spin variable
  double reflect(spin &V) {
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

  void reflect_fast(spin &V) {
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

  double parallel(spin &V) {
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

  void parallel_fast(spin &V) {
    double vNorm = V.norm();

    a1 = V.a1 / vNorm;
    a2 = V.a2 / vNorm;
    a3 = V.a3 / vNorm;
  }

  bool IsUnit() { return (this->norm() - 1.) > 1e-6 ? false : true; }

  // vector which spin variable is multiplied on
  // it's contribution from single neighbour site in one direction
  spin contribution(const su2 &A) const {
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
  spin contribution_conj(const su2 &A) const {
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

  void contribution1(const su2 &A, const spin &b) {
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

    a1 +=
        b.a1 * (A5 + 2 * q1) + 2 * b.a2 * (q12 + q03) + 2 * b.a3 * (q13 - q02);
    a2 +=
        2 * b.a1 * (q12 - q03) + b.a2 * (A5 + 2 * q2) + 2 * b.a3 * (q23 + q01);
    a3 +=
        2 * b.a1 * (q13 + q02) + 2 * b.a2 * (q23 - q01) + b.a3 * (A5 + 2 * q3);
  }

  void contribution1_conj(const su2 &A, const spin &b) {
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

    a1 +=
        b.a1 * (A5 + 2 * q1) + 2 * b.a2 * (q12 - q03) + 2 * b.a3 * (q13 + q02);
    a2 +=
        2 * b.a1 * (q12 + q03) + b.a2 * (A5 + 2 * q2) + 2 * b.a3 * (q23 - q01);
    a3 +=
        2 * b.a1 * (q13 - q02) + 2 * b.a2 * (q23 + q01) + b.a3 * (A5 + 2 * q3);
  }

  // Calculation have been done for specific element of gauge matrix G with g_3
  // = 0 Here:
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
  su2 GetGaugeMatrix() {
    double sq = sqrt(2. - 2 * a3);
    return su2(-a1 / sq, 0., -sq / 2., -a2 / sq);
  }
  friend spin operator*(const double &x, const spin &A) {
    return spin(A.a1 * x, A.a2 * x, A.a3 * x);
  }
  friend spin operator*(const spin &A, const double &x) {
    return spin(A.a1 * x, A.a2 * x, A.a3 * x);
  }
  friend spin operator+(const spin &A, const spin &B) {
    return spin(A.a1 + B.a1, A.a2 + B.a2, A.a3 + B.a3);
  }
  friend spin operator-(const spin &A, const spin &B) {
    return spin(A.a1 - B.a1, A.a2 - B.a2, A.a3 - B.a3);
  }
  friend double operator*(const spin &A, const spin &B) {
    return A.a1 * B.a1 + A.a2 * B.a2 + A.a3 * B.a3;
  }

  friend std::ostream &operator<<(std::ostream &os, const spin &A) {
    os << "a1 = " << A.a1 << " " << "a2 = " << A.a2 << " " << "a3 = " << A.a3
       << " ";
    return os;
  }
};

inline std::vector<su3> get_generators_su3() {
  std::vector<su3> generators(8);

  Eigen::Matrix3cd matrix = Eigen::Matrix3cd::Identity(3, 3);

  // lambda1
  matrix(0, 1) = std::complex<double>(1, 0);
  matrix(1, 0) = std::complex<double>(1, 0);
  generators[0] = su3(matrix);

  // lambda2
  matrix(0, 1) = std::complex<double>(0, -1);
  matrix(1, 0) = std::complex<double>(0, 1);
  generators[1] = su3(matrix);

  // lambda3
  matrix(0, 1) = std::complex<double>(0, 0);
  matrix(1, 0) = std::complex<double>(0, 0);
  matrix(0, 0) = std::complex<double>(1, 0);
  matrix(1, 1) = std::complex<double>(-1, 0);
  generators[2] = su3(matrix);

  // lambda4
  matrix(0, 0) = std::complex<double>(0, 0);
  matrix(1, 1) = std::complex<double>(0, 0);
  matrix(0, 2) = std::complex<double>(1, 0);
  matrix(2, 0) = std::complex<double>(1, 0);
  generators[3] = su3(matrix);

  // lambda5
  matrix(0, 2) = std::complex<double>(0, -1);
  matrix(2, 0) = std::complex<double>(0, 1);
  generators[4] = su3(matrix);

  // lambda6
  matrix(0, 2) = std::complex<double>(0, 0);
  matrix(2, 0) = std::complex<double>(0, 0);
  matrix(1, 2) = std::complex<double>(1, 0);
  matrix(2, 1) = std::complex<double>(1, 0);
  generators[5] = su3(matrix);

  // lambda7
  matrix(1, 2) = std::complex<double>(0, -1);
  matrix(2, 1) = std::complex<double>(0, 1);
  generators[6] = su3(matrix);

  // lambda8
  matrix(1, 2) = std::complex<double>(0, 0);
  matrix(2, 1) = std::complex<double>(0, 0);
  matrix(0, 0) = std::complex<double>(1 / sqrt(3), 0);
  matrix(1, 1) = std::complex<double>(1 / sqrt(3), 0);
  matrix(2, 2) = std::complex<double>(-2 / sqrt(3), 0);
  generators[7] = su3(matrix);

  return generators;
}

// class to convert MatrixType1 into MatrixType2
template <class MatrixType1, class MatrixType2> class MatrixConverter {
public:
  MatrixType2 convert_matrix(const MatrixType1 &A) {
    std::cout << "conversion from " << typeid(MatrixType1).name() << " to "
              << typeid(MatrixType2).name() << " is not implemented"
              << std::endl;
    return MatrixType2();
  }
};

template <>
inline abelian MatrixConverter<su2, abelian>::convert_matrix(const su2 &A) {
  return abelian(1, atan2(A.a3, A.a0));
}

template <>
inline su3_abelian
MatrixConverter<su3, su3_abelian>::convert_matrix(const su3 &A) {
  std::vector<double> angles(3);
  su3_abelian B;
  double sum = 0;
  for (int c = 0; c < 3; c++) {
    angles[c] = atan2(A.matrix(c, c).imag(), A.matrix(c, c).real());
    sum += angles[c];
  }
  while (sum >= M_PI) {
    sum -= 2 * M_PI;
  }
  while (sum < -M_PI) {
    sum += 2 * M_PI;
  }
  for (int c = 0; c < 3; c++) {
    B.matrix[c] = std::complex<double>(cos(angles[c] - sum / 3),
                                       sin(angles[c] - sum / 3));
  }
  return B;
}

template <> inline su2 MatrixConverter<su2, su2>::convert_matrix(const su2 &A) {
  double module = sqrt(A.a0 * A.a0 + A.a3 * A.a3);
  su2 B = su2(A.a0 / module, 0, 0, A.a3 / module);
  return A ^ B;
}

template <> inline su3 MatrixConverter<su3, su3>::convert_matrix(const su3 &A) {
  su3 B;
  for (int c = 0; c < 3; c++) {
    B.matrix(c, c) = A.matrix(c, c) / std::sqrt(std::norm(A.matrix(c, c)));
  }
  return A ^ B;
}