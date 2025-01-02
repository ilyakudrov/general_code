#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <iostream>

class complex_test {
public:
  double real;
  double imag;

  complex_test(const double real1, const double imag1) {
    real = real1;
    imag = imag1;
  }

  complex_test() {
    real = 0;
    imag = 0;
  }

  double module() { return sqrt(real * real + imag * imag); }

  double norm2() { return real * real + imag * imag; }

  double angle() { return atan2(imag, real); }

  complex_test conj() { return complex_test(real, -imag); }

  complex_test negative() { return complex_test(-real, -imag); }

  complex_test mult_by_imag(double x) {
    return complex_test(-imag * x, real * x);
  }

  void add(double r, double i) {
    real += r;
    imag += i;
  }

  complex_test sqrt_complex() {
    double module = sqrt(sqrt(real * real + imag * imag));

    double phi = atan2(imag, real) / 2;

    return complex_test(module * cos(phi), module * sin(phi));
  }

  complex_test &operator+=(const complex_test &a) {
    real += a.real;
    imag += a.imag;

    return *this;
  }

  complex_test &operator-=(const complex_test &a) {
    real -= a.real;
    imag -= a.imag;

    return *this;
  }

  complex_test &operator*=(const complex_test &a) {
    double tmp1 = real * a.real - imag * a.imag;
    imag = real * a.imag + imag * a.real;
    real = tmp1;

    return *this;
  }

  complex_test &operator*=(const double &a) {
    real *= a;
    imag *= a;

    return *this;
  }

  complex_test &operator/=(const complex_test &a) {
    double norm_factor = 1. / (a.real * a.real + a.imag * a.imag);
    double tmp = (real * a.real + imag * a.imag) * norm_factor;

    imag = (imag * a.real - real * a.imag) * norm_factor;
    real = tmp;

    return *this;
  }

  complex_test &operator/=(const double &a) {
    real /= a;
    imag /= a;

    return *this;
  }

  friend complex_test operator+(const complex_test &a, const complex_test &b) {
    return complex_test(a.real + b.real, a.imag + b.imag);
  }

  friend complex_test operator-(const complex_test &a, const complex_test &b) {
    return complex_test(a.real - b.real, a.imag - b.imag);
  }

  friend complex_test operator*(const complex_test &a, const complex_test &b) {
    return complex_test(a.real * b.real - a.imag * b.imag,
                        a.real * b.imag + a.imag * b.real);
  }

  friend complex_test operator*(const double &a, const complex_test &b) {
    return complex_test(a * b.real, a * b.imag);
  }

  friend complex_test operator*(const complex_test &a, const double &b) {
    return complex_test(a.real * b, a.imag * b);
  }

  friend complex_test operator^(const complex_test &a, const complex_test &b) {
    return complex_test(a.real * b.real + a.imag * b.imag,
                        a.imag * b.real - a.real * b.imag);
  }

  friend complex_test operator%(const complex_test &a, const complex_test &b) {
    return complex_test(a.real * b.real + a.imag * b.imag,
                        a.real * b.imag - a.imag * b.real);
  }

  friend complex_test operator&(const complex_test &a, const complex_test &b) {
    return complex_test(a.real * b.real - a.imag * b.imag,
                        -a.imag * b.real - a.real * b.imag);
  }

  friend complex_test operator/(const complex_test &a, const double &b) {
    return complex_test(a.real / b, a.imag / b);
  }

  friend complex_test operator/(const complex_test &a, const complex_test &b) {
    double norm_factor = 1. / (b.real * b.real + b.imag * b.imag);
    return complex_test((a.real * b.real + a.imag * b.imag) * norm_factor,
                        (a.imag * b.real - a.real * b.imag) * norm_factor);
  }

  friend std::ostream &operator<<(std::ostream &os, const complex_test &a) {
    os << "(" << a.real << ", " << a.imag << ")";
    return os;
  }
};

class su3_test {
public:
  complex_test matrix[3][3];
  su3_test() {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (i != j)
          matrix[i][j] = complex_test(0, 0);
        else
          matrix[i][j] = complex_test(1., 0);
      }
    }
  }

  su3_test(complex_test B[3][3]) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        matrix[i][j] = B[i][j];
      }
    }
  }

  double tr() {
    return (matrix[0][0].real + matrix[1][1].real + matrix[2][2].real) / 3;
  }

  complex_test tr_complex() {
    return (matrix[0][0] + matrix[1][1] + matrix[2][2]) / 3;
  }

  double multiply_conj_tr(const su3_test &B) {
    double trace = 0;
    for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
        trace += matrix[i][k].real * B.matrix[i][k].real +
                 matrix[i][k].imag * B.matrix[i][k].imag;
      }
    }
    return trace / 3;
  }

  double multiply_tr(const su3_test &B) {
    double trace = 0;
    for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
        trace += matrix[i][k].real * B.matrix[k][i].real -
                 matrix[i][k].imag * B.matrix[k][i].imag;
      }
    }
    return trace / 3;
  }

  double multiply_conj_tr_adjoint(const su3_test &B) {
    complex_test trace = complex_test(0, 0);
    for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
        trace += matrix[i][k] ^ B.matrix[i][k];
      }
    }

    return trace.norm2() - 1;
  }

  su3_test inverse() {
    su3_test B;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        B.matrix[i][j].real = matrix[j][i].real;
        B.matrix[i][j].imag = -matrix[j][i].imag;
      }
    }
    return B;
  }

  complex_test determinant_part(complex_test B[3][3], int i, int j, int k) {
    return B[0][i] * (B[1][j] * B[2][k] - B[2][j] * B[1][k]);
  }

  double module() {
    complex_test determinant;
    determinant = determinant_part(matrix, 0, 1, 2);
    determinant = determinant + determinant_part(matrix, 1, 0, 2);
    determinant = determinant + determinant_part(matrix, 2, 0, 1);
    return determinant.real;
  }

  complex_test determinant() {
    complex_test determinant;
    determinant = determinant_part(matrix, 0, 1, 2);
    determinant = determinant - determinant_part(matrix, 1, 0, 2);
    determinant = determinant + determinant_part(matrix, 2, 0, 1);
    return determinant;
  }

  complex_test unitarity_check() {
    su3_test A = this->conj();
    A = A * *this;

    complex_test diagonal = complex_test(0, 0),
                 non_diagonal = complex_test(0, 0);
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

  su3_test conj() const {
    su3_test B;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        B.matrix[i][j].real = matrix[j][i].real;
        B.matrix[i][j].imag = -matrix[j][i].imag;
      }
    }
    return B;
  }

  su3_test mult_by_imag(double x) {
    su3_test C;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        C.matrix[i][j] = matrix[i][j].mult_by_imag(x);
      }
    }
    return C;
  }

  void project_su2(complex_test v[][2]) {
    complex_test u[2][2];

    // v*v^dagger
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        u[i][j] = complex_test(0, 0);
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

    complex_test w[2][2];

    w[0][0] = complex_test((y1 - u[1][1].real), 0) / u[1][0];
    w[0][1] = complex_test((y2 - u[1][1].real), 0) / u[1][0];

    double w_norm = 1 / sqrt(w[0][0].norm2() + 1);
    w[0][0] *= w_norm;
    w[1][0] = complex_test(w_norm, 0);
    w_norm = 1 / sqrt(w[0][1].norm2() + 1);
    w[0][1] *= w_norm;
    w[1][1] = complex_test(w_norm, 0);

    // v^dagger * v
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        u[i][j] = complex_test(0, 0);
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

    complex_test r[2][2];

    r[0][0] = complex_test((y1 - u[1][1].real), 0) / u[1][0];
    r[0][1] = complex_test((y2 - u[1][1].real), 0) / u[1][0];

    double r_norm = 1 / sqrt(r[0][0].norm2() + 1);
    r[0][0] *= r_norm;
    r[1][0] = complex_test(r_norm, 0);
    r_norm = 1 / sqrt(r[0][1].norm2() + 1);
    r[0][1] *= r_norm;
    r[1][1] = complex_test(r_norm, 0);

    complex_test dw = w[0][0] * w[1][1] - w[1][0] * w[0][1];
    complex_test dr = r[0][0] * r[1][1] - r[1][0] * r[0][1];

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        w[i][j] /= dw.sqrt_complex();
        r[i][j] /= dr.sqrt_complex();
      }
    }

    dw = w[0][0] * w[1][1] - w[1][0] * w[0][1];
    dr = r[0][0] * r[1][1] - r[1][0] * r[0][1];

    complex_test ax1[2][2];

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        ax1[i][j] = complex_test(0, 0);
        for (int k = 0; k < 2; k++) {
          ax1[i][j] += (w[k][i] % v[k][j]);
        }
      }
    }

    complex_test ax2[2][2];

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        ax2[i][j] = complex_test(0, 0);
        for (int k = 0; k < 2; k++) {
          ax2[i][j] += ax1[i][k] * r[k][j];
        }
      }
    }

    complex_test tmp = ax2[0][0] + ax2[1][1].conj();

    ax2[0][0] = tmp / tmp.module();
    ax2[1][1] = ax2[0][0].conj();
    ax2[0][1] = complex_test(0, 0);
    ax2[1][0] = complex_test(0, 0);

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        ax1[i][j] = complex_test(0, 0);
        for (int k = 0; k < 2; k++) {
          ax1[i][j] += w[i][k] * ax2[k][j];
        }
      }
    }

    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        v[i][j] = complex_test(0, 0);
        for (int k = 0; k < 2; k++) {
          v[i][j] += ax1[i][k] ^ r[j][k];
        }
      }
    }
  }

  su3_test proj() {
    su3_test A;
    A = (3. / 2) * *this - 0.5 * *this * this->conj() * *this;
    A = (complex_test(1, 0) -
         (1. / 3) * (A.determinant() - complex_test(1, 0))) *
        A;
    for (int i = 0; i < 5; i++) {
      A = (3. / 2) * A - 0.5 * A * A.conj() * A;
      A = (complex_test(1, 0) -
           (1. / 3) * (A.determinant() - complex_test(1, 0))) *
          A;
    }
    return A;
  }

  su3_test proj1() {

    su3_test A;

    complex_test tmp1(0, 0);
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
    tmp1 = complex_test(0, 0);
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
    tmp1 = complex_test(0, 0);
    complex_test tmp2(0, 0);
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
    complex_test determinant = A.determinant();

    double ph = atan2(determinant.imag, determinant.real) / 3;

    tmp1 = complex_test(cos(ph), -sin(ph));

    A = A * tmp1;

    su3_test v = A % *this;

    complex_test v1[2][2];

    for (int iter = 0; iter < 3; iter++) {

      // first subgroup

      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
          v1[i][j] = v.matrix[i][j];
        }
      }

      project_su2(v1);

      su3_test s3;

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

      s3 = su3_test();

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

      s3 = su3_test();

      s3.matrix[1][1] = v1[0][0];
      s3.matrix[2][2] = v1[1][1];
      s3.matrix[2][1] = v1[1][0];
      s3.matrix[1][2] = v1[0][1];

      v = s3 % v;
      A = A * s3;
    }

    return A;
  }

  su3_test lambda3_mult() {

    su3_test A;

    A.matrix[0][0] = matrix[0][0];
    A.matrix[0][1] = matrix[0][1].negative();
    A.matrix[1][0] = matrix[1][0];
    A.matrix[1][1] = matrix[1][1].negative();
    A.matrix[2][0] = A.matrix[2][0];
    A.matrix[2][1] = A.matrix[2][1].negative();

    A.matrix[2][2] = complex_test(0, 0);

    return su3_test(A);
  }

  su3_test lambda8_mult() {

    su3_test A(matrix);

    A.matrix[0][0] = A.matrix[0][0] / 3;
    A.matrix[0][1] = A.matrix[0][1] / 3;
    A.matrix[1][0] = A.matrix[1][0] / 3;
    A.matrix[1][1] = A.matrix[1][1] / 3;

    A.matrix[0][2] = A.matrix[0][2] * (-2. / 3);
    A.matrix[1][2] = A.matrix[1][2] * (-2. / 3);
    A.matrix[2][0] = A.matrix[2][0] * (-2. / 3);
    A.matrix[2][1] = A.matrix[2][1] * (-2. / 3);

    A.matrix[2][2] = A.matrix[2][2] * (4. / 3);

    return su3_test(A);
  }

  friend su3_test operator+(const su3_test &A, const su3_test &B) {
    su3_test C;
    complex_test a;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        C.matrix[i][j] = A.matrix[i][j] + B.matrix[i][j];
      }
    }
    return C;
  }
  friend su3_test operator-(const su3_test &A, const su3_test &B) {
    su3_test C;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        C.matrix[i][j] = A.matrix[i][j] - B.matrix[i][j];
      }
    }
    return C;
  }
  friend su3_test operator*(const double &x, const su3_test &A) {
    su3_test C;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        C.matrix[i][j] = x * A.matrix[i][j];
      }
    }
    return C;
  }
  friend su3_test operator*(const su3_test &A, const double &x) {
    su3_test C;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        C.matrix[i][j] = x * A.matrix[i][j];
      }
    }
    return C;
  }
  friend su3_test operator*(const complex_test &x, const su3_test &A) {
    su3_test C;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        C.matrix[i][j] = x * A.matrix[i][j];
      }
    }
    return C;
  }
  friend su3_test operator*(const su3_test &A, const complex_test &x) {
    su3_test C;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        C.matrix[i][j] = x * A.matrix[i][j];
      }
    }
    return C;
  }

  friend su3_test operator*(const su3_test &A, const su3_test &B) {
    su3_test C;
    complex_test a;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        a = complex_test(0, 0);
        for (int k = 0; k < 3; k++) {
          a = a + A.matrix[i][k] * B.matrix[k][j];
        }
        C.matrix[i][j] = a;
      }
    }
    return C;
  }

  friend su3_test operator^(const su3_test &A, const su3_test &B) {
    su3_test C;
    complex_test a;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        a = complex_test(0, 0);
        for (int k = 0; k < 3; k++) {
          a = a + (A.matrix[i][k] ^ B.matrix[j][k]);
        }
        C.matrix[i][j] = a;
      }
    }
    return C;
  }

  friend su3_test operator%(const su3_test &A, const su3_test &B) {
    su3_test C;
    complex_test a;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        a = complex_test(0, 0);
        for (int k = 0; k < 3; k++) {
          a = a + (A.matrix[k][i] % B.matrix[k][j]);
        }
        C.matrix[i][j] = a;
      }
    }
    return C;
  }

  friend std::ostream &operator<<(std::ostream &os, const su3_test &A) {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        os << "(" << A.matrix[i][j].real << ", " << A.matrix[i][j].imag << ") ";
      }
      os << std::endl;
    }
    return os;
  }
};

class su3_eigen {
public:
  Eigen::Matrix3cd matrix;
  su3_eigen() {
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (i != j)
          matrix(i, j) = std::complex<double>(0, 0);
        else
          matrix(i, j) = std::complex<double>(1., 0);
      }
    }
  }

  su3_eigen(const Eigen::Matrix3cd &B) { matrix = B; }

  double tr() {
    return (matrix(0, 0).real() + matrix(1, 1).real() + matrix(2, 2).real()) /
           3;
  }

  std::complex<double> tr_complex() {
    return (matrix(0, 0) + matrix(1, 1) + matrix(2, 2)) / 3.;
  }

  // double multiply_conj_tr(const su3_eigen &B) {
  //   double trace = 0;
  //   for (int i = 0; i < 3; i++) {
  //     for (int k = 0; k < 3; k++) {
  //       trace += matrix(i, k).real() * B.matrix(i, k).real() +
  //                matrix(i, k).imag() * B.matrix(i, k).imag();
  //     }
  //   }
  //   return trace / 3;
  // }

  double multiply_conj_tr(const su3_eigen &B) const {
    return (matrix * B.matrix.adjoint()).trace().real() / 3;
  }

  double multiply_tr(const su3_eigen &B) {
    double trace = 0;
    for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
        trace += matrix(i, k).real() * B.matrix(k, i).real() -
                 matrix(i, k).imag() * B.matrix(k, i).imag();
      }
    }
    return trace / 3;
  }

  double multiply_conj_tr_adjoint(const su3_eigen &B) {
    std::complex<double> trace = std::complex<double>(0, 0);
    for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
        trace += matrix(i, k) * std::conj(B.matrix(i, k));
      }
    }
    return std::norm(trace) - 1;
  }

  // su3_eigen conj() const { return su3_eigen(matrix.adjoint().eval()); }

  friend su3_eigen operator+(const su3_eigen &A, const su3_eigen &B) {
    return su3_eigen(A.matrix + B.matrix);
  }
  friend su3_eigen operator-(const su3_eigen &A, const su3_eigen &B) {
    return su3_eigen(A.matrix - B.matrix);
  }
  friend su3_eigen operator*(const double &x, const su3_eigen &A) {
    return su3_eigen(x * A.matrix);
  }
  friend su3_eigen operator*(const su3_eigen &A, const double &x) {
    return su3_eigen(A.matrix * x);
  }
  friend su3_eigen operator*(const std::complex<double> &x,
                             const su3_eigen &A) {
    return su3_eigen(x * A.matrix);
  }
  friend su3_eigen operator*(const su3_eigen &A,
                             const std::complex<double> &x) {
    return su3_eigen(A.matrix * x);
  }

  friend su3_eigen operator*(const su3_eigen &A, const su3_eigen &B) {
    return su3_eigen(A.matrix * B.matrix);
  }

  friend su3_eigen operator^(const su3_eigen &A, const su3_eigen &B) {
    return su3_eigen(A.matrix * B.matrix.adjoint().eval());
  }

  friend su3_eigen operator%(const su3_eigen &A, const su3_eigen &B) {
    return su3_eigen(A.matrix.adjoint().eval() * B.matrix);
  }

  friend std::ostream &operator<<(std::ostream &os, const su3_eigen &A) {
    os << A.matrix << std::endl;
    return os;
  }
};