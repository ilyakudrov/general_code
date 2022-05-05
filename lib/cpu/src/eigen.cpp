#include "../../../lib/cpu/include/eigen.h"

#define SPACE_ITER_START                                                       \
  for (int t = 0; t < t_size; t++) {                                           \
    for (int z = 0; z < z_size; z++) {                                         \
      for (int y = 0; y < y_size; y++) {                                       \
        for (int x = 0; x < x_size; x++) {                                     \
          link.go(x, y, z, t);                                                 \
          link.update(0);                                                      \
          link.update(1);                                                      \
          link.update(2);                                                      \
          link.update(3);

#define SPACE_ITER_END                                                         \
  }                                                                            \
  }                                                                            \
  }                                                                            \
  }

std::vector<su2> make_matrix_staggered(const std::vector<su2> &conf,
                                       double mu_q) {
  int matrix_size = x_size * y_size * z_size * t_size * 8;
  link1 link(x_size, y_size, z_size, t_size);

  std::vector<su2> matrix(matrix_size);

  int delta_4;
  int sign;
  int border_sign;

  int place_matrix;

  SPACE_ITER_START

  place_matrix = link.place * 2;

  for (int mu = 0; mu < 4; mu++) {
    if (mu == 3)
      delta_4 = 1;
    else
      delta_4 = 0;

    link.move_dir(mu);

    sign = eta_sign(mu, link);

    // positive direction
    border_sign = 1;
    if (mu == 3 && link.coordinate[mu] == link.lattice_size[mu] - 1)
      border_sign = -1;

    matrix[place_matrix + mu * 2] =
        exp(mu_q * delta_4) / 2 * border_sign * sign * conf[link.place + mu];

    // negative direction
    border_sign = 1;
    if (mu == 3 && link.coordinate[mu] == 0)
      border_sign = -1;

    link.move(mu, -1);

    matrix[place_matrix + mu * 2 + 1] = -exp(-mu_q * delta_4) / 2 *
                                        border_sign * sign *
                                        (conf[link.place + mu]).conj();

    link.move(mu, 1);
  }

  SPACE_ITER_END

  return matrix;
}

std::vector<complex>
matrix_multiplication_staggered(std::vector<su2> &matrix,
                                const std::vector<complex> &vec_input) {
  link1 link(x_size, y_size, z_size, t_size);
  int vec_size = x_size * y_size * z_size * t_size * 2;

  std::vector<complex> vec_output(vec_size);
  for (int i = 0; i < vec_size; i++) {
    vec_output[i].re = 0;
    vec_output[i].im = 0;
  }

  int place_vector_center;
  int place_matrix;

  SPACE_ITER_START

  place_vector_center = link.place / 2;
  place_matrix = link.place * 2;

  for (int mu = 0; mu < 4; mu++) {

    // positive direction
    link.move(mu, 1);

    matrix_multiplication_su2(&vec_input[link.place / 2],
                              &vec_output[place_vector_center],
                              matrix[place_matrix + 2 * mu]);

    // negative direction
    link.move(mu, -2);

    matrix_multiplication_su2(&vec_input[link.place / 2],
                              &vec_output[place_vector_center],
                              matrix[place_matrix + 2 * mu + 1]);

    link.move(mu, 1);
  }

  SPACE_ITER_END

  return vec_output;
}

void matrix_multiplication_su2(const complex *vec_input, complex *vec_output,
                               su2 &A) {
  vec_output[0].re += A.a0 * vec_input[0].re - A.a3 * vec_input[0].im +
                      A.a2 * vec_input[1].re - A.a1 * vec_input[1].im;
  vec_output[0].im += A.a0 * vec_input[0].im + A.a3 * vec_input[0].re +
                      A.a2 * vec_input[1].im + A.a1 * vec_input[1].re;
  vec_output[1].re += -A.a2 * vec_input[0].re - A.a1 * vec_input[0].im +
                      A.a0 * vec_input[1].re + A.a3 * vec_input[1].im;
  vec_output[1].im += -A.a2 * vec_input[0].im + A.a1 * vec_input[0].re +
                      A.a0 * vec_input[1].im - A.a3 * vec_input[1].re;
}

double eta_sign(int mu, link1 &link) {
  int n = 0;
  for (int i = 0; i < mu; i++) {
    n += (link.coordinate[i]);
  }
  return 1 - (n % 2) * 2;
}