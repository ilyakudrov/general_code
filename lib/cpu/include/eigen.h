#pragma once

#include "data.h"
#include "link.h"

template <class DataPattern>
std::vector<su2>
make_matrix_staggered(const Data::LatticeData<DataPattern, su2> &conf,
                      double mu_q) {
  DataPattern data_pattern(conf.lat_dim);
  int matrix_size = data_pattern.get_lattice_size * 8;
  std::vector<su2> matrix(matrix_size);
  int delta_4;
  int sign;
  int border_sign;
  int place_matrix;
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          place_matrix = data_pattern.get_index_site * 8;
          for (int mu = 0; mu < 4; mu++) {
            if (mu == 3)
              delta_4 = 1;
            else
              delta_4 = 0;
            sign = eta_sign(mu, data_pattern);
            // positive direction
            border_sign = 1;
            if (mu == 3 &&
                data_pattern.lat_coord[mu] == data_pattern.lat_dim[mu] - 1)
              border_sign = -1;
            matrix[place_matrix + mu * 2] =
                exp(mu_q * delta_4) / 2 * border_sign * sign *
                conf[data_pattern.get_index_link(mu)];
            // negative direction
            border_sign = 1;
            if (mu == 3 && data_pattern.lat_coord[mu] == 0)
              border_sign = -1;
            data_pattern.move_backward(1, mu);
            matrix[place_matrix + mu * 2 + 1] =
                -exp(-mu_q * delta_4) / 2 * border_sign * sign *
                (conf[data_pattern.get_index_link(mu)]).conj();
            data_pattern.move_forward(1, mu);
          }
        }
      }
    }
  }
  return matrix;
}

template <class DataPattern>
std::vector<std::complex<double>> matrix_multiplication_staggered(
    std::vector<su2> &matrix,
    const std::vector<std::complex<double>> &vec_input,
    DataPattern &data_pattern) {
  int vec_size = data_pattern.get_lattice_size * 2;
  std::vector<std::complex<double>> vec_output(vec_size);
  int place_vector_center;
  int place_matrix;
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          place_vector_center = data_pattern.get_index_site * 2;
          place_matrix = data_pattern.get_index_site * 8;
          for (int mu = 0; mu < 4; mu++) {
            // positive direction
            data_pattern.move_forward(1, mu);
            matrix_multiplication_su2(
                &vec_input[data_pattern.get_index_site * 2],
                &vec_output[place_vector_center],
                matrix[place_matrix + 2 * mu]);
            // negative direction
            data_pattern.move_backward(2, mu);
            matrix_multiplication_su2(
                &vec_input[data_pattern.get_index_site * 2],
                &vec_output[place_vector_center],
                matrix[place_matrix + 2 * mu + 1]);
            data_pattern.move_forward(1, mu);
          }
        }
      }
    }
  }
  return vec_output;
}

template <class DataPattern>
inline double eta_sign(int mu, DataPattern &data_pattern) {
  int n = 0;
  for (int i = 0; i < mu; i++) {
    n += (data_pattern.lat_coord[i]);
  }
  return 1 - (n % 2) * 2;
}

void matrix_multiplication_su2(const std::complex<double> *vec_input,
                               std::complex<double> *vec_output, su2 &A);