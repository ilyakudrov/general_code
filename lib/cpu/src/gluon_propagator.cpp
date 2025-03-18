#include "../include/indexing.h"
#include "../include/matrix.h"

#include <Eigen/Dense>
#include <array>
#include <complex>
#include <vector>

#pragma omp declare reduction(                                                 \
        arr_double_plus : std::array<std::complex<double>, 144> : std::        \
        transform(omp_out.begin(), omp_out.end(), omp_in.begin(),              \
                      omp_out.begin(), std::plus<std::complex<double>>()))     \
    initializer(omp_priv = decltype(omp_orig)())

Eigen::Matrix2cd get_vector_potential_matrix(su2 &a) {
  Eigen::Matrix2cd U;
  U(0, 0) = std::complex<double>(a.a0, a.a3);
  U(0, 1) = std::complex<double>(a.a2, a.a1);
  U(1, 0) = std::complex<double>(-a.a2, a.a1);
  U(1, 1) = std::complex<double>(a.a0, -a.a3);
  Eigen::Matrix2cd A = U - U.adjoint();
  return std::complex<double>(0, -0.5) *
         (A - 0.5 * A.trace() * Eigen::Matrix2cd::Identity(2, 2));
}

std::array<Eigen::Matrix2cd, 3> get_pauli_matrices() {
  std::array<Eigen::Matrix2cd, 3> pauli_matrices;
  pauli_matrices[0] << 0, 1, 1, 0;
  pauli_matrices[1] << 0, std::complex<double>(0, -1),
      std::complex<double>(0, 1), 0;
  pauli_matrices[2] << 1, 0, 0, -1;
  return pauli_matrices;
}

std::array<double, 3> get_vector_potential_coefficients(
    su2 &a, std::array<Eigen::Matrix2cd, 3> &pauli_matrices) {
  std::array<double, 3> coefficients;
  for (int i = 0; i < 3; i++) {
    coefficients[i] =
        -2 *
        (get_vector_potential_matrix(a) * pauli_matrices[i]).trace().real();
  }
  return coefficients;
}

std::vector<std::array<double, 12>>
get_vector_potential(std::vector<su2> &conf) {
  std::vector<std::array<double, 12>> vector_potential(x_size * y_size *
                                                       z_size * t_size);
  std::array<Eigen::Matrix2cd, 3> pauli_matrices = get_pauli_matrices();
  std::array<int, 4> lat_coord;
  std::array<int, 4> lat_dim = {x_size, y_size, z_size, t_size};
  std::array<std::array<double, 3>, 4> vector_potential_tmp1;
  std::array<double, 12> vector_potential_tmp2;
  int index;
  // #pragma omp parallel for private(vector_potential_tmp, index,                  \
//                                      vector_potential_tmp2)                    \
//     firstprivate(pauli_matrices)
  for (int t = 0; t < lat_dim[3]; t++) {
    for (int z = 0; z < lat_dim[2]; z++) {
      for (int y = 0; y < lat_dim[1]; y++) {
        for (int x = 0; x < lat_dim[0]; x++) {
          for (int mu = 0; mu < 4; mu++) {
            lat_coord = {x, y, z, t};
            vector_potential_tmp1[mu] = get_vector_potential_coefficients(
                conf[get_index_matrix(lat_coord, mu)], pauli_matrices);
          }
          index = get_index_site(lat_coord);
          for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 3; j++) {
              vector_potential_tmp2[i * 3 + j] = vector_potential_tmp1[i][j];
            }
          }
          vector_potential[get_index_site(lat_coord)] = vector_potential_tmp2;
        }
      }
    }
  }
  return vector_potential;
}

std::vector<std::complex<double>>
get_furier_coefficients(std::array<double, 4> &momenta) {
  std::vector<std::complex<double>> furier_coefficients(x_size * y_size *
                                                        z_size * t_size);
  std::array<int, 4> lat_coord;
  std::array<int, 4> lat_dim = {x_size, y_size, z_size, t_size};
  for (int t = 0; t < lat_dim[3]; t++) {
    for (int z = 0; z < lat_dim[2]; z++) {
      for (int y = 0; y < lat_dim[1]; y++) {
        for (int x = 0; x < lat_dim[0]; x++) {
          lat_coord = {x, y, z, t};
          furier_coefficients[get_index_site(lat_coord)] =
              std::exp(std::complex<double>(0, 1) *
                       (momenta[0] * lat_coord[0] + momenta[1] * lat_coord[1] +
                        momenta[2] * lat_coord[2] + momenta[3] * lat_coord[3]));
        }
      }
    }
  }
  return furier_coefficients;
}

std::array<std::complex<double>, 144> calculate_gluon_propagator(
    std::vector<std::array<double, 12>> &vector_potential,
    std::vector<std::complex<double>> &furier_coefficients) {
  std::array<std::complex<double>, 144> gluon_propagator;
  std::array<std::complex<double>, 12> furier_coefficient_x;
#pragma omp parallel
  for (int i = 0; i < vector_potential.size(); i++) {
    // std::cout << i << std::endl;
    for (int m = 0; m < 12; m++) {
      furier_coefficient_x[m] =
          vector_potential[i][m] * std::conj(furier_coefficients[i]);
    }
#pragma omp for firstprivate(furier_coefficient_x)                             \
    reduction(arr_double_plus : gluon_propagator)
    for (int j = 0; j < vector_potential.size(); j++) {
      for (int m = 0; m < 12; m++) {
        for (int n = 0; n < 12; n++) {
          gluon_propagator[m * 12 + n] += furier_coefficient_x[m] *
                                          vector_potential[j][n] *
                                          furier_coefficients[j];
        }
      }
    }
  }
  int lattice_volume = x_size * y_size * z_size * t_size;
  for (int i = 0; i < 144; i++) {
    gluon_propagator[i] /= lattice_volume;
  }
  return gluon_propagator;
}