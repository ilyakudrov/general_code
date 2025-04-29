#include "../include/indexing.h"
#include "../include/matrix.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <complex>
#include <vector>

#pragma omp declare reduction(                                                 \
        arr_double_plus_144 : std::array<std::complex<double>, 144> : std::    \
        transform(omp_out.begin(), omp_out.end(), omp_in.begin(),              \
                      omp_out.begin(), std::plus<std::complex<double>>()))     \
    initializer(omp_priv = decltype(omp_orig)())

#pragma omp declare reduction(                                                 \
        arr_double_plus_12 : std::array<std::complex<double>, 12> : std::      \
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
        (get_vector_potential_matrix(a) * pauli_matrices[i]).trace().real() *
        0.5;
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
    std::vector<std::complex<double>> &furier_coefficients, double multiplier) {
  std::array<std::complex<double>, 144> gluon_propagator;
  std::array<std::complex<double>, 12> vector_potential_sum;
  double a = multiplier / (x_size * y_size * z_size * t_size);
  for (int i = 0; i < vector_potential.size(); i++) {
    for (int m = 0; m < 12; m++) {
      vector_potential_sum[m] +=
          vector_potential[i][m] * furier_coefficients[i];
    }
  }
  for (int m = 0; m < 12; m++) {
    for (int n = 0; n < 12; n++) {
      gluon_propagator[m * 12 + n] =
          vector_potential_sum[m] * std::conj(vector_potential_sum[n]) * a;
    }
  }
  return gluon_propagator;
}

std::array<std::complex<double>, 12> calculate_gluon_propagator_diagonal(
    std::vector<std::array<double, 12>> &vector_potential,
    std::vector<std::complex<double>> &furier_coefficients, double multiplier) {
  std::array<std::complex<double>, 12> gluon_propagator;
  // #pragma omp parallel for reduction(arr_double_plus_12 : gluon_propagator)      \
//     schedule(dynamic)
  for (int j = 0; j < vector_potential.size(); j++) {
    for (int m = 0; m < 12; m++) {
      gluon_propagator[m] += vector_potential[j][m] * furier_coefficients[j];
    }
  }
  int lattice_volume = x_size * y_size * z_size * t_size;
  for (int m = 0; m < 12; m++) {
    gluon_propagator[m] = gluon_propagator[m] * std::conj(gluon_propagator[m]) *
                          (multiplier / lattice_volume);
  }
  return gluon_propagator;
}

std::vector<std::array<double, 3>>
get_vector_potential_longitudinal(std::vector<su2> &conf) {
  std::vector<std::array<double, 3>> vector_potential(x_size * y_size * z_size *
                                                      t_size);
  std::array<Eigen::Matrix2cd, 3> pauli_matrices = get_pauli_matrices();
  std::array<int, 4> lat_coord;
  std::array<int, 4> lat_dim = {x_size, y_size, z_size, t_size};
  int mu = 3;
  for (int t = 0; t < lat_dim[3]; t++) {
    for (int z = 0; z < lat_dim[2]; z++) {
      for (int y = 0; y < lat_dim[1]; y++) {
        for (int x = 0; x < lat_dim[0]; x++) {
          lat_coord = {x, y, z, t};
          vector_potential[get_index_site(lat_coord)] =
              get_vector_potential_coefficients(
                  conf[get_index_matrix(lat_coord, mu)], pauli_matrices);
        }
      }
    }
  }
  return vector_potential;
}

double calculate_gluon_propagator_longitudinal_zero_momentum(
    std::vector<std::array<double, 3>> &vector_potential, double multiplier) {
  double result = 0;
  result = 0;
  for (int j = 0; j < 3; j++) {
    double tmp = 0;
#pragma omp parallel for reduction(+ : tmp) schedule(dynamic)
    for (int i = 0; i < vector_potential.size(); i++) {
      tmp += vector_potential[i][j];
    }
    result += tmp * tmp / 3;
  }
  int lattice_volume = x_size * y_size * z_size * t_size;
  return result / lattice_volume * multiplier;
}

std::complex<double> calculate_gluon_propagator_single(
    std::vector<std::array<double, 12>> &vector_potential,
    std::vector<std::complex<double>> &furier_coefficients, double multiplier,
    int m, int n) {
  std::complex<double> gluon_propagator{0, 0};
  for (int i = 0; i < vector_potential.size(); i++) {
    gluon_propagator += vector_potential[i][m] * furier_coefficients[i];
  }
  int lattice_volume = x_size * y_size * z_size * t_size;
  gluon_propagator = gluon_propagator * std::conj(gluon_propagator) *
                     (multiplier / lattice_volume);
  return gluon_propagator;
}

std::vector<std::vector<std::array<double, 4>>> generate_momenta(int Ns,
                                                                 int Nt) {
  double multiplyer_t = 2 * M_PI / Ns;
  double multiplyer_s = 2 * M_PI / Nt;
  std::vector<std::vector<std::array<double, 4>>> momenta;
  std::array<int, 3> momentum1;
  std::array<int, 3> momentum2;
  for (int qx = 0; qx <= Ns / 2; qx++) {
    for (int qy = qx; qy <= Ns / 2; qy++) {
      for (int qz = qy; qz <= Ns / 2; qz++) {
        momentum1 = {qx, qy, qz};
        std::vector<std::array<double, 4>> momenta_tmp;
        momenta_tmp.push_back({momentum1[0] * multiplyer_s,
                               momentum1[1] * multiplyer_s,
                               momentum1[2] * multiplyer_s, 0});
        momentum2 = momentum1;
        while (std::next_permutation(momentum2.begin(), momentum2.end())) {
          momenta_tmp.push_back({momentum2[0] * multiplyer_s,
                                 momentum2[1] * multiplyer_s,
                                 momentum2[2] * multiplyer_s, 0});
        }
        momentum2 = momentum1;
        while (std::prev_permutation(momentum2.begin(), momentum2.end())) {
          momenta_tmp.push_back({momentum2[0] * multiplyer_s,
                                 momentum2[1] * multiplyer_s,
                                 momentum2[2] * multiplyer_s, 0});
        }
        if (qx != qy && qx != Ns / 2 && qx != 0 && qy != 0 && qz != 0) {
          momentum2 = momentum1;
          momentum2[0] = -momentum2[0];
          momenta_tmp.push_back({momentum2[0] * multiplyer_s,
                                 momentum2[1] * multiplyer_s,
                                 momentum2[2] * multiplyer_s, 0});
          while (std::next_permutation(momentum2.begin(), momentum2.end())) {
            momenta_tmp.push_back({momentum2[0] * multiplyer_s,
                                   momentum2[1] * multiplyer_s,
                                   momentum2[2] * multiplyer_s, 0});
          }
          momentum2 = momentum1;
          momentum2[0] = -momentum2[0];
          while (std::prev_permutation(momentum2.begin(), momentum2.end())) {
            momenta_tmp.push_back({momentum2[0] * multiplyer_s,
                                   momentum2[1] * multiplyer_s,
                                   momentum2[2] * multiplyer_s, 0});
          }
        }
        if (qy != qz && qy != Ns / 2 && qy != 0 && qz != 0) {
          momentum2 = momentum1;
          momentum2[1] = -momentum2[1];
          momenta_tmp.push_back({momentum2[0] * multiplyer_s,
                                 momentum2[1] * multiplyer_s,
                                 momentum2[2] * multiplyer_s, 0});
          while (std::next_permutation(momentum2.begin(), momentum2.end())) {
            momenta_tmp.push_back({momentum2[0] * multiplyer_s,
                                   momentum2[1] * multiplyer_s,
                                   momentum2[2] * multiplyer_s, 0});
          }
          momentum2 = momentum1;
          momentum2[1] = -momentum2[1];
          while (std::prev_permutation(momentum2.begin(), momentum2.end())) {
            momenta_tmp.push_back({momentum2[0] * multiplyer_s,
                                   momentum2[1] * multiplyer_s,
                                   momentum2[2] * multiplyer_s, 0});
          }
        }
        if (qz != Ns / 2 && qy != 0 && qz != 0) {
          momentum2 = momentum1;
          momentum2[2] = -momentum2[2];
          momenta_tmp.push_back({momentum2[0] * multiplyer_s,
                                 momentum2[1] * multiplyer_s,
                                 momentum2[2] * multiplyer_s, 0});
          while (std::next_permutation(momentum2.begin(), momentum2.end())) {
            momenta_tmp.push_back({momentum2[0] * multiplyer_s,
                                   momentum2[1] * multiplyer_s,
                                   momentum2[2] * multiplyer_s, 0});
          }
          momentum2 = momentum1;
          momentum2[2] = -momentum2[2];
          while (std::prev_permutation(momentum2.begin(), momentum2.end())) {
            momenta_tmp.push_back({momentum2[0] * multiplyer_s,
                                   momentum2[1] * multiplyer_s,
                                   momentum2[2] * multiplyer_s, 0});
          }
        }
        momenta.push_back(momenta_tmp);
      }
    }
  }
  return momenta;
}