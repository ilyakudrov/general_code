#include "../include/gluon_propagator.h"
#include "../include/data.h"
#include "../include/indexing.h"
#include "../include/matrix.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <complex>
#include <omp.h>
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

#pragma omp declare reduction(                                                 \
        arr_double_plus_4 : std::array<std::complex<double>, 4> : std::        \
        transform(omp_out.begin(), omp_out.end(), omp_in.begin(),              \
                      omp_out.begin(), std::plus<std::complex<double>>()))     \
    initializer(omp_priv = decltype(omp_orig)())

Eigen::Matrix2cd get_vector_potential_matrix(const su2 &a) {
  Eigen::Matrix2cd U;
  U(0, 0) = std::complex<double>(a.a0, a.a3);
  U(0, 1) = std::complex<double>(a.a2, a.a1);
  U(1, 0) = std::complex<double>(-a.a2, a.a1);
  U(1, 1) = std::complex<double>(a.a0, -a.a3);
  Eigen::Matrix2cd A = U - U.adjoint();
  return std::complex<double>(0, -0.5) *
         (A - 0.5 * A.trace() * Eigen::Matrix2cd::Identity(2, 2));
}

Eigen::Matrix2cd get_vector_potential_matrix(const abelian &a) {
  Eigen::Matrix2cd U;
  U(0, 0) = std::complex<double>(cos(a.phi), sin(a.phi));
  U(0, 1) = std::complex<double>(0, 0);
  U(1, 0) = std::complex<double>(0, 0);
  U(1, 1) = std::complex<double>(cos(a.phi), -sin(a.phi));
  Eigen::Matrix2cd A = U - U.adjoint();
  return std::complex<double>(0, -0.5) *
         (A - 0.5 * A.trace() * Eigen::Matrix2cd::Identity(2, 2));
}

std::array<double, 3> get_vector_potential_coefficients(
    const su2 &a, const std::array<Eigen::Matrix2cd, 3> &pauli_matrices) {
  std::array<double, 3> coefficients;
  for (int i = 0; i < 3; i++) {
    coefficients[i] =
        (get_vector_potential_matrix(a) * pauli_matrices[i]).trace().real() *
        0.5;
  }
  return coefficients;
}

double get_vector_potential_coefficients(const abelian &a) {
  Eigen::Matrix2cd A;
  A << 1, 0, 0, -1;
  return (get_vector_potential_matrix(a) * A).trace().real() * 0.5;
}

std::array<Eigen::Matrix2cd, 3> get_pauli_matrices() {
  std::array<Eigen::Matrix2cd, 3> pauli_matrices;
  pauli_matrices[0] << 0, 1, 1, 0;
  pauli_matrices[1] << 0, std::complex<double>(0, -1),
      std::complex<double>(0, 1), 0;
  pauli_matrices[2] << 1, 0, 0, -1;
  return pauli_matrices;
}

std::vector<std::array<double, 12>> get_vector_potential(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf) {
  DataPatternLexicographical data_pattern(conf.lat_dim);
  std::vector<std::array<double, 12>> vector_potential(
      data_pattern.get_lattice_size());
  std::array<Eigen::Matrix2cd, 3> pauli_matrices = get_pauli_matrices();
  std::array<std::array<double, 3>, 4> vector_potential_tmp1;
  std::array<double, 12> vector_potential_tmp2;
  int index;
  // #pragma omp parallel for private(vector_potential_tmp, index,                  \
//                                      vector_potential_tmp2)                    \
//     firstprivate(pauli_matrices)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          for (int mu = 0; mu < 4; mu++) {
            data_pattern.lat_coord = {x, y, z, t};
            vector_potential_tmp1[mu] = get_vector_potential_coefficients(
                conf[data_pattern.get_index_link(mu)], pauli_matrices);
          }
          index = data_pattern.get_index_site();
          for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 3; j++) {
              vector_potential_tmp2[i * 3 + j] = vector_potential_tmp1[i][j];
            }
          }
          vector_potential[data_pattern.get_index_site()] =
              vector_potential_tmp2;
        }
      }
    }
  }
  return vector_potential;
}

std::vector<std::array<double, 4>> get_vector_potential(
    const Data::LatticeData<DataPatternLexicographical, abelian> &conf) {
  DataPatternLexicographical data_pattern(conf.lat_dim);
  std::vector<std::array<double, 4>> vector_potential(
      data_pattern.get_lattice_size());
  std::array<double, 4> vector_potential_tmp;
  int index;
  // #pragma omp parallel for private(vector_potential_tmp, index,                  \
  //                                    vector_potential_tmp2)                    \
  //   firstprivate(pauli_matrices)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          for (int mu = 0; mu < 4; mu++) {
            data_pattern.lat_coord = {x, y, z, t};
            vector_potential_tmp[mu] = get_vector_potential_coefficients(
                conf[data_pattern.get_index_link(mu)]);
          }
          vector_potential[data_pattern.get_index_site()] =
              vector_potential_tmp;
        }
      }
    }
  }
  return vector_potential;
}

std::vector<std::complex<double>>
get_furier_coefficients(const std::array<double, 4> &momenta,
                        DataPatternLexicographical &data_pattern) {
  std::vector<std::complex<double>> furier_coefficients(
      data_pattern.get_lattice_size());
  double power;
#pragma omp parallel for collapse(3) private(power)                            \
    firstprivate(data_pattern, momenta)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          power = momenta[0] * data_pattern.lat_coord[0] +
                  momenta[1] * data_pattern.lat_coord[1] +
                  momenta[2] * data_pattern.lat_coord[2] +
                  momenta[3] * data_pattern.lat_coord[3];
          furier_coefficients[data_pattern.get_index_site()] =
              std::complex<double>(cos(power), sin(power));
        }
      }
    }
  }
  return furier_coefficients;
}

std::array<std::complex<double>, 144> calculate_gluon_propagator(
    const std::vector<std::array<double, 12>> &vector_potential,
    const std::vector<std::complex<double>> &furier_coefficients,
    double multiplier, DataPatternLexicographical &data_pattern) {
  std::array<std::complex<double>, 144> gluon_propagator;
  std::array<std::complex<double>, 12> vector_potential_sum;
  double a = multiplier / data_pattern.get_lattice_size();
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

std::array<std::complex<double>, 144> calculate_gluon_propagator_group(
    const std::vector<std::array<double, 12>> &vector_potential,
    const std::vector<std::array<double, 4>> &momenta, double multiplier,
    DataPatternLexicographical &data_pattern) {
  std::array<std::complex<double>, 144> gluon_propagator;
  double a = multiplier / data_pattern.get_lattice_size() / momenta.size();
  std::vector<std::complex<double>> furier_coefficients;
  for (int p = 0; p < momenta.size(); p++) {
    furier_coefficients = get_furier_coefficients(momenta[p], data_pattern);
    std::array<std::complex<double>, 12> vector_potential_sum;
    for (int i = 0; i < vector_potential.size(); i++) {
      for (int m = 0; m < 12; m++) {
        vector_potential_sum[m] +=
            vector_potential[i][m] * furier_coefficients[i];
      }
    }
    for (int m = 0; m < 12; m++) {
      for (int n = 0; n < 12; n++) {
        gluon_propagator[m * 12 + n] +=
            vector_potential_sum[m] * std::conj(vector_potential_sum[n]) * a;
      }
    }
  }
  return gluon_propagator;
}

std::array<std::complex<double>, 4> calculate_gluon_propagator_diagonal_group(
    const std::vector<std::array<double, 4>> &vector_potential,
    const std::vector<std::array<double, 4>> &momenta, double multiplier,
    DataPatternLexicographical &data_pattern) {
  std::array<std::complex<double>, 4> gluon_propagator;
  double a = multiplier / data_pattern.get_lattice_size() / momenta.size();
  std::vector<std::complex<double>> furier_coefficients;
  for (int p = 0; p < momenta.size(); p++) {
    furier_coefficients = get_furier_coefficients(momenta[p], data_pattern);
    std::array<std::complex<double>, 4> vector_potential_sum;
#pragma omp parallel for reduction(arr_double_plus_4 : vector_potential_sum)
    for (int i = 0; i < vector_potential.size(); i++) {
      for (int m = 0; m < 4; m++) {
        vector_potential_sum[m] +=
            vector_potential[i][m] * furier_coefficients[i];
      }
    }
    for (int m = 0; m < 4; m++) {
      gluon_propagator[m] +=
          vector_potential_sum[m] * std::conj(vector_potential_sum[m]) * a;
    }
  }
  return gluon_propagator;
}

std::array<std::complex<double>, 144> calculate_gluon_propagator_lattice(
    const std::vector<std::array<double, 12>> &vector_potential,
    const std::array<double, 4> &momenta, double multiplier,
    DataPatternLexicographical &data_pattern) {
  std::array<std::complex<double>, 144> gluon_propagator;
  std::array<std::complex<double>, 12> vector_potential_sum;
  double a = multiplier / data_pattern.get_lattice_size();
  int index;
  std::complex<double> furier_coefficient;
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_index_site();
          double power = momenta[0] * data_pattern.lat_coord[0] +
                         momenta[1] * data_pattern.lat_coord[1] +
                         momenta[2] * data_pattern.lat_coord[2] +
                         momenta[3] * data_pattern.lat_coord[3];
          furier_coefficient = std::complex<double>(cos(power), sin(power));
          for (int m = 0; m < 12; m++) {
            vector_potential_sum[m] +=
                vector_potential[index][m] * furier_coefficient;
          }
        }
      }
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
    const std::vector<std::array<double, 12>> &vector_potential,
    const std::vector<std::complex<double>> &furier_coefficients,
    double multiplier, DataPatternLexicographical &data_pattern) {
  std::array<std::complex<double>, 12> gluon_propagator;
  // #pragma omp parallel for reduction(arr_double_plus_12 : gluon_propagator)      \
//     schedule(dynamic)
  for (int j = 0; j < vector_potential.size(); j++) {
    for (int m = 0; m < 12; m++) {
      gluon_propagator[m] += vector_potential[j][m] * furier_coefficients[j];
    }
  }
  int lattice_volume = data_pattern.get_lattice_size();
  for (int m = 0; m < 12; m++) {
    gluon_propagator[m] = gluon_propagator[m] * std::conj(gluon_propagator[m]) *
                          (multiplier / lattice_volume);
  }
  return gluon_propagator;
}

std::vector<std::array<double, 3>> get_vector_potential_longitudinal(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf) {
  DataPatternLexicographical data_pattern(conf.lat_dim);
  std::vector<std::array<double, 3>> vector_potential(
      data_pattern.get_lattice_size());
  std::array<Eigen::Matrix2cd, 3> pauli_matrices = get_pauli_matrices();
  int mu = 3;
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          vector_potential[data_pattern.get_index_site()] =
              get_vector_potential_coefficients(
                  conf[data_pattern.get_index_link(mu)], pauli_matrices);
        }
      }
    }
  }
  return vector_potential;
}

double calculate_gluon_propagator_longitudinal_zero_momentum(
    const std::vector<std::array<double, 3>> &vector_potential,
    double multiplier, DataPatternLexicographical &data_pattern) {
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
  int lattice_volume = data_pattern.get_lattice_size();
  return result / lattice_volume * multiplier;
}

std::complex<double> calculate_gluon_propagator_single(
    const std::vector<std::array<double, 12>> &vector_potential,
    const std::vector<std::complex<double>> &furier_coefficients,
    double multiplier, int m, int n, DataPatternLexicographical &data_pattern) {
  std::complex<double> gluon_propagator{0, 0};
  for (int i = 0; i < vector_potential.size(); i++) {
    gluon_propagator += vector_potential[i][m] * furier_coefficients[i];
  }
  int lattice_volume = data_pattern.get_lattice_size();
  gluon_propagator = gluon_propagator * std::conj(gluon_propagator) *
                     (multiplier / lattice_volume);
  return gluon_propagator;
}

std::vector<std::vector<std::array<double, 4>>> generate_momenta(int Ns,
                                                                 int Nt) {
  double multiplyer_t = 2 * M_PI / Nt;
  double multiplyer_s = 2 * M_PI / Ns;
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