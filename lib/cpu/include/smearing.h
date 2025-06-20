#pragma once

#include "data.h"
#include "indexing.h"
#include "link.h"

#include <map>
#include <unordered_map>
#include <vector>

template <class DataPattern>
su2 staples_first(const std::vector<su2> &vec, DataPattern data_pattern, int nu,
                  int eta) {
  su2 A;
  su2 B;
  A = vec[data_pattern.get_index_link(eta)];
  data_pattern.move_forward(1, eta);
  A = A * vec[data_pattern.get_index_link(nu)];
  data_pattern.move_forward(1, nu);
  data_pattern.move_backward(1, eta);
  A = A ^ vec[data_pattern.get_index_link(eta)];
  data_pattern.move_backward(1, nu);
  data_pattern.move_backward(1, eta);
  B = vec[data_pattern.get_index_link(eta)].conj();
  B = B * vec[data_pattern.get_index_link(nu)];
  data_pattern.move_backward(1, nu);
  B = B * vec[data_pattern.get_index_link(eta)];
  data_pattern.move_forward(1, eta);
  data_pattern.move_backward(1, nu);
  return (A + B);
}

template <class DataPattern>
su2 stout_omega(const Data::LatticeData<DataPattern, su2> &conf,
                DataPattern data_pattern, int nu, double rho) {
  su2 A;
  su2 B(0., 0., 0., 0.);
  A = conf.array[data_pattern.get_index_link(nu)].inverse();
  for (int i = 0; i < 4; i++) {
    if (i != nu) {
      B = B + staples_first(conf.array, data_pattern, nu, i);
    }
  }
  A = B * A * rho;
  return A;
}

template <class DataPattern>
su2 stout_factor(const Data::LatticeData<DataPattern, su2> &conf,
                 DataPattern data_pattern, int nu, double rho) {
  su2 A;
  su2 B;
  su2 C;
  su2 C1;
  C = stout_omega(conf, data_pattern, rho);
  C1.a0 = C.a0;
  C1.a1 = -C.a1;
  C1.a2 = -C.a2;
  C1.a3 = -C.a3;
  A = (-1.) / 2 * (C1 + (-1.) * C + (-1.) / 2 * A * (C1.tr() + (-1.) * C.tr()));
  B.a0 = exp(A.a0) * cos(powf(A.a1 * A.a1 + A.a2 * A.a2 + A.a3 * A.a3, 0.5));
  B.a1 = exp(A.a0) * A.a1 *
         sin(powf(A.a1 * A.a1 + A.a2 * A.a2 + A.a3 * A.a3, 0.5)) /
         powf(A.a1 * A.a1 + A.a2 * A.a2 + A.a3 * A.a3, 0.5);
  B.a2 = exp(A.a0) * A.a2 *
         sin(powf(A.a1 * A.a1 + A.a2 * A.a2 + A.a3 * A.a3, 0.5)) /
         powf(A.a1 * A.a1 + A.a2 * A.a2 + A.a3 * A.a3, 0.5);
  B.a3 = exp(A.a0) * A.a3 *
         sin(powf(A.a1 * A.a1 + A.a2 * A.a2 + A.a3 * A.a3, 0.5)) /
         powf(A.a1 * A.a1 + A.a2 * A.a2 + A.a3 * A.a3, 0.5);
  return B;
}

template <class DataPattern>
void smearing_stout(Data::LatticeData<DataPattern, su2> &conf, double rho) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<su2> vec(conf.array.size());
#pragma omp parallel for collapse(4) firstprivate(data_pattern, rho)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          for (int mu = 0; mu < 4; mu++) {
            vec[data_pattern.get_index_link(mu)] =
                stout_factor(conf, data_pattern, mu, rho) *
                conf.array[data_pattern.get_index_link(mu)];
          }
        }
      }
    }
  }
  conf.array = vec;
}

template <class DataPattern, class MatrixType>
void smearing_APE(Data::LatticeData<DataPattern, MatrixType> &conf,
                  double alpha) {
  Data::LatticeData<DataPattern, MatrixType> smeared = conf;
  DataPattern data_pattern(conf.lat_dim);
  MatrixType staple;
  MatrixType staple_sum;
  int index;
#pragma omp parallel for collapse(4) private(index, staple, staple_sum)        \
    firstprivate(data_pattern, alpha)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          for (int mu = 0; mu < 3; mu++) {
            staple_sum = (1 - alpha) * smeared[data_pattern.get_index_link(mu)];
            for (int nu = 0; nu < 3; nu++) {
              if (mu != nu) {
                index = data_pattern.get_index_link(nu);
                data_pattern.move_forward(1, nu);
                staple = conf[index] * conf[data_pattern.get_index_link(mu)];
                data_pattern.move_backward(1, nu);
                data_pattern.move_forward(1, mu);
                staple_sum =
                    staple_sum + ((alpha / 4) * staple ^
                                  conf[data_pattern.get_index_link(nu)]);
                data_pattern.move_backward(1, nu);
                data_pattern.move_backward(1, mu);
                staple = conf[data_pattern.get_index_link(nu)] %
                         conf[data_pattern.get_index_link(mu)];
                data_pattern.move_forward(1, mu);
                staple_sum =
                    staple_sum + ((alpha / 4) * staple *
                                  conf[data_pattern.get_index_link(nu)]);
                data_pattern.move_forward(1, nu);
                data_pattern.move_backward(1, mu);
              }
            }
            smeared[data_pattern.get_index_link(mu)] = staple_sum.proj();
          }
        }
      }
    }
  }
  conf = smeared;
}

inline su3_angles
multiply_proj_su3_angles(const std::array<su3_angles, 4> &staple,
                         const su3_angles &link, const double &a,
                         const double &b) {
  su3_angles A;
  for (int i = 0; i < 3; i++) {
    A.matrix[i] =
        atan2(a * sin(link.matrix[i]) +
                  b * (sin(staple[0].matrix[i]) + sin(staple[1].matrix[i]) +
                       sin(staple[2].matrix[i]) + sin(staple[3].matrix[i])),
              a * cos(link.matrix[i]) +
                  b * (cos(staple[0].matrix[i]) + cos(staple[1].matrix[i]) +
                       cos(staple[2].matrix[i]) + cos(staple[3].matrix[i])));
  }
  double sum = 0;
  for (int i = 0; i < 3; i++) {
    sum += A.matrix[i];
  }
  while (sum >= M_PI) {
    sum -= 2 * M_PI;
  }
  while (sum < -M_PI) {
    sum += 2 * M_PI;
  }
  for (int i = 0; i < 3; i++) {
    A.matrix[i] = A.matrix[i] - sum / 3;
  }
  return A;
}

template <class DataPattern>
void smearing_APE(Data::LatticeData<DataPattern, su3_angles> &conf,
                  double alpha) {
  Data::LatticeData<DataPattern, su3_angles> smeared = conf;
  DataPattern data_pattern(conf.lat_dim);
  int index1;
  int index2;
#pragma omp parallel for collapse(4) private(index1, index2)                   \
    firstprivate(data_pattern, alpha)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          for (int mu = 0; mu < 3; mu++) {
            int count = 0;
            std::array<su3_angles, 4> tmp;
            for (int nu = 0; nu < 3; nu++) {
              if (mu != nu) {
                index1 = data_pattern.get_index_link(nu);
                data_pattern.move_forward(1, nu);
                index2 = data_pattern.get_index_link(mu);
                data_pattern.move_backward(1, nu);
                data_pattern.move_forward(1, mu);
                tmp[count] = (conf[index1] * conf[index2]) ^
                             conf[data_pattern.get_index_link(nu)];
                count++;
                data_pattern.move_backward(1, nu);
                data_pattern.move_backward(1, mu);
                index1 = data_pattern.get_index_link(nu);
                index2 = data_pattern.get_index_link(mu);
                data_pattern.move_forward(1, mu);
                tmp[count] = (conf[index1] % conf[index2]) *
                             conf[data_pattern.get_index_link(nu)];
                count++;
                data_pattern.move_forward(1, nu);
                data_pattern.move_backward(1, mu);
              }
            }
            index1 = data_pattern.get_index_link(mu);
            smeared[index1] = multiply_proj_su3_angles(tmp, conf[index1],
                                                       1 - alpha, alpha / 4);
          }
        }
      }
    }
  }
  conf = smeared;
}

template <class DataPattern, class MatrixType>
void smearing_APE_2d(std::vector<MatrixType> &conf1,
                     std::vector<MatrixType> &conf2, DataPattern &data_pattern,
                     int mu, int nu, double alpha) {
  std::vector<MatrixType> smeared1(conf1.size());
  std::vector<MatrixType> smeared2(conf2.size());
  MatrixType staple;
  MatrixType A;
  MatrixType B;
  int index1;
  int index2;
  int index3;
#pragma omp parallel for collapse(4) private(index1, index2, index3, staple,   \
                                                 A, B)                         \
    firstprivate(data_pattern, mu, nu, alpha)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index1 = data_pattern.get_index_site();
          data_pattern.move_forward(1, nu);
          index2 = data_pattern.get_index_site();
          data_pattern.move_backward(1, nu);
          data_pattern.move_forward(1, mu);
          index3 = data_pattern.get_index_site();
          A = conf2[index3] ^ conf1[index2];
          B = conf2[index1] ^ A;
          A = conf1[index1] * A;
          data_pattern.move_backward(1, nu);
          data_pattern.move_backward(1, mu);
          index2 = data_pattern.get_index_site();
          staple = conf2[index2] % conf1[index2];
          data_pattern.move_forward(1, mu);
          smeared1[index1] =
              ((1 - alpha) * conf1[index1] +
               (alpha / 2) *
                   (B + (staple * conf2[data_pattern.get_index_site()])))
                  .proj();
          data_pattern.move_forward(1, nu);
          data_pattern.move_backward(2, mu);
          index2 = data_pattern.get_index_site();
          B = conf1[index2] % conf2[index2];
          data_pattern.move_forward(1, nu);
          smeared2[index1] =
              ((1 - alpha) * conf2[index1] +
               (alpha / 2) * (A + (B * conf1[data_pattern.get_index_site()])))
                  .proj();
        }
      }
    }
  }
  conf1 = std::move(smeared1);
  conf2 = std::move(smeared2);
}

std::vector<std::vector<std::tuple<int, int, int>>> make_indices3();
std::vector<std::tuple<int, int, int>>
indices_to_delete(std::vector<std::vector<std::tuple<int, int, int>>> &indices,
                  std::vector<std::tuple<int, int, int>> &deleted_indices,
                  int dir);

template <class DataPattern, class MatrixType>
void smearing_plane_HYP(std::vector<MatrixType> &smeared,
                        const Data::LatticeData<DataPattern, MatrixType> &conf,
                        int mu, int nu, double alpha, double divisor) {
  MatrixType staple;
  int index1;
  int index2;
  DataPattern data_pattern(conf.lat_dim);
#pragma omp parallel for collapse(4) private(index1, index2, staple)           \
    firstprivate(data_pattern, mu, nu, alpha, divisor)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index1 = data_pattern.get_index_site();
          index2 = data_pattern.get_index_link(nu);
          data_pattern.move_forward(1, nu);
          staple = conf[index2] * conf[data_pattern.get_index_link(mu)];
          data_pattern.move_backward(1, nu);
          data_pattern.move_forward(1, mu);
          smeared[index1] =
              smeared[index1] + ((alpha / divisor) * staple ^
                                 conf[data_pattern.get_index_link(nu)]);
          data_pattern.move_backward(1, nu);
          data_pattern.move_backward(1, mu);
          staple = conf[data_pattern.get_index_link(nu)] %
                   conf[data_pattern.get_index_link(mu)];
          data_pattern.move_forward(1, mu);
          smeared[index1] =
              smeared[index1] + ((alpha / divisor) * staple *
                                 conf[data_pattern.get_index_link(nu)]);
        }
      }
    }
  }
}

template <class DataPattern, class MatrixType>
void smearing_plane_HYP(std::vector<MatrixType> &smeared,
                        const std::vector<MatrixType> &conf_mu,
                        const std::vector<MatrixType> &conf_nu,
                        DataPattern &data_pattern, int mu, int nu, double alpha,
                        double divisor) {
  MatrixType staple;
  int index;
#pragma omp parallel for collapse(4) private(index, staple)                    \
    firstprivate(data_pattern, mu, nu, alpha, divisor)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_index_site();
          data_pattern.move_forward(1, nu);
          staple = conf_nu[index] * conf_mu[data_pattern.get_index_site()];
          data_pattern.move_backward(1, nu);
          data_pattern.move_forward(1, mu);
          smeared[index] =
              smeared[index] + ((alpha / divisor) * staple ^
                                conf_nu[data_pattern.get_index_site()]);
          data_pattern.move_backward(1, nu);
          data_pattern.move_backward(1, mu);
          staple = conf_nu[data_pattern.get_index_site()] %
                   conf_mu[data_pattern.get_index_site()];
          data_pattern.move_forward(1, mu);
          smeared[index] =
              smeared[index] + ((alpha / divisor) * staple *
                                conf_nu[data_pattern.get_index_site()]);
        }
      }
    }
  }
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType>
HYP_initialize_vector(const Data::LatticeData<DataPattern, MatrixType> &conf,
                      int mu, double alpha) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> links(data_pattern.get_lattice_size());
#pragma omp parallel for collapse(4) firstprivate(data_pattern, mu, alpha)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          links[data_pattern.get_index_site()] =
              (1 - alpha) * conf[data_pattern.get_index_link(mu)];
        }
      }
    }
  }
  return links;
}

template <class DataPattern, class MatrixType>
void make_step1(std::vector<std::vector<MatrixType>> &links1,
                const Data::LatticeData<DataPattern, MatrixType> &conf,
                std::vector<std::tuple<int, int, int>> &indices,
                std::map<std::tuple<int, int, int>, int> &indices_map,
                double alpha) {
  for (int i = 0; i < indices.size(); i++) {
    if (indices_map.count(indices[i]) == 0) {
      indices_map[indices[i]] = links1.size();
      links1.push_back(
          HYP_initialize_vector(conf, std::get<0>(indices[i]), alpha));
      for (int eta = 0; eta < 3; eta++) {
        if (eta != std::get<0>(indices[i]) && eta != std::get<1>(indices[i]) &&
            eta != std::get<2>(indices[i])) {
          smearing_plane_HYP(links1.back(), conf, std::get<0>(indices[i]), eta,
                             alpha, 2);
        }
      }
#pragma omp parallel for
      for (int i = 0; i < links1.back().size(); i++) {
        links1.back()[i] = links1.back()[i].proj();
      }
    }
  }
}

template <class DataPattern, class MatrixType>
void make_step2(std::vector<std::vector<MatrixType>> &links2,
                std::vector<std::vector<MatrixType>> &links1,
                const Data::LatticeData<DataPattern, MatrixType> &conf,
                std::map<std::tuple<int, int, int>, int> &indices_map1, int nu,
                double alpha2) {
  links2[0] = HYP_initialize_vector(conf, 3, alpha2);
  links2[1] = HYP_initialize_vector(conf, nu, alpha2);
  int index = 0;
  DataPattern data_pattern(conf.lat_dim);
  for (int rho = 0; rho < 3; rho++) {
    if (rho != nu) {
      if (nu < rho) {
        index = indices_map1[std::tuple<int, int, int>(3, nu, rho)];
      } else {
        index = indices_map1[std::tuple<int, int, int>(3, rho, nu)];
      }
      smearing_plane_HYP(
          links2[0], links1[index],
          links1[indices_map1[std::tuple<int, int, int>(rho, nu, 3)]],
          data_pattern, 3, rho, alpha2, 4);
    }
  }
#pragma omp parallel for
  for (int i = 0; i < links2[0].size(); i++) {
    links2[0][i] = links2[0][i].proj();
  }
  for (int rho = 0; rho < 3; rho++) {
    if (rho != nu) {
      smearing_plane_HYP(
          links2[1],
          links1[indices_map1[std::tuple<int, int, int>(nu, rho, 3)]],
          links1[indices_map1[std::tuple<int, int, int>(rho, nu, 3)]],
          data_pattern, nu, rho, alpha2, 4);
    }
  }
#pragma omp parallel for
  for (int i = 0; i < links2[1].size(); i++) {
    links2[1][i] = links2[1][i].proj();
  }
}

template <class T>
inline void
delete_vectors(std::vector<std::tuple<int, int, int>> &indices_to_delete,
               std::map<std::tuple<int, int, int>, int> indices_map,
               std::vector<std::vector<T>> &links) {
  for (int i = 0; i < indices_to_delete.size(); i++) {
    links[indices_map[indices_to_delete[i]]].clear();
    links[indices_map[indices_to_delete[i]]].shrink_to_fit();
  }
}

template <class DataPattern, class MatrixType>
void smearing_HYP(Data::LatticeData<DataPattern, MatrixType> &conf,
                  double alpha1, double alpha2, double alpha3) {
  std::vector<std::vector<std::tuple<int, int, int>>> indices3 =
      make_indices3();
  std::vector<MatrixType> smeared = HYP_initialize_vector(conf, 3, alpha1);
  std::vector<std::vector<MatrixType>> links2(2);
  std::vector<std::vector<MatrixType>> links1;
  std::map<std::tuple<int, int, int>, int> indices_map3;
  std::vector<std::tuple<int, int, int>> deleted_indices;
  std::vector<std::tuple<int, int, int>> delete_indices;
  DataPattern data_pattern(conf.lat_dim);
  for (int nu3 = 0; nu3 < 3; nu3++) {
    make_step1(links1, conf, indices3[nu3], indices_map3, alpha3);
    make_step2(links2, links1, conf, indices_map3, nu3, alpha2);
    delete_indices = indices_to_delete(indices3, deleted_indices, nu3);
    delete_vectors(delete_indices, indices_map3, links1);
    smearing_plane_HYP(smeared, links2[0], links2[1], data_pattern, 3, nu3,
                       alpha1, 6);
  }
#pragma omp parallel for
  for (int i = 0; i < smeared.size(); i++) {
    smeared[i] = smeared[i].proj();
  }
#pragma omp parallel for collapse(4) firstprivate(data_pattern)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          conf[data_pattern.get_index_link(3)] =
              smeared[data_pattern.get_index_site()];
        }
      }
    }
  }
}

std::vector<su2> smearing_stout(Data::data<su2> &conf, double rho);
template <class T>
std::vector<std::vector<T>> separate_smearing(std::vector<T> &conf);
template <class T> void smearing_APE(std::vector<T> &conf, double alpha);
template <class T>
void smearing_APE_2d(std::vector<T> &conf1, std::vector<T> &conf2, int mu,
                     int nu, double alpha);
template <class T>
void smearing_APE_2d_test(std::vector<T> &conf1, std::vector<T> &conf2, int mu,
                          int nu, double alpha);
template <class T>
void smearing_HYP(std::vector<T> &array, double alpha1, double alpha2,
                  double alpha3);
