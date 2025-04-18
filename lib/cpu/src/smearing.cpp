#include "../include/data.h"
#include "../include/indexing.h"
#include "../include/link.h"
#include "../include/matrix.h"

#include "stdlib.h"
#include <algorithm>
#include <map>
#include <omp.h>
#include <utility>
#include <vector>

#define data_size 4 * x_size *y_size *z_size *t_size
#define PLACE1_LINK_NODIR                                                      \
  (link.coordinate[3]) * link.lattice_size[0] * link.lattice_size[1] *         \
          link.lattice_size[2] +                                               \
      (link.coordinate[2]) * link.lattice_size[0] * link.lattice_size[1] +     \
      (link.coordinate[1]) * link.lattice_size[0] + (link.coordinate[0])
#define PLACE4_LINK_DIR                                                        \
  (link.coordinate[3]) * 4 * link.lattice_size[0] * link.lattice_size[1] *     \
          link.lattice_size[2] +                                               \
      (link.coordinate[2]) * 4 * link.lattice_size[0] * link.lattice_size[1] + \
      (link.coordinate[1]) * 4 * link.lattice_size[0] +                        \
      (link.coordinate[0]) * 4 + link.direction

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

#define SPACE_ITER_START_3D                                                    \
  for (int z = 0; z < z_size; z++) {                                           \
    for (int y = 0; y < y_size; y++) {                                         \
      for (int x = 0; x < x_size; x++) {                                       \
        link.go(x, y, z, t);                                                   \
        link.update(0);                                                        \
        link.update(1);                                                        \
        link.update(2);                                                        \
        link.update(3);

#define SPACE_ITER_END_3D                                                      \
  }                                                                            \
  }                                                                            \
  }

template <class T>
T staples_first(const std::vector<T> &vec, link1 &link, int eta) {
  T A;
  T B;
  int dir = link.direction;
  A = vec[link.place + eta];
  link.move(eta, 1);
  A = A * vec[link.place + dir];
  link.move(dir, 1);
  link.move(eta, -1);
  A = A ^ vec[link.place + eta];
  link.move(dir, -1);
  link.move(eta, -1);
  B = vec[link.place + eta].conj();
  B = B * vec[link.place + dir];
  link.move(dir, 1);
  B = B * vec[link.place + eta];
  link.move(eta, 1);
  link.move(dir, -1);

  return (A + B);
}

su2 stout_omega(Data::data<su2> &conf, link1 &link, double rho) {
  int dir = link.direction;
  su2 A;
  su2 B(0., 0., 0., 0.);
  A = conf.array[PLACE4_LINK_DIR].inverse();
  for (int i = 0; i < 4; i++) {
    if (i != dir) {
      B = B + staples_first(conf.array, link, i);
    }
  }
  A = B * A * rho;
  return A;
}

su2 stout_factor(Data::data<su2> &conf, link1 &link, double rho) {
  su2 A;
  su2 B;
  su2 C;
  su2 C1;
  C = stout_omega(conf, link, rho);
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

std::vector<su2> smearing_stout(Data::data<su2> &conf, double rho) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<su2> vec(data_size);
  vec = conf.array;
  SPACE_ITER_START
  for (int i = 0; i < 4; i++) {
    link.move_dir(i);
    vec[PLACE4_LINK_DIR] =
        stout_factor(conf, link, rho) * conf.array[PLACE4_LINK_DIR];
  }
  SPACE_ITER_END
  return vec;
}

template <class T>
std::vector<std::vector<T>> separate_smearing(std::vector<T> &conf) {
  int data_size1 = x_size * y_size * z_size * t_size;
  std::vector<std::vector<T>> result(4, std::vector<T>(data_size1));

  link1 link(x_size, y_size, z_size, t_size);

  for (int mu = 0; mu < 4; ++mu) {

    SPACE_ITER_START

    result[mu][link.place / 4] = conf[link.place + mu];

    SPACE_ITER_END
  }

  return result;
}

template <class T> void smearing_APE(std::vector<T> &conf, double alpha) {
  std::vector<T> smeared = conf;
  std::array<int, 4> lat_dim = {x_size, y_size, z_size, t_size};
  std::array<int, 4> lat_coord;
  T staple;
  T staple_sum;
  int index;
#pragma omp parallel for collapse(4) private(                                  \
        lat_coord, index, staple, staple_sum) firstprivate(lat_dim, alpha)
  for (int t = 0; t < lat_dim[3]; t++) {
    for (int z = 0; z < lat_dim[2]; z++) {
      for (int y = 0; y < lat_dim[1]; y++) {
        for (int x = 0; x < lat_dim[0]; x++) {
          lat_coord = {x, y, z, t};
          for (int mu = 0; mu < 3; mu++) {
            staple_sum = (1 - alpha) * smeared[get_index_matrix(lat_coord, mu)];
            for (int nu = 0; nu < 3; nu++) {
              if (mu != nu) {
                index = get_index_matrix(lat_coord, nu);
                lat_coord[nu] = (lat_coord[nu] + 1) % lat_dim[nu];
                staple = conf[index] * conf[get_index_matrix(lat_coord, mu)];
                lat_coord[nu] = (lat_coord[nu] + lat_dim[nu] - 1) % lat_dim[nu];
                lat_coord[mu] = (lat_coord[mu] + 1) % lat_dim[mu];
                staple_sum =
                    staple_sum + ((alpha / 4) * staple ^
                                  conf[get_index_matrix(lat_coord, nu)]);
                lat_coord[nu] = (lat_coord[nu] + lat_dim[nu] - 1) % lat_dim[nu];
                lat_coord[mu] = (lat_coord[mu] + lat_dim[mu] - 1) % lat_dim[mu];
                staple = conf[get_index_matrix(lat_coord, nu)] %
                         conf[get_index_matrix(lat_coord, mu)];
                lat_coord[mu] = (lat_coord[mu] + 1) % lat_dim[mu];
                staple_sum =
                    staple_sum + ((alpha / 4) * staple *
                                  conf[get_index_matrix(lat_coord, nu)]);
                lat_coord[nu] = (lat_coord[nu] + 1) % lat_dim[nu];
                lat_coord[mu] = (lat_coord[mu] + lat_dim[mu] - 1) % lat_dim[mu];
              }
            }
            smeared[get_index_matrix(lat_coord, mu)] = staple_sum.proj();
          }
        }
      }
    }
  }
  conf = smeared;
}

su3_angles multiply_proj_su3_angles(const std::array<su3_angles, 4> &staple,
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

template <> void smearing_APE(std::vector<su3_angles> &conf, double alpha) {
  std::vector<su3_angles> smeared = conf;
  std::array<int, 4> lat_dim = {x_size, y_size, z_size, t_size};
  std::array<int, 4> lat_coord;
  int index1;
  int index2;
#pragma omp parallel for collapse(4) private(lat_coord, index1, index2)        \
    firstprivate(lat_dim, alpha)
  for (int t = 0; t < lat_dim[3]; t++) {
    for (int z = 0; z < lat_dim[2]; z++) {
      for (int y = 0; y < lat_dim[1]; y++) {
        for (int x = 0; x < lat_dim[0]; x++) {
          lat_coord = {x, y, z, t};
          for (int mu = 0; mu < 3; mu++) {
            int count = 0;
            std::array<su3_angles, 4> tmp;
            for (int nu = 0; nu < 3; nu++) {
              if (mu != nu) {
                index1 = get_index_matrix(lat_coord, nu);
                lat_coord[nu] = (lat_coord[nu] + 1) % lat_dim[nu];
                index2 = get_index_matrix(lat_coord, mu);
                lat_coord[nu] = (lat_coord[nu] + lat_dim[nu] - 1) % lat_dim[nu];
                lat_coord[mu] = (lat_coord[mu] + 1) % lat_dim[mu];
                tmp[count] = (conf[index1] * conf[index2]) ^
                             conf[get_index_matrix(lat_coord, nu)];
                count++;
                lat_coord[nu] = (lat_coord[nu] + lat_dim[nu] - 1) % lat_dim[nu];
                lat_coord[mu] = (lat_coord[mu] + lat_dim[mu] - 1) % lat_dim[mu];
                index1 = get_index_matrix(lat_coord, nu);
                index2 = get_index_matrix(lat_coord, mu);
                lat_coord[mu] = (lat_coord[mu] + 1) % lat_dim[mu];
                tmp[count] = (conf[index1] % conf[index2]) *
                             conf[get_index_matrix(lat_coord, nu)];
                count++;
                lat_coord[nu] = (lat_coord[nu] + 1) % lat_dim[nu];
                lat_coord[mu] = (lat_coord[mu] + lat_dim[mu] - 1) % lat_dim[mu];
              }
            }
            index1 = get_index_matrix(lat_coord, mu);
            smeared[index1] = multiply_proj_su3_angles(tmp, conf[index1],
                                                       1 - alpha, alpha / 4);
          }
        }
      }
    }
  }
  conf = smeared;
}

template <class T>
void smearing_APE_2d(std::vector<T> &conf1, std::vector<T> &conf2, int mu,
                     int nu, double alpha) {
  std::vector<T> smeared1(conf1.size());
  std::vector<T> smeared2(conf2.size());
  std::array<int, 4> lat_dim = {x_size, y_size, z_size, t_size};
  std::array<int, 4> lat_coord;
  T staple;
  T A;
  T B;
  int index1;
  int index2;
  int index3;
#pragma omp parallel for collapse(4) private(lat_coord, index1, index2,        \
                                                 index3, staple, A, B)         \
    firstprivate(lat_dim, mu, nu, alpha)
  for (int t = 0; t < lat_dim[3]; t++) {
    for (int z = 0; z < lat_dim[2]; z++) {
      for (int y = 0; y < lat_dim[1]; y++) {
        for (int x = 0; x < lat_dim[0]; x++) {
          lat_coord = {x, y, z, t};
          index1 = get_index_site(lat_coord);
          lat_coord[nu] = (lat_coord[nu] + 1) % lat_dim[nu];
          index2 = get_index_site(lat_coord);
          lat_coord[nu] = (lat_coord[nu] + lat_dim[nu] - 1) % lat_dim[nu];
          lat_coord[mu] = (lat_coord[mu] + 1) % lat_dim[mu];
          index3 = get_index_site(lat_coord);
          A = conf2[index3] ^ conf1[index2];
          B = conf2[index1] ^ A;
          A = conf1[index1] * A;
          lat_coord[nu] = (lat_coord[nu] + lat_dim[nu] - 1) % lat_dim[nu];
          lat_coord[mu] = (lat_coord[mu] + lat_dim[mu] - 1) % lat_dim[mu];
          index2 = get_index_site(lat_coord);
          staple = conf2[index2] % conf1[index2];
          lat_coord[mu] = (lat_coord[mu] + 1) % lat_dim[mu];
          smeared1[index1] =
              ((1 - alpha) * conf1[index1] +
               (alpha / 2) * (B + (staple * conf2[get_index_site(lat_coord)])))
                  .proj();
          lat_coord[nu] = (lat_coord[nu] + 1) % lat_dim[nu];
          lat_coord[mu] = (lat_coord[mu] + lat_dim[mu] - 2) % lat_dim[mu];
          index2 = get_index_site(lat_coord);
          B = conf1[index2] % conf2[index2];
          lat_coord[nu] = (lat_coord[nu] + 1) % lat_dim[nu];
          smeared2[index1] =
              ((1 - alpha) * conf2[index1] +
               (alpha / 2) * (A + (B * conf1[get_index_site(lat_coord)])))
                  .proj();
        }
      }
    }
  }
  conf1 = std::move(smeared1);
  conf2 = std::move(smeared2);
}

std::vector<std::vector<std::tuple<int, int, int>>> make_indices3() {
  std::vector<std::vector<std::tuple<int, int, int>>> indices(3);
  for (int nu = 0; nu < 3; nu++) {
    // for V_{nu;3}
    for (int rho = 0; rho < 3; rho++) {
      if (rho != nu) {
        indices[nu].push_back(std::tuple<int, int, int>(rho, nu, 3));
        indices[nu].push_back(std::tuple<int, int, int>(nu, rho, 3));
        if (rho < nu) {
          indices[nu].push_back(std::tuple<int, int, int>(3, rho, nu));
        } else {
          indices[nu].push_back(std::tuple<int, int, int>(3, nu, rho));
        }
      }
    }
  }
  return indices;
}

bool if_not_contained(std::vector<std::tuple<int, int, int>> &indices,
                      std::tuple<int, int, int> index) {
  bool exists = false;
  for (int i = 0; i < indices.size(); i++) {
    exists = exists || (indices[i] == index);
  }
  return !exists;
}

std::vector<std::tuple<int, int, int>>
indices_to_delete(std::vector<std::vector<std::tuple<int, int, int>>> &indices,
                  std::vector<std::tuple<int, int, int>> &deleted_indices,
                  int dir) {
  std::vector<std::tuple<int, int, int>> indices_to_delete;
  bool add_index = false;
  // up to dir direction of last step of HYP try to find
  // indices, which will not be used later
  for (int i = 0; i <= dir; i++) {
    for (int j = 0; j < indices[i].size(); j++) {
      add_index = if_not_contained(indices_to_delete, indices[i][j]) &&
                  if_not_contained(indices_to_delete, indices[i][j]);
      for (int k = dir; k <= 3; k++) {
        add_index = add_index && if_not_contained(indices[k], indices[i][j]);
      }
    }
  }
  return indices_to_delete;
}

template <class T>
std::vector<T> HYP_initialize_vector(std::vector<T> &conf, double alpha) {
  std::vector<T> links(conf.size());
  for (int i = 0; i < conf.size(); i++) {
    links[i] = (1 - alpha) * conf[i];
  }
  return links;
}

template <class T>
void delete_vectors(std::vector<std::tuple<int, int, int>> &indices_to_delete,
                    std::map<std::tuple<int, int, int>, int> indices_map,
                    std::vector<std::vector<T>> &links) {
  for (int i = 0; i < indices_to_delete.size(); i++) {
    links[indices_map[indices_to_delete[i]]].clear();
    links[indices_map[indices_to_delete[i]]].shrink_to_fit();
  }
}

template <class T>
void smearing_plane_HYP(std::vector<T> &smeared, std::vector<T> &conf_mu,
                        std::vector<T> &conf_nu, int mu, int nu, double alpha,
                        double divisor) {
  T staple;
  int index;
  std::array<int, 4> lat_dim = {x_size, y_size, z_size, t_size};
  std::array<int, 4> lat_coord;
#pragma omp parallel for collapse(4) private(lat_coord, index, staple)         \
    firstprivate(lat_dim)
  for (int t = 0; t < lat_dim[3]; t++) {
    for (int z = 0; z < lat_dim[2]; z++) {
      for (int y = 0; y < lat_dim[1]; y++) {
        for (int x = 0; x < lat_dim[0]; x++) {
          lat_coord = {x, y, z, t};
          index = get_index_site(lat_coord);
          lat_coord[nu] = (lat_coord[nu] + 1) % lat_dim[nu];
          staple = conf_nu[index] * conf_mu[get_index_site(lat_coord)];
          lat_coord[nu] = (lat_coord[nu] + lat_dim[nu] - 1) % lat_dim[nu];
          lat_coord[mu] = (lat_coord[mu] + 1) % lat_dim[mu];
          smeared[index] =
              smeared[index] +
              ((alpha / divisor) * staple ^ conf_nu[get_index_site(lat_coord)]);
          lat_coord[nu] = (lat_coord[nu] + lat_dim[nu] - 1) % lat_dim[nu];
          lat_coord[mu] = (lat_coord[mu] + lat_dim[mu] - 1) % lat_dim[mu];
          staple = conf_nu[get_index_site(lat_coord)] %
                   conf_mu[get_index_site(lat_coord)];
          lat_coord[mu] = (lat_coord[mu] + 1) % lat_dim[mu];
          smeared[index] =
              smeared[index] +
              ((alpha / divisor) * staple * conf_nu[get_index_site(lat_coord)]);
        }
      }
    }
  }
}

template <class T>
void make_step1(std::vector<std::vector<T>> &links1,
                std::vector<std::vector<T>> &conf,
                std::vector<std::tuple<int, int, int>> &indices,
                std::map<std::tuple<int, int, int>, int> &indices_map,
                double alpha) {
  std::vector<int> steps = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};
  for (int i = 0; i < indices.size(); i++) {
    if (indices_map.count(indices[i]) == 0) {
      indices_map[indices[i]] = links1.size();
      links1.push_back(
          HYP_initialize_vector(conf[std::get<0>(indices[i])], alpha));
      for (int eta = 0; eta < 3; eta++) {
        if (eta != std::get<0>(indices[i]) && eta != std::get<1>(indices[i]) &&
            eta != std::get<2>(indices[i])) {
          smearing_plane_HYP(links1.back(), conf[std::get<0>(indices[i])],
                             conf[eta], std::get<0>(indices[i]), eta, alpha, 2);
        }
      }
#pragma omp parallel for
      for (int i = 0; i < links1.back().size(); i++) {
        links1.back()[i] = links1.back()[i].proj();
      }
    }
  }
}

template <class T>
void make_step2(std::vector<std::vector<T>> &links2,
                std::vector<std::vector<T>> &links1,
                std::vector<std::vector<T>> &conf,
                std::map<std::tuple<int, int, int>, int> &indices_map1, int nu,
                double alpha2) {
  std::vector<int> steps = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};
  links2[0] = HYP_initialize_vector(conf[3], alpha2);
  links2[1] = HYP_initialize_vector(conf[nu], alpha2);
  int index = 0;
  for (int rho = 0; rho < 3; rho++) {
    if (rho != nu) {
      if (nu < rho) {
        index = indices_map1[std::tuple<int, int, int>(3, nu, rho)];
      } else {
        index = indices_map1[std::tuple<int, int, int>(3, rho, nu)];
      }
      smearing_plane_HYP(
          links2[0], links1[index],
          links1[indices_map1[std::tuple<int, int, int>(rho, nu, 3)]], 3, rho,
          alpha2, 4);
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
          links1[indices_map1[std::tuple<int, int, int>(rho, nu, 3)]], nu, rho,
          alpha2, 4);
    }
  }
#pragma omp parallel for
  for (int i = 0; i < links2[1].size(); i++) {
    links2[1][i] = links2[1][i].proj();
  }
}

template <class T>
void smearing_HYP(std::vector<T> &array, double alpha1, double alpha2,
                  double alpha3) {
  std::vector<int> steps = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};
  std::vector<std::vector<T>> conf = separate_smearing(array);
  array.clear();
  array.shrink_to_fit();
  std::vector<std::vector<std::tuple<int, int, int>>> indices3 =
      make_indices3();
  std::vector<T> smeared = HYP_initialize_vector(conf[3], alpha1);
  std::vector<std::vector<T>> links2(2);
  std::vector<std::vector<T>> links1;
  std::map<std::tuple<int, int, int>, int> indices_map3;
  std::vector<std::tuple<int, int, int>> deleted_indices;
  std::vector<std::tuple<int, int, int>> delete_indices;
  for (int nu3 = 0; nu3 < 3; nu3++) {
    make_step1(links1, conf, indices3[nu3], indices_map3, alpha3);
    make_step2(links2, links1, conf, indices_map3, nu3, alpha2);
    delete_indices = indices_to_delete(indices3, deleted_indices, nu3);
    delete_vectors(delete_indices, indices_map3, links1);
    smearing_plane_HYP(smeared, links2[0], links2[1], 3, nu3, alpha1, 6);
  }
#pragma omp parallel for
  for (int i = 0; i < smeared.size(); i++) {
    smeared[i] = smeared[i].proj();
  }
  conf[3] = smeared;
  for (int i = 0; i < links2.size(); i++) {
    links2[i].clear();
    links2[i].shrink_to_fit();
  }
  for (int i = 0; i < links1.size(); i++) {
    links1[i].clear();
    links1[i].shrink_to_fit();
  }
  smeared.clear();
  smeared.shrink_to_fit();
  array = std::vector<T>(4 * conf[3].size());
  std::array<int, 4> lat_coord;
#pragma omp parallel for collapse(4) private(lat_coord)
  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
          lat_coord = {x, y, z, t};
          for (int mu = 0; mu < 4; mu++) {
            array[get_index_matrix(lat_coord, mu)] =
                conf[mu][get_index_site(lat_coord)];
          }
        }
      }
    }
  }
}

// specialications

// su2
template std::vector<std::vector<su2>>
separate_smearing(std::vector<su2> &conf);
template void smearing_APE(std::vector<su2> &conf, double alpha);
template void smearing_APE_2d(std::vector<su2> &conf1, std::vector<su2> &conf2,
                              int mu, int nu, double alpha);
template void smearing_HYP(std::vector<su2> &array, double alpha1,
                           double alpha2, double alpha3);

// abelian
template std::vector<std::vector<abelian>>
separate_smearing(std::vector<abelian> &conf);
template void smearing_APE(std::vector<abelian> &conf, double alpha);
template void smearing_APE_2d(std::vector<abelian> &conf1,
                              std::vector<abelian> &conf2, int mu, int nu,
                              double alpha);
template void smearing_HYP(std::vector<abelian> &array, double alpha1,
                           double alpha2, double alpha3);

// su3
template std::vector<std::vector<su3>>
separate_smearing(std::vector<su3> &conf);
template void smearing_APE(std::vector<su3> &conf, double alpha);
template void smearing_APE_2d(std::vector<su3> &conf1, std::vector<su3> &conf2,
                              int mu, int nu, double alpha);
template void smearing_HYP(std::vector<su3> &array, double alpha1,
                           double alpha2, double alpha3);

// su3_abelian
template std::vector<std::vector<su3_abelian>>
separate_smearing(std::vector<su3_abelian> &conf);
template void smearing_APE(std::vector<su3_abelian> &conf, double alpha);
template void smearing_APE_2d(std::vector<su3_abelian> &conf1,
                              std::vector<su3_abelian> &conf2, int mu, int nu,
                              double alpha);
template void smearing_HYP(std::vector<su3_abelian> &array, double alpha1,
                           double alpha2, double alpha3);

// su3_angles
template std::vector<std::vector<su3_angles>>
separate_smearing(std::vector<su3_angles> &conf);
template void smearing_APE_2d(std::vector<su3_angles> &conf1,
                              std::vector<su3_angles> &conf2, int mu, int nu,
                              double alpha);
template void smearing_HYP(std::vector<su3_angles> &array, double alpha1,
                           double alpha2, double alpha3);