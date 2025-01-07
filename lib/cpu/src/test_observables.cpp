#include "../include/indexing.h"
#include "../include/matrix.h"
#include <Eigen/Dense>

#include <map>
#include <omp.h>
#include <vector>

extern int x_size;
extern int y_size;
extern int z_size;
extern int t_size;

template <class T>
double plaket_site(const std::vector<T> &conf, std::vector<int> &lat_dim,
                   std::vector<int> &lat_coord) {
  su3 A;
  double plaket = 0;
  std::vector<int> lat_coord1 = lat_coord;
  for (int mu = 0; mu < 3; mu++) {
    for (int nu = mu + 1; nu < 4; nu++) {
      A = conf[get_index_matrix(lat_coord1, mu)];
      lat_coord1[mu] = (lat_coord1[mu] + 1) % lat_dim[mu];
      A = A * conf[get_index_matrix(lat_coord1, nu)];
      lat_coord1[mu] = (lat_coord1[mu] + lat_dim[mu] - 1) % lat_dim[mu];
      lat_coord1[nu] = (lat_coord1[nu] + 1) % lat_dim[nu];
      A = A ^ conf[get_index_matrix(lat_coord1, mu)];
      lat_coord1[nu] = (lat_coord1[nu] + lat_dim[nu] - 1) % lat_dim[nu];
      plaket += A.multiply_conj_tr(conf[get_index_matrix(lat_coord1, nu)]);
    }
  }
  return plaket / 6;
}

template <class T> double plaket_indexed(const std::vector<su3> &conf) {
  double plaket;
  std::vector<int> lat_coord(4);
  std::vector<int> lat_dim = {x_size, y_size, z_size, t_size};
#pragma omp parallel for collapse(4) private(lat_coord) firstprivate(lat_dim)  \
    reduction(+ : plaket)
  for (int t = 0; t < lat_dim[3]; t++) {
    for (int z = 0; z < lat_dim[2]; z++) {
      for (int y = 0; y < lat_dim[1]; y++) {
        for (int x = 0; x < lat_dim[0]; x++) {
          lat_coord = {x, y, z, t};
          plaket += plaket_site(conf, lat_dim, lat_coord);
        }
      }
    }
  }
  return plaket / (lat_dim[0] * lat_dim[1] * lat_dim[2] * lat_dim[3]);
}

template <class T>
void wilson_lines_prolong(const std::vector<T> &conf,
                          std::vector<T> &wilson_lines, int length, int mu) {
  std::vector<int> lat_coord(4);
  std::vector<int> lat_dim = {x_size, y_size, z_size, t_size};
  int index;
#pragma omp parallel for collapse(4) private(lat_coord, index)                 \
    firstprivate(lat_dim, mu)
  for (int t = 0; t < lat_dim[3]; t++) {
    for (int z = 0; z < lat_dim[2]; z++) {
      for (int y = 0; y < lat_dim[1]; y++) {
        for (int x = 0; x < lat_dim[0]; x++) {
          lat_coord = {x, y, z, t};
          index = get_index_site(lat_coord);
          lat_coord[mu] = (lat_coord[mu] + length) % lat_dim[mu];
          wilson_lines[index] =
              wilson_lines[index] * conf[get_index_matrix(lat_coord, mu)];
        }
      }
    }
  }
}

template <>
void wilson_lines_prolong(const std::vector<Eigen::Matrix3cd> &conf,
                          std::vector<Eigen::Matrix3cd> &wilson_lines,
                          int length, int mu) {
  std::vector<int> lat_coord(4);
  std::vector<int> lat_dim = {x_size, y_size, z_size, t_size};
  int index;
#pragma omp parallel for collapse(4) private(lat_coord, index)                 \
    firstprivate(lat_dim, mu)
  for (int t = 0; t < lat_dim[3]; t++) {
    for (int z = 0; z < lat_dim[2]; z++) {
      for (int y = 0; y < lat_dim[1]; y++) {
        for (int x = 0; x < lat_dim[0]; x++) {
          lat_coord = {x, y, z, t};
          index = get_index_site(lat_coord);
          lat_coord[mu] = (lat_coord[mu] + length) % lat_dim[mu];
          wilson_lines[index] *= conf[get_index_matrix(lat_coord, mu)];
        }
      }
    }
  }
}

template <class T>
std::vector<T> wilson_lines_get_length_one(const std::vector<T> &conf, int mu) {
  std::vector<T> wilson_lines(conf.size() / 4);
  std::vector<int> lat_coord(4);
  std::vector<int> lat_dim = {x_size, y_size, z_size, t_size};
#pragma omp parallel for collapse(4) private(lat_coord)                        \
    firstprivate(mu, lat_dim)
  for (int t = 0; t < lat_dim[3]; t++) {
    for (int z = 0; z < lat_dim[2]; z++) {
      for (int y = 0; y < lat_dim[1]; y++) {
        for (int x = 0; x < lat_dim[0]; x++) {
          lat_coord = {x, y, z, t};
          wilson_lines[get_index_site(lat_coord)] =
              conf[get_index_matrix(lat_coord, mu)];
        }
      }
    }
  }
  return wilson_lines;
}

template <>
std::vector<Eigen::Matrix3cd>
wilson_lines_get_length_one(const std::vector<Eigen::Matrix3cd> &conf, int mu) {
  std::vector<Eigen::Matrix3cd> wilson_lines(conf.size() / 4);
  std::vector<int> lat_coord(4);
  std::vector<int> lat_dim = {x_size, y_size, z_size, t_size};
#pragma omp parallel for collapse(4) private(lat_coord)                        \
    firstprivate(mu, lat_dim)
  for (int t = 0; t < lat_dim[3]; t++) {
    for (int z = 0; z < lat_dim[2]; z++) {
      for (int y = 0; y < lat_dim[1]; y++) {
        for (int x = 0; x < lat_dim[0]; x++) {
          lat_coord = {x, y, z, t};
          wilson_lines[get_index_site(lat_coord)] =
              conf[get_index_matrix(lat_coord, mu)];
        }
      }
    }
  }
  return wilson_lines;
}

template <class T>
std::vector<T> wilson_lines_indexed(const std::vector<T> &conf, int length,
                                    int mu) {
  std::vector<T> wilson_lines = wilson_lines_get_length_one(conf, mu);
  for (int i = 1; i < length; i++) {
    wilson_lines_prolong(conf, wilson_lines, length, mu);
  }
  return wilson_lines;
}

std::vector<Eigen::Matrix3cd>
wilson_lines_indexed(const std::vector<Eigen::Matrix3cd> &conf, int length,
                     int mu) {
  std::vector<Eigen::Matrix3cd> wilson_lines =
      wilson_lines_get_length_one(conf, mu);
  for (int i = 1; i < length; i++) {
    wilson_lines_prolong(conf, wilson_lines, length, mu);
  }
  return wilson_lines;
}

template <class T>
double wilson_plane_indexed(const std::vector<T> &wilson_lines_mu,
                            const std::vector<T> &wilson_lines_nu, int mu,
                            int nu, int length_mu, int length_nu) {
  std::vector<int> lat_dim = {x_size, y_size, z_size, t_size};
  T wilson_loop;
  double result = 0;
  std::vector<int> lat_coord(4);
  int index;
#pragma omp parallel for collapse(4) private(lat_coord, wilson_loop, index)    \
    firstprivate(lat_dim, mu, nu, length_mu, length_nu) reduction(+ : result)
  for (int t = 0; t < lat_dim[3]; t++) {
    for (int z = 0; z < lat_dim[2]; z++) {
      for (int y = 0; y < lat_dim[1]; y++) {
        for (int x = 0; x < lat_dim[0]; x++) {
          lat_coord = {x, y, z, t};
          index = get_index_site(lat_coord);
          lat_coord[mu] = (lat_coord[mu] + length_mu) % lat_dim[mu];
          wilson_loop = wilson_lines_mu[index] *
                        wilson_lines_nu[get_index_site(lat_coord)];
          lat_coord[nu] = (lat_coord[nu] + length_nu) % lat_dim[nu];
          lat_coord[mu] =
              (lat_coord[mu] + lat_dim[mu] - length_mu) % lat_dim[mu];
          wilson_loop =
              wilson_loop ^ wilson_lines_mu[get_index_site(lat_coord)];
          lat_coord[nu] =
              (lat_coord[nu] + lat_dim[nu] - length_nu) % lat_dim[nu];
          result += wilson_loop.multiply_conj_tr(
              wilson_lines_nu[get_index_site(lat_coord)]);
        }
      }
    }
  }
  return result / (x_size * y_size * z_size * t_size);
}

template <class T>
double wilson_plane_indexed_single_rxt(
    const std::vector<T> &wilson_lines_mu,
    const std::vector<std::vector<T>> &wilson_lines_nu, int mu, int length_mu,
    int length_nu) {
  std::vector<int> lat_dim = {x_size, y_size, z_size, t_size};
  T wilson_loop;
  double result = 0;
  std::vector<int> lat_coord(4);
  int index1;
  int index2;
#pragma omp parallel for collapse(4) private(lat_coord, wilson_loop, index1,   \
                                                 index2)                       \
    firstprivate(lat_dim, mu, length_mu, length_nu) reduction(+ : result)
  for (int t = 0; t < lat_dim[3]; t++) {
    for (int z = 0; z < lat_dim[2]; z++) {
      for (int y = 0; y < lat_dim[1]; y++) {
        for (int x = 0; x < lat_dim[0]; x++) {
          lat_coord = {x, y, z, t};
          index1 = get_index_site(lat_coord);
          lat_coord[mu] = (lat_coord[mu] + length_mu) % lat_dim[mu];
          index2 = get_index_site(lat_coord);
          lat_coord[mu] =
              (lat_coord[mu] + lat_dim[mu] - length_mu) % lat_dim[mu];
          for (int nu = 0; nu < 3; nu++) {
            wilson_loop = wilson_lines_mu[index1] * wilson_lines_nu[nu][index2];
            lat_coord[nu] = (lat_coord[nu] + length_nu) % lat_dim[nu];
            wilson_loop =
                wilson_loop ^ wilson_lines_mu[get_index_site(lat_coord)];
            lat_coord[nu] =
                (lat_coord[nu] + lat_dim[nu] - length_nu) % lat_dim[nu];
            result += wilson_loop.multiply_conj_tr(
                wilson_lines_nu[nu][get_index_site(lat_coord)]);
          }
        }
      }
    }
  }
  return result / (x_size * y_size * z_size * t_size * 3);
}

template <>
double wilson_plane_indexed_single_rxt(
    const std::vector<su3> &wilson_lines_mu,
    const std::vector<std::vector<su3>> &wilson_lines_nu, int mu, int length_mu,
    int length_nu) {
  std::vector<int> lat_dim = {x_size, y_size, z_size, t_size};
  su3 wilson_loop;
  double result = 0;
  std::vector<int> lat_coord(4);
  int index1;
  int index2;
#pragma omp parallel for collapse(4) private(lat_coord, wilson_loop, index1,   \
                                                 index2)                       \
    firstprivate(lat_dim, mu, length_mu, length_nu) reduction(+ : result)
  for (int t = 0; t < lat_dim[3]; t++) {
    for (int z = 0; z < lat_dim[2]; z++) {
      for (int y = 0; y < lat_dim[1]; y++) {
        for (int x = 0; x < lat_dim[0]; x++) {
          lat_coord = {x, y, z, t};
          index1 = get_index_site(lat_coord);
          lat_coord[mu] = (lat_coord[mu] + length_mu) % lat_dim[mu];
          index2 = get_index_site(lat_coord);
          lat_coord[mu] =
              (lat_coord[mu] + lat_dim[mu] - length_mu) % lat_dim[mu];
          for (int nu = 0; nu < 3; nu++) {
            wilson_loop = wilson_lines_mu[index1] * wilson_lines_nu[nu][index2];
            lat_coord[nu] = (lat_coord[nu] + length_nu) % lat_dim[nu];
            wilson_loop =
                wilson_loop ^ wilson_lines_mu[get_index_site(lat_coord)];
            lat_coord[nu] =
                (lat_coord[nu] + lat_dim[nu] - length_nu) % lat_dim[nu];
            result += wilson_loop.multiply_conj_tr(
                wilson_lines_nu[nu][get_index_site(lat_coord)]);
          }
        }
      }
    }
  }
  return result / (x_size * y_size * z_size * t_size * 3);
}

double wilson_plane_indexed_single_rxt(
    const std::vector<Eigen::Matrix3cd> &wilson_lines_mu,
    const std::vector<std::vector<Eigen::Matrix3cd>> &wilson_lines_nu, int mu,
    int length_mu, int length_nu) {
  std::vector<int> lat_dim = {x_size, y_size, z_size, t_size};
  Eigen::Matrix3cd wilson_loop;
  double result = 0;
  std::vector<int> lat_coord(4);
  int index1;
  int index2;
#pragma omp parallel for collapse(4) private(lat_coord, wilson_loop, index1,   \
                                                 index2)                       \
    firstprivate(lat_dim, mu, length_mu, length_nu) reduction(+ : result)
  for (int t = 0; t < lat_dim[3]; t++) {
    for (int z = 0; z < lat_dim[2]; z++) {
      for (int y = 0; y < lat_dim[1]; y++) {
        for (int x = 0; x < lat_dim[0]; x++) {
          lat_coord = {x, y, z, t};
          index1 = get_index_site(lat_coord);
          lat_coord[mu] = (lat_coord[mu] + length_mu) % lat_dim[mu];
          index2 = get_index_site(lat_coord);
          lat_coord[mu] =
              (lat_coord[mu] + lat_dim[mu] - length_mu) % lat_dim[mu];
          for (int nu = 0; nu < 3; nu++) {
            wilson_loop = wilson_lines_mu[index1] * wilson_lines_nu[nu][index2];
            lat_coord[nu] = (lat_coord[nu] + length_nu) % lat_dim[nu];
            wilson_loop =
                wilson_loop *
                wilson_lines_mu[get_index_site(lat_coord)].adjoint().eval();
            lat_coord[nu] =
                (lat_coord[nu] + lat_dim[nu] - length_nu) % lat_dim[nu];
            result += (wilson_loop *
                       wilson_lines_nu[nu][get_index_site(lat_coord)].adjoint())
                          .trace()
                          .real() /
                      3;
          }
        }
      }
    }
  }
  return result / (x_size * y_size * z_size * t_size * 3);
}

template <class T>
std::map<std::tuple<int, int>, double>
wilson_parallel_indexed(const std::vector<T> &conf, int r_min, int r_max,
                        int time_min, int time_max) {
  std::map<std::tuple<int, int>, double> wilson_loops;
  std::vector<std::vector<T>> time_lines(time_max - time_min + 1);
  std::vector<T> space_lines;
  time_lines[0] = wilson_lines_indexed(conf, time_min, 3);
  for (int t = time_min + 1; t <= time_max; t++) {
    time_lines[t - time_min] = time_lines[t - time_min - 1];
    wilson_lines_prolong(conf, time_lines[t - time_min], t - 1, 3);
  }
  for (int mu = 0; mu < 3; mu++) {
    space_lines = wilson_lines_indexed(conf, r_min, mu);
    for (int r = r_min; r <= r_max; r++) {
      for (int t = time_min; t <= time_max; t++) {
        wilson_loops[std::tuple<int, int>(t, r)] += wilson_plane_indexed(
            space_lines, time_lines[t - time_min], mu, 3, r, t);
      }
      wilson_lines_prolong(conf, space_lines, r, mu);
    }
  }
  for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
    it->second = it->second / 3;
  }
  return wilson_loops;
}

template <class T>
std::map<std::tuple<int, int>, double>
wilson_parallel_indexed_single_rxt(const std::vector<T> &conf, int r_min,
                                   int r_max, int time_min, int time_max) {
  double start_time1;
  double end_time1;
  double search_time1 = 0;
  double start_time2;
  double end_time2;
  double search_time2 = 0;
  std::map<std::tuple<int, int>, double> wilson_loops;
  std::vector<std::vector<T>> time_lines(time_max - time_min + 1);
  std::vector<std::vector<T>> space_lines(3);
  start_time1 = omp_get_wtime();
  time_lines[0] = wilson_lines_indexed(conf, time_min, 3);
  for (int t = time_min + 1; t <= time_max; t++) {
    time_lines[t - time_min] = time_lines[t - time_min - 1];
    wilson_lines_prolong(conf, time_lines[t - time_min], t - 1, 3);
  }
  for (int mu = 0; mu < 3; mu++) {
    space_lines[mu] = wilson_lines_indexed(conf, r_min, mu);
  }
  end_time1 = omp_get_wtime();
  search_time1 += end_time1 - start_time1;
  for (int r = r_min; r <= r_max; r++) {
    for (int t = time_min; t <= time_max; t++) {
      start_time2 = omp_get_wtime();
      wilson_loops[std::tuple<int, int>(t, r)] +=
          wilson_plane_indexed_single_rxt(time_lines[t - time_min], space_lines,
                                          3, t, r);
      end_time2 = omp_get_wtime();
      search_time2 += end_time2 - start_time2;
    }
    start_time1 = omp_get_wtime();
    for (int mu = 0; mu < 3; mu++) {
      wilson_lines_prolong(conf, space_lines[mu], r, mu);
    }
    end_time1 = omp_get_wtime();
    search_time1 += end_time1 - start_time1;
  }
  std::cout << "wilson_lines time: " << search_time1 << std::endl;
  std::cout << "wilson_plane time: " << search_time2 << std::endl;
  return wilson_loops;
}

template <>
std::map<std::tuple<int, int>, double>
wilson_parallel_indexed_single_rxt(const std::vector<Eigen::Matrix3cd> &conf,
                                   int r_min, int r_max, int time_min,
                                   int time_max) {
  double start_time1;
  double end_time1;
  double search_time1 = 0;
  double start_time2;
  double end_time2;
  double search_time2 = 0;
  std::map<std::tuple<int, int>, double> wilson_loops;
  std::vector<std::vector<Eigen::Matrix3cd>> time_lines(time_max - time_min +
                                                        1);
  std::vector<std::vector<Eigen::Matrix3cd>> space_lines(3);
  start_time1 = omp_get_wtime();
  time_lines[0] = wilson_lines_indexed(conf, time_min, 3);
  for (int t = time_min + 1; t <= time_max; t++) {
    time_lines[t - time_min] = time_lines[t - time_min - 1];
    wilson_lines_prolong(conf, time_lines[t - time_min], t - 1, 3);
  }
  for (int mu = 0; mu < 3; mu++) {
    space_lines[mu] = wilson_lines_indexed(conf, r_min, mu);
  }
  end_time1 = omp_get_wtime();
  search_time1 += end_time1 - start_time1;
  for (int r = r_min; r <= r_max; r++) {
    for (int t = time_min; t <= time_max; t++) {
      start_time2 = omp_get_wtime();
      wilson_loops[std::tuple<int, int>(t, r)] +=
          wilson_plane_indexed_single_rxt(time_lines[t - time_min], space_lines,
                                          3, t, r);
      end_time2 = omp_get_wtime();
      search_time2 += end_time2 - start_time2;
    }
    start_time1 = omp_get_wtime();
    for (int mu = 0; mu < 3; mu++) {
      wilson_lines_prolong(conf, space_lines[mu], r, mu);
    }
    end_time1 = omp_get_wtime();
    search_time1 += end_time1 - start_time1;
  }
  std::cout << "wilson_lines time: " << search_time1 << std::endl;
  std::cout << "wilson_plane time: " << search_time2 << std::endl;
  return wilson_loops;
}

template std::map<std::tuple<int, int>, double>
wilson_parallel_indexed(const std::vector<su3> &conf, int r_min, int r_max,
                        int time_min, int time_max);
template std::map<std::tuple<int, int>, double>
wilson_parallel_indexed_single_rxt(const std::vector<su3> &conf, int r_min,
                                   int r_max, int time_min, int time_max);