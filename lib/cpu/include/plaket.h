#pragma once

#include "data.h"
#include "indexing.h"

#include <vector>

template <class DataPattern, class MatrixType>
double plaket_plane(const Data::LatticeData<DataPattern, MatrixType> &conf,
                    DataPattern &data_pattern, int mu, int nu) {
  MatrixType A;
  double plaket = 0;
  A = conf[data_pattern.get_index_link(mu)];
  data_pattern.move_forward(1, mu);
  A = A * conf[data_pattern.get_index_link(nu)];
  data_pattern.move_backward(1, mu);
  data_pattern.move_forward(1, nu);
  A = A ^ conf[data_pattern.get_index_link(mu)];
  data_pattern.move_backward(1, nu);
  plaket += A.multiply_conj_tr(conf[data_pattern.get_index_link(nu)]);
  return plaket;
}

template <class DataPattern, class MatrixType>
double plaket_site(const Data::LatticeData<DataPattern, MatrixType> &conf,
                   DataPattern &data_pattern) {
  MatrixType A;
  double plaket = 0;
  for (int mu = 0; mu < 3; mu++) {
    for (int nu = mu + 1; nu < 4; nu++) {
      A = conf[data_pattern.get_index_link(mu)];
      data_pattern.move_forward(1, mu);
      A = A * conf[data_pattern.get_index_link(nu)];
      data_pattern.move_backward(1, mu);
      data_pattern.move_forward(1, nu);
      A = A ^ conf[data_pattern.get_index_link(mu)];
      data_pattern.move_backward(1, nu);
      plaket += A.multiply_conj_tr(conf[data_pattern.get_index_link(nu)]);
    }
  }
  return plaket / 6;
}

template <class DataPattern, class MatrixType>
double plaket_site_time(const Data::LatticeData<DataPattern, MatrixType> &conf,
                        DataPattern &data_pattern) {
  MatrixType A;
  double plaket = 0;
  int nu = 3;
  for (int mu = 0; mu < 3; mu++) {
    A = conf[data_pattern.get_index_link(mu)];
    data_pattern.move_forward(1, mu);
    A = A * conf[data_pattern.get_index_link(nu)];
    data_pattern.move_backward(1, mu);
    data_pattern.move_forward(1, nu);
    A = A ^ conf[data_pattern.get_index_link(mu)];
    data_pattern.move_backward(1, nu);
    plaket += A.multiply_conj_tr(conf[data_pattern.get_index_link(nu)]);
  }
  return plaket / 3;
}

template <class DataPattern, class MatrixType>
double plaket_site_space(const Data::LatticeData<DataPattern, MatrixType> &conf,
                         DataPattern &data_pattern) {
  MatrixType A;
  double plaket = 0;
  for (int mu = 0; mu < 2; mu++) {
    for (int nu = mu + 1; nu < 3; nu++) {
      A = conf[data_pattern.get_index_link(mu)];
      data_pattern.move_forward(1, mu);
      A = A * conf[data_pattern.get_index_link(nu)];
      data_pattern.move_backward(1, mu);
      data_pattern.move_forward(1, nu);
      A = A ^ conf[data_pattern.get_index_link(mu)];
      data_pattern.move_backward(1, nu);
      plaket += A.multiply_conj_tr(conf[data_pattern.get_index_link(nu)]);
    }
  }
  return plaket / 3;
}

template <class DataPattern, class MatrixType>
double plaket(const Data::LatticeData<DataPattern, MatrixType> &conf) {
  double plaket;
  DataPattern data_pattern(conf.lat_dim);
#pragma omp parallel for collapse(4) firstprivate(data_pattern)                \
    reduction(+ : plaket)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          plaket += plaket_site(conf, data_pattern);
        }
      }
    }
  }
  return plaket / (data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                   data_pattern.lat_dim[2] * data_pattern.lat_dim[3]);
}

template <class DataPattern, class MatrixType>
double plaket_time(const Data::LatticeData<DataPattern, MatrixType> &conf) {
  double plaket;
  DataPattern data_pattern(conf.lat_dim);
#pragma omp parallel for collapse(4) firstprivate(data_pattern)                \
    reduction(+ : plaket)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          plaket += plaket_site_time(conf, data_pattern);
        }
      }
    }
  }
  return plaket / (data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                   data_pattern.lat_dim[2] * data_pattern.lat_dim[3]);
}

template <class DataPattern, class MatrixType>
double plaket_space(const Data::LatticeData<DataPattern, MatrixType> &conf) {
  double plaket;
  DataPattern data_pattern(conf.lat_dim);
#pragma omp parallel for collapse(4) firstprivate(data_pattern)                \
    reduction(+ : plaket)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          plaket += plaket_site_space(conf, data_pattern);
        }
      }
    }
  }
  return plaket / (data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                   data_pattern.lat_dim[2] * data_pattern.lat_dim[3]);
}

// preserves data_pattern coordinates
// counterclockwise direction in (mu,nu) plane
template <class DataPattern, class MatrixType>
inline MatrixType
plaket_left_down(const Data::LatticeData<DataPattern, MatrixType> &conf,
                 DataPattern &data_pattern, int mu, int nu) {
  MatrixType A = conf[data_pattern.get_index_link(mu)];
  data_pattern.move_forward(1, mu);
  A = A * conf[data_pattern.get_index_link(nu)];
  data_pattern.move_backward(1, mu);
  data_pattern.move_forward(1, nu);
  A = A ^ conf[data_pattern.get_index_link(mu)];
  data_pattern.move_backward(1, nu);
  A = A ^ conf[data_pattern.get_index_link(nu)];
  return A;
}

// preserves data_pattern coordinates
template <class DataPattern, class MatrixType>
inline MatrixType
plaket_left_up(const Data::LatticeData<DataPattern, MatrixType> &conf,
               DataPattern &data_pattern, int mu, int nu) {
  data_pattern.move_backward(1, nu);
  MatrixType A = conf[data_pattern.get_index_link(nu)].conj();
  A = A * conf[data_pattern.get_index_link(mu)];
  data_pattern.move_forward(1, mu);
  A = A * conf[data_pattern.get_index_link(nu)];
  data_pattern.move_backward(1, mu);
  data_pattern.move_forward(1, nu);
  A = A ^ conf[data_pattern.get_index_link(mu)];
  return A;
}

// preserves data_pattern coordinates
template <class DataPattern, class MatrixType>
inline MatrixType
plaket_right_down(const Data::LatticeData<DataPattern, MatrixType> &conf,
                  DataPattern &data_pattern, int mu, int nu) {
  MatrixType A = conf[data_pattern.get_index_link(nu)];
  data_pattern.move_forward(1, nu);
  data_pattern.move_backward(1, mu);
  A = A ^ conf[data_pattern.get_index_link(mu)];
  data_pattern.move_backward(1, nu);
  A = A ^ conf[data_pattern.get_index_link(nu)];
  A = A * conf[data_pattern.get_index_link(mu)];
  data_pattern.move_forward(1, mu);
  return A;
}

// preserves data_pattern coordinates
template <class DataPattern, class MatrixType>
inline MatrixType
plaket_right_up(const Data::LatticeData<DataPattern, MatrixType> &conf,
                DataPattern &data_pattern, int mu, int nu) {
  data_pattern.move_backward(1, mu);
  MatrixType A = conf[data_pattern.get_index_link(mu)].conj();
  data_pattern.move_backward(1, nu);
  A = A ^ conf[data_pattern.get_index_link(nu)];
  A = A * conf[data_pattern.get_index_link(mu)];
  data_pattern.move_forward(1, mu);
  A = A * conf[data_pattern.get_index_link(nu)];
  data_pattern.move_forward(1, nu);
  return A;
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType> calculate_plaket_time_left_down(
    const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> vec(data_pattern.get_lattice_size() * 3);
  int index;
#pragma omp parallel for collapse(4) firstprivate(data_pattern) private(index)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_index_site() * 3;
          for (int dir = 0; dir < 3; dir++) {
            vec[index + dir] = plaket_left_down(conf, data_pattern, dir, 3);
          }
        }
      }
    }
  }
  return vec;
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType> calculate_plaket_time_left_up(
    const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> vec(data_pattern.get_lattice_size() * 3);
  int index;
#pragma omp parallel for collapse(4) firstprivate(data_pattern) private(index)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_index_site() * 3;
          for (int dir = 0; dir < 3; dir++) {
            vec[index + dir] = plaket_left_up(conf, data_pattern, dir, 3);
          }
        }
      }
    }
  }
  return vec;
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType> calculate_plaket_time_right_down(
    const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> vec(data_pattern.get_lattice_size() * 3);
  int index;
#pragma omp parallel for collapse(4) firstprivate(data_pattern) private(index)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_index_site() * 3;
          for (int dir = 0; dir < 3; dir++) {
            vec[index + dir] = plaket_right_down(conf, data_pattern, dir, 3);
          }
        }
      }
    }
  }
  return vec;
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType> calculate_plaket_time_right_up(
    const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> vec(data_pattern.get_lattice_size() * 3);
  int index;
#pragma omp parallel for collapse(4) firstprivate(data_pattern) private(index)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_index_site() * 3;
          for (int dir = 0; dir < 3; dir++) {
            vec[index + dir] = plaket_right_up(conf, data_pattern, dir, 3);
          }
        }
      }
    }
  }
  return vec;
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType> calculate_plaket_schwinger_time_left(
    const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> vec(data_pattern.get_lattice_size() * 3);
  int index;
#pragma omp parallel for collapse(4) firstprivate(data_pattern) private(index)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_index_site() * 3;
          for (int dir = 0; dir < 3; dir++) {
            vec[index + dir] = (plaket_left_down(conf, data_pattern, dir, 3) +
                                plaket_left_up(conf, data_pattern, dir, 3)) *
                               0.5;
          }
        }
      }
    }
  }
  return vec;
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType> calculate_plaket_schwinger_time_right(
    const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> vec(data_pattern.get_lattice_size() * 3);
  int index;
#pragma omp parallel for collapse(4) firstprivate(data_pattern) private(index)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_index_site() * 3;
          for (int dir = 0; dir < 3; dir++) {
            vec[index + dir] = (plaket_right_down(conf, data_pattern, dir, 3) +
                                plaket_right_up(conf, data_pattern, dir, 3)) *
                               0.5;
          }
        }
      }
    }
  }
  return vec;
}

template <class DataPattern, class MatrixType>
std::vector<double> calculate_plaket_schwinger_time_left_tr(
    const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<double> vec(data_pattern.get_lattice_size() * 3);
  int index;
#pragma omp parallel for collapse(4) firstprivate(data_pattern) private(index)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_index_site() * 3;
          for (int dir = 0; dir < 3; dir++) {
            vec[index + dir] =
                (plaket_left_down(conf, data_pattern, dir, 3).tr() +
                 plaket_left_up(conf, data_pattern, dir, 3).tr()) *
                0.5;
          }
        }
      }
    }
  }
  return vec;
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType> calculate_plaket_schwinger_time_tr(
    const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> vec(data_pattern.get_lattice_size() * 3);
  int index;
#pragma omp parallel for collapse(4) firstprivate(data_pattern) private(index)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          index = data_pattern.get_lattice_size() * 3;
          for (int dir = 0; dir < 3; dir++) {
            for (int mu = 0; mu < 3; mu++) {
              if (dir != mu) {
                vec[index + dir] =
                    vec[index + dir] +
                    (plaket_right_down(conf, data_pattern, mu, 3).conj() +
                     plaket_right_up(conf, data_pattern, mu, 3).conj() +
                     plaket_left_down(conf, data_pattern, mu, 3) +
                     plaket_left_up(conf, data_pattern, mu, 3)) *
                        0.125;
              }
            }
          }
        }
      }
    }
  }
  return vec;
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType> calculate_plaket_schwinger_space_l(
    const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> vec(data_pattern.get_lattice_size() * 3);
  int index;
#pragma omp parallel for collapse(4) firstprivate(data_pattern) private(index)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          index = data_pattern.get_lattice_size() * 3;
          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
              if (i != j) {
                for (int k = 0; k < 3; k++) {
                  if (k != i && k != j) {
                    vec[index + i] =
                        vec[index + i] +
                        (plaket_right_down(conf, data_pattern, j, k) +
                         plaket_right_up(conf, data_pattern, j, k) +
                         plaket_left_down(conf, data_pattern, j, k) +
                         plaket_left_up(conf, data_pattern, j, k)) *
                            0.25;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return vec;
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType> calculate_plaket_schwinger_space_tr(
    const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> vec(data_pattern.get_lattice_size() * 3);
  int index;
#pragma omp parallel for collapse(4) firstprivate(data_pattern) private(index)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          index = data_pattern.get_lattice_size() * 3;
          for (int dir = 0; dir < 3; dir++) {
            for (int i = 0; i < 3; i++) {
              if (i != dir) {
                for (int j = 0; j < 3; j++) {
                  if (j != i) {
                    for (int k = 0; k < 3; k++) {
                      if (k != i && k != j) {
                        vec[index + dir] =
                            vec[index + dir] +
                            (plaket_right_down(conf, data_pattern, j, k) +
                             plaket_right_up(conf, data_pattern, j, k) +
                             plaket_left_down(conf, data_pattern, j, k) +
                             plaket_left_up(conf, data_pattern, j, k)) *
                                0.125;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return vec;
}

template <class DataPattern, class MatrixType>
std::vector<double> calculate_plaket_time_tr(
    const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<double> vec(data_pattern.get_lattice_size() * 3);
  int index;
#pragma omp parallel for collapse(4) firstprivate(data_pattern) private(index)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_index_site() * 3;
          for (int mu = 0; mu < 3; mu++) {
            vec[index + mu] =
                plaket_right_up(conf, data_pattern, mu, 3).tr_real();
          }
        }
      }
    }
  }
  return vec;
}

template <class DataPattern, class MatrixType>
std::vector<double> calculate_plaket_space_tr(
    const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<double> vec(data_pattern.get_lattice_size() * 3);
  int index;
  int count;
#pragma omp parallel for collapse(4)                                           \
    firstprivate(data_pattern) private(index, count)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_index_site() * 3;
          count = 0;
          for (int mu = 0; mu < 2; mu++) {
            for (int nu = mu + 1; nu < 3; nu++) {
              vec[index + count] =
                  plaket_left_down(conf, data_pattern, mu, nu).tr_real();
              count++;
            }
          }
        }
      }
    }
  }
  return vec;
}

template <class DataPattern, class MatrixType>
std::vector<double>
plaket_aver_time_tr(const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<double> plaket_tr = calculate_plaket_time_tr(conf);
  std::vector<double> vec(data_pattern.get_lattice_size() * 3);
  int index;
  double trace_aver;
#pragma omp parallel for collapse(4)                                           \
    firstprivate(data_pattern) private(index, trace_aver)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_index_site() * 3;
          trace_aver = 0;
          for (int mu = 0; mu < 3; mu++) {
            trace_aver += plaket_tr[index + mu];
            data_pattern.move_backward(1, mu);
            trace_aver += plaket_tr[data_pattern.get_index_site() * 3 + mu];
            data_pattern.move_backward(1, 3);
            trace_aver += plaket_tr[data_pattern.get_index_site() * 3 + mu];
            data_pattern.move_forward(1, mu);
            trace_aver += plaket_tr[data_pattern.get_index_site() * 3 + mu];
            vec[index + mu] = trace_aver / 4;
            data_pattern.move_forward(1, 3);
          }
        }
      }
    }
  }
  return vec;
}

template <class DataPattern, class MatrixType>
std::vector<double>
plaket_aver_space_tr(const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<double> plaket_tr = calculate_plaket_space_tr(conf);
  std::vector<double> vec(data_pattern.get_lattice_size() * 3);
  int index;
  int count;
  double trace_aver;
#pragma omp parallel for collapse(4)                                           \
    firstprivate(data_pattern) private(index, count, trace_aver)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_index_site() * 3;
          count = 0;
          trace_aver = 0;
          for (int mu = 0; mu < 2; mu++) {
            for (int nu = mu + 1; nu < 3; nu++) {
              trace_aver += plaket_tr[index + count];
              data_pattern.move_backward(1, mu);
              trace_aver +=
                  plaket_tr[data_pattern.get_index_site() * 3 + count];
              data_pattern.move_backward(1, 3);
              trace_aver +=
                  plaket_tr[data_pattern.get_index_site() * 3 + count];
              data_pattern.move_forward(1, mu);
              trace_aver +=
                  plaket_tr[data_pattern.get_index_site() * 3 + count];
              vec[index + count] = trace_aver / 4;
              data_pattern.move_forward(1, 3);
              count++;
            }
          }
        }
      }
    }
  }
  return vec;
}

// average over all directions for trace of temporal wilson loop
// is used for flux tube
template <class DataPattern, class MatrixType>
std::vector<double> plaket_time_site_average_tr(
    const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<double> vec(data_pattern.get_lattice_size());
  double trace_aver;
#pragma omp parallel for collapse(4)                                           \
    firstprivate(data_pattern) private(trace_aver)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          trace_aver = 0;
          for (int mu = 0; mu < 3; mu++) {
            trace_aver += plaket_plane(conf, data_pattern, mu, 3);
            data_pattern.move_backward(1, mu);
            trace_aver += plaket_plane(conf, data_pattern, mu, 3);
            data_pattern.move_backward(1, 3);
            trace_aver += plaket_plane(conf, data_pattern, mu, 3);
            data_pattern.move_forward(1, mu);
            trace_aver += plaket_plane(conf, data_pattern, mu, 3);
            data_pattern.move_forward(1, 3);
          }
          vec[data_pattern.get_index_site()] = trace_aver / 12;
        }
      }
    }
  }
  return vec;
}

// average over all directions for trace of spacial wilson loop
// is used for flux tube
template <class DataPattern, class MatrixType>
std::vector<double> plaket_space_site_average_tr(
    const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<double> vec(data_pattern.get_lattice_size());
  double trace_aver;
#pragma omp parallel for collapse(4)                                           \
    firstprivate(data_pattern) private(trace_aver)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          trace_aver = 0;
          for (int mu = 0; mu < 2; mu++) {
            for (int nu = mu + 1; nu < 3; nu++) {
              trace_aver += plaket_plane(conf, data_pattern, mu, nu);
              data_pattern.move_backward(1, mu);
              trace_aver += plaket_plane(conf, data_pattern, mu, nu);
              data_pattern.move_backward(1, nu);
              trace_aver += plaket_plane(conf, data_pattern, mu, nu);
              data_pattern.move_forward(1, mu);
              trace_aver += plaket_plane(conf, data_pattern, mu, nu);
              data_pattern.move_forward(1, nu);
            }
          }
          vec[data_pattern.get_index_site()] = trace_aver / 12;
        }
      }
    }
  }
  return vec;
}