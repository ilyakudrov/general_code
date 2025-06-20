#pragma once

#include "data.h"
#include "indexing.h"
#include "link.h"
#include "matrix.h"
#include "wilson_loops.h"

#include <map>
#include <vector>

#pragma omp declare reduction(                                                 \
        vec_double_plus : std::vector<double> : std::transform(                \
                omp_out.begin(), omp_out.end(), omp_in.begin(),                \
                    omp_out.begin(), std::plus<double>()))                     \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

template <class DataPattern, class MatrixType>
std::vector<MatrixType> calculate_schwinger_lines_short(
    const Data::LatticeData<DataPattern, MatrixType> &conf, int d) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> vec(data_pattern.get_lattice_size);
#pragma omp parallel for collapse(4) firstprivate(data_pattern)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          for (int mu = 0; mu < 3; mu++) {
            vec[data_pattern.get_lattice_size * 3 + mu] =
                wilson_line(conf, data_pattern, mu, d);
          }
        }
      }
    }
  }
  return vec;
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType> calculate_wilson_loops_schwinger(
    const Data::LatticeData<DataPattern, MatrixType> &conf, int r, int time) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> vec(data_pattern.get_lattice_size * 3);
  int index;
#pragma omp parallel for collapse(4) firstprivate(data_pattern) private(index)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_lattice_size * 3;
          for (int mu = 0; mu < 3; mu++) {
            vec[index + mu] =
                wilson_loop_schwinger(conf, data_pattern, r, time, mu);
          }
        }
      }
    }
  }
  return vec;
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType> calculate_wilson_loops_schwinger_opposite(
    const Data::LatticeData<DataPattern, MatrixType> &conf, int r, int time) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> vec(data_pattern.get_lattice_size * 3);
  int index;
#pragma omp parallel for collapse(4) firstprivate(data_pattern) private(index)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_lattice_size * 3;
          for (int mu = 0; mu < 3; mu++) {
            vec[index + mu] =
                wilson_loop_schwinger_opposite(conf, data_pattern, r, time, mu);
          }
        }
      }
    }
  }
  return vec;
}

template <class DataPattern, class MatrixType>
std::map<int, double> wilson_plaket_schwinger_longitudinal_l(
    const std::vector<MatrixType> &wilson_loops,
    const std::vector<MatrixType> &wilson_loops_opposite,
    const std::vector<MatrixType> &plaket_left,
    const std::vector<MatrixType> &plaket_right,
    const std::vector<std::vector<MatrixType>> &schwinger_lines,
    DataPattern &data_pattern, int d_ouside, int time, int r) {
  std::vector<double> correlator(2 * d_ouside + r);
  MatrixType W;
  MatrixType A;
  MatrixType S;
  int d;
  int place;
#pragma omp parallel for collapse(4) private(W, A, S, d, place)                \
    firstprivate(d_ouside, data_pattern)                                       \
    reduction(vec_double_plus : correlator)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          for (int dir = 0; dir < 3; dir++) {
            data_pattern.lat_coord = {x, y, z, t};
            W = wilson_loops[data_pattern.get_lattice_size * 3 + dir];
            d = -d_ouside;
            data_pattern.move_backward(d_ouside - 1, dir);
            while (d < -1) {
              place = data_pattern.get_lattice_size * 3;
              S = schwinger_lines[abs(d) - 2][place + dir];
              A = W ^ S;
              A = A ^ plaket_right[place + dir];
              correlator[d + d_ouside] += A.multiply_tr(S);
              data_pattern.move_forward(1, dir);
              d++;
            }
            correlator[d + d_ouside] += W.multiply_conj_tr(
                plaket_right[data_pattern.get_lattice_size * 3 + dir]);
            d++;
            correlator[d + d_ouside] += W.multiply_tr(
                plaket_left[data_pattern.get_lattice_size * 3 + dir]);
            d++;
            while (d < (r + 1) / 2) {
              S = schwinger_lines[d - 1]
                                 [data_pattern.get_lattice_size * 3 + dir];
              A = W * S;
              data_pattern.move_forward(d, dir);
              A = A * plaket_left[data_pattern.get_lattice_size * 3 + dir];
              correlator[d + d_ouside] += A.multiply_conj_tr(S);
              data_pattern.move_backward(d, dir);
              d++;
            }
            data_pattern.move_forward(r, dir);
            W = wilson_loops_opposite[data_pattern.get_lattice_size * 3 + dir];
            data_pattern.move_backward(r / 2 - 1, dir);
            while (d < r - 1) {
              S = schwinger_lines[r - d - 2]
                                 [data_pattern.get_lattice_size * 3 + dir];
              A = W ^ S;
              A = A * plaket_right[data_pattern.get_lattice_size * 3 + dir];
              correlator[d + d_ouside] += A.multiply_tr(S);
              data_pattern.move_forward(1, dir);
              d++;
            }
            correlator[d + d_ouside] += W.multiply_tr(
                plaket_right[data_pattern.get_lattice_size * 3 + dir]);
            d++;
            correlator[d + d_ouside] += W.multiply_conj_tr(
                plaket_left[data_pattern.get_lattice_size * 3 + dir]);
            d++;
            while (d < r + d_ouside) {
              S = schwinger_lines[d - r - 1]
                                 [data_pattern.get_lattice_size * 3 + dir];
              A = W * S;
              data_pattern.move_forward(d - r, dir);
              A = A ^ plaket_left[data_pattern.get_lattice_size * 3 + dir];
              correlator[d + d_ouside] += A.multiply_conj_tr(S);
              data_pattern.move_backward(d - r, dir);
              d++;
            }
          }
        }
      }
    }
  }
  std::map<int, double> result;
  for (int i = 0; i < correlator.size(); i++) {
    result[i - d_ouside - r / 2] =
        correlator[i] / (data_pattern.get_lattice_size * 3);
  }
  return result;
}

template <class DataPattern, class MatrixType>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_longitudinal_l(
    const Data::LatticeData<DataPattern, MatrixType> &conf,
    const std::vector<std::vector<MatrixType>> &schwinger_lines_short,
    int T_min, int T_max, int R_min, int R_max, int d_ouside) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> plaket_time_left =
      calculate_plaket_schwinger_time_left(conf);
  std::vector<MatrixType> plaket_time_right =
      calculate_plaket_schwinger_time_right(conf);
  std::vector<MatrixType> wilson_loops;
  std::vector<MatrixType> wilson_loops_opposite;
  std::map<std::tuple<int, int, int>, double> result;
  std::map<int, double> schwinger_electric;
  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {
      wilson_loops = calculate_wilson_loops_schwinger(conf, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(conf, r, t);
      schwinger_electric = wilson_plaket_schwinger_longitudinal_l(
          wilson_loops, wilson_loops_opposite, plaket_time_left,
          plaket_time_right, schwinger_lines_short, data_pattern, d_ouside, t,
          r);
      for (auto it = schwinger_electric.begin(); it != schwinger_electric.end();
           ++it) {
        result[std::tuple<int, int, double>(t, r, it->first)] = it->second;
      }
    }
  }
  return result;
}

template <class DataPattern, class MatrixType>
std::map<int, double> wilson_plaket_schwinger_longitudinal_tr(
    const std::vector<MatrixType> &wilson_loops,
    const std::vector<MatrixType> &wilson_loops_opposite,
    const std::vector<MatrixType> &plaket,
    const std::vector<std::vector<MatrixType>> &schwinger_lines,
    DataPattern &data_pattern, int d_ouside, int time, int r) {
  std::vector<double> correlator(2 * d_ouside + r + 1);
  MatrixType W;
  MatrixType A;
  MatrixType S;
  int d;
  int place;
#pragma omp parallel for collapse(4) private(W, A, S, d, place)                \
    firstprivate(d_ouside, data_pattern)                                       \
    reduction(vec_double_plus : correlator)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          for (int dir = 0; dir < 3; dir++) {
            data_pattern.lat_coord = {x, y, z, t};
            W = wilson_loops[data_pattern.get_lattice_size * 3 + dir];
            d = -d_ouside;
            data_pattern.move_backward(d_ouside, dir);
            while (d < 0) {
              place = data_pattern.get_lattice_size * 3;
              S = schwinger_lines[abs(d) - 1][place + dir];
              A = W ^ S;
              A = A * plaket[place + dir];
              correlator[d + d_ouside] += A.multiply_tr(S);
              data_pattern.move_forward(1, dir);
              d++;
            }
            correlator[d + d_ouside] +=
                W.multiply_tr(plaket[data_pattern.get_lattice_size * 3 + dir]);
            d++;
            while (d <= r / 2) {
              S = schwinger_lines[d - 1]
                                 [data_pattern.get_lattice_size * 3 + dir];
              A = W * S;
              data_pattern.move_forward(d, dir);
              A = A * plaket[data_pattern.get_lattice_size * 3 + dir];
              correlator[d + d_ouside] += A.multiply_conj_tr(S);
              data_pattern.move_backward(d, dir);
              d++;
            }
            data_pattern.move_forward(r, dir);
            W = wilson_loops_opposite[data_pattern.get_lattice_size * 3 + dir];
            data_pattern.move_backward((r + 1) / 2 - 1, dir);
            while (d < r) {
              S = schwinger_lines[r - d - 1]
                                 [data_pattern.get_lattice_size * 3 + dir];
              A = W ^ S;
              A = A ^ plaket[data_pattern.get_lattice_size * 3 + dir];
              correlator[d + d_ouside] += A.multiply_tr(S);
              data_pattern.move_forward(1, dir);
              d++;
            }
            correlator[d + d_ouside] += W.multiply_conj_tr(
                plaket[data_pattern.get_lattice_size * 3 + dir]);
            d++;
            while (d <= r + d_ouside) {
              S = schwinger_lines[d - r - 1]
                                 [data_pattern.get_lattice_size * 3 + dir];
              A = W * S;
              data_pattern.move_forward(d - r, dir);
              A = A ^ plaket[data_pattern.get_lattice_size * 3 + dir];
              correlator[d + d_ouside] += A.multiply_conj_tr(S);
              data_pattern.move_backward(d - r, dir);
              d++;
            }
          }
        }
      }
    }
  }
  std::map<int, double> result;
  for (int i = 0; i < correlator.size(); i++) {
    result[i - d_ouside - r / 2] =
        correlator[i] / (data_pattern.get_lattice_size * 3);
  }
  return result;
}

template <class DataPattern, class MatrixType>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_longitudinal_tr(
    const Data::LatticeData<DataPattern, MatrixType> &conf,
    const std::vector<std::vector<MatrixType>> &schwinger_lines_short,
    int T_min, int T_max, int R_min, int R_max, int d_ouside) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> plaket_time =
      calculate_plaket_schwinger_time_tr(conf);
  std::vector<MatrixType> wilson_loops;
  std::vector<MatrixType> wilson_loops_opposite;
  std::map<std::tuple<int, int, int>, double> result;
  std::map<int, double> schwinger_electric;
  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {
      wilson_loops = calculate_wilson_loops_schwinger(conf, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(conf, r, t);
      schwinger_electric = wilson_plaket_schwinger_longitudinal_tr(
          wilson_loops, wilson_loops_opposite, plaket_time,
          schwinger_lines_short, data_pattern, d_ouside, t, r);
      for (auto it = schwinger_electric.begin(); it != schwinger_electric.end();
           ++it) {
        result[std::tuple<int, int, double>(t, r, it->first)] = it->second;
      }
    }
  }
  return result;
}

template <class DataPattern, class MatrixType>
std::map<int, double> wilson_plaket_schwinger_electric_transversal_l_even(
    const std::vector<MatrixType> &wilson_loops,
    const std::vector<MatrixType> &wilson_loops_opposite,
    const std::vector<MatrixType> &plaket_left,
    const std::vector<MatrixType> &plaket_right,
    const std::vector<std::vector<MatrixType>> &schwinger_lines,
    DataPattern &data_pattern, int d_max, int time, int r) {
  std::vector<double> correlator(2 * d_max + 1);
  MatrixType W;
  MatrixType A1;
  MatrixType A2;
  MatrixType S1;
  MatrixType S2;
  int d;
  int place;
#pragma omp parallel for collapse(4) private(W, A1, A2, S1, S2, d, place)      \
    firstprivate(d_max, data_pattern) reduction(vec_double_plus : correlator)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          for (int dir = 0; dir < 3; dir++) {
            data_pattern.lat_coord = {x, y, z, t};
            W = wilson_loops[data_pattern.get_lattice_size * 3 + dir];
            S1 = schwinger_lines[r / 2 - 2]
                                [data_pattern.get_lattice_size * 3 + dir];
            A1 = W * S1;
            data_pattern.move_forward(r / 2 - 1, dir);
            for (int mu = 0; mu < 3; mu++) {
              if (dir != mu) {
                data_pattern.move_backward(d_max, mu);
                for (int d = -d_max; d < 0; d++) {
                  S2 = schwinger_lines[abs(d) - 1]
                                      [data_pattern.get_lattice_size * 3 + mu];
                  A2 = A1 ^ S2;
                  A2 =
                      A2 * plaket_left[data_pattern.get_lattice_size * 3 + dir];
                  A2 = A2 * S2;
                  correlator[d + d_max] += A2.multiply_conj_tr(S1);
                  data_pattern.move_forward(1, mu);
                }
                A2 = A1 * plaket_left[data_pattern.get_lattice_size * 3 + dir];
                correlator[d_max] += A2.multiply_conj_tr(S1);
                for (int d = 1; d <= d_max; d++) {
                  S2 = schwinger_lines[d - 1]
                                      [data_pattern.get_lattice_size * 3 + mu];
                  data_pattern.move_forward(d, mu);
                  A2 = A1 * S2;
                  A2 =
                      A2 * plaket_left[data_pattern.get_lattice_size * 3 + dir];
                  A2 = A2 ^ S2;
                  correlator[d + d_max] += A2.multiply_conj_tr(S1);
                  data_pattern.move_backward(d, mu);
                }
              }
            }
            data_pattern.lat_coord = {x, y, z, t};
            data_pattern.move_forward(r, dir);
            W = wilson_loops_opposite[data_pattern.get_lattice_size * 3 + dir];
            data_pattern.move_backward(r / 2 - 1, dir);
            S1 = schwinger_lines[r / 2 - 2]
                                [data_pattern.get_lattice_size * 3 + dir]
                                    .conj();
            A1 = W * S1;
            for (int mu = 0; mu < 3; mu++) {
              if (dir != mu) {
                data_pattern.move_backward(d_max, mu);
                for (int d = -d_max; d < 0; d++) {
                  S2 = schwinger_lines[abs(d) - 1]
                                      [data_pattern.get_lattice_size * 3 + mu];
                  A2 = A1 ^ S2;
                  A2 = A2 *
                       plaket_right[data_pattern.get_lattice_size * 3 + dir];
                  A2 = A2 * S2;
                  correlator[d + d_max] += A2.multiply_conj_tr(S1);
                  data_pattern.move_forward(1, mu);
                }
                A2 = A1 * plaket_right[data_pattern.get_lattice_size * 3 + dir];
                correlator[d_max] += A2.multiply_conj_tr(S1);
                for (int d = 1; d <= d_max; d++) {
                  S2 = schwinger_lines[d - 1]
                                      [data_pattern.get_lattice_size * 3 + mu];
                  data_pattern.move_forward(d, mu);
                  A2 = A1 * S2;
                  A2 = A2 *
                       plaket_right[data_pattern.get_lattice_size * 3 + dir];
                  A2 = A2 ^ S2;
                  correlator[d + d_max] += A2.multiply_conj_tr(S1);
                  data_pattern.move_backward(d, mu);
                }
              }
            }
          }
        }
      }
    }
  }
  std::map<int, double> result;
  for (int i = 0; i < correlator.size(); i++) {
    result[i - d_max] = correlator[i] / (data_pattern.get_lattice_size * 12);
  }
  return result;
}

template <class DataPattern, class MatrixType>
std::map<int, double> wilson_plaket_schwinger_electric_transversal_l_odd(
    const std::vector<MatrixType> &wilson_loops,
    const std::vector<MatrixType> &plaket_left,
    const std::vector<MatrixType> &plaket_right,
    const std::vector<std::vector<MatrixType>> &schwinger_lines,
    DataPattern &data_pattern, int d_max, int time, int r) {
  std::vector<double> correlator(2 * d_max + 1);
  MatrixType W;
  MatrixType A1;
  MatrixType A2;
  MatrixType S1;
  MatrixType S2;
  int d;
  int place;
#pragma omp parallel for collapse(4) private(W, A1, A2, S1, S2, d, place)      \
    firstprivate(d_max, data_pattern) reduction(vec_double_plus : correlator)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          for (int dir = 0; dir < 3; dir++) {
            data_pattern.lat_coord = {x, y, z, t};
            W = wilson_loops[data_pattern.get_lattice_size * 3 + dir];
            S1 = schwinger_lines[r / 2 - 1]
                                [data_pattern.get_lattice_size * 3 + dir];
            A1 = W * S1;
            data_pattern.move_forward(r / 2, dir);
            for (int mu = 0; mu < 3; mu++) {
              if (dir != mu) {
                data_pattern.move_backward(d_max, mu);
                for (int d = -d_max; d < 0; d++) {
                  S2 = schwinger_lines[abs(d) - 1]
                                      [data_pattern.get_lattice_size * 3 + mu];
                  A2 = A1 ^ S2;
                  A2 =
                      A2 * plaket_left[data_pattern.get_lattice_size * 3 + dir];
                  A2 = A2 * S2;
                  correlator[d + d_max] += A2.multiply_conj_tr(S1);
                  data_pattern.move_forward(1, mu);
                }
                A2 = A1 * plaket_left[data_pattern.get_lattice_size * 3 + dir];
                correlator[d_max] += A2.multiply_conj_tr(S1);
                for (int d = 1; d <= d_max; d++) {
                  S2 = schwinger_lines[d - 1]
                                      [data_pattern.get_lattice_size * 3 + mu];
                  data_pattern.move_forward(d, mu);
                  A2 = A1 * S2;
                  A2 =
                      A2 * plaket_left[data_pattern.get_lattice_size * 3 + dir];
                  A2 = A2 ^ S2;
                  correlator[d + d_max] += A2.multiply_conj_tr(S1);
                  data_pattern.move_backward(d, mu);
                }
              }
            }
          }
        }
      }
    }
  }
  std::map<int, double> result;
  for (int i = 0; i < correlator.size(); i++) {
    result[i - d_max] = correlator[i] / (data_pattern.get_lattice_size * 6);
  }
  return result;
}

template <class DataPattern, class MatrixType>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_transversal_l(
    const Data::LatticeData<DataPattern, MatrixType> &conf,
    const std::vector<std::vector<MatrixType>> &schwinger_lines_short,
    int T_min, int T_max, int R_min, int R_max, int d_max) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> plaket_time_left =
      calculate_plaket_schwinger_time_left(conf);
  // std::vector<MatrixType> plaket_time_left =
  // calculate_plaket_time_left_up(conf); std::vector<MatrixType>
  // plaket_time_left = calculate_plaket_time_left_down(conf);
  std::vector<MatrixType> plaket_time_right =
      calculate_plaket_schwinger_time_right(conf);
  std::vector<MatrixType> wilson_loops;
  std::vector<MatrixType> wilson_loops_opposite;
  std::map<std::tuple<int, int, int>, double> result;
  std::map<int, double> schwinger_electric;
  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {
      wilson_loops = calculate_wilson_loops_schwinger(conf, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(conf, r, t);
      if (r % 2 == 0) {
        schwinger_electric =
            wilson_plaket_schwinger_electric_transversal_l_even(
                wilson_loops, wilson_loops_opposite, plaket_time_left,
                plaket_time_right, schwinger_lines_short, data_pattern, d_max,
                t, r);
      } else {
        schwinger_electric = wilson_plaket_schwinger_electric_transversal_l_odd(
            wilson_loops, plaket_time_left, plaket_time_right,
            schwinger_lines_short, data_pattern, d_max, t, r);
      }
      for (auto it = schwinger_electric.begin(); it != schwinger_electric.end();
           ++it) {
        result[std::tuple<int, int, double>(t, r, it->first)] = it->second;
      }
    }
  }
  return result;
}

template <class DataPattern, class MatrixType>
std::map<int, double> wilson_plaket_schwinger_electric_transversal_tr_even(
    const std::vector<MatrixType> &wilson_loops,
    const std::vector<MatrixType> &plaket,
    const std::vector<std::vector<MatrixType>> &schwinger_lines,
    DataPattern &data_pattern, int d_max, int time, int r) {
  std::vector<double> correlator(2 * d_max + 1);
  MatrixType W;
  MatrixType A1;
  MatrixType A2;
  MatrixType S1;
  MatrixType S2;
  int d;
  int place;
#pragma omp parallel for collapse(4) private(W, A1, A2, S1, S2, d, place)      \
    firstprivate(d_max, data_pattern) reduction(vec_double_plus : correlator)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          for (int dir = 0; dir < 3; dir++) {
            data_pattern.lat_coord = {x, y, z, t};
            W = wilson_loops[data_pattern.get_lattice_size * 3 + dir];
            S1 = schwinger_lines[r / 2 - 2]
                                [data_pattern.get_lattice_size * 3 + dir];
            A1 = W * S1;
            data_pattern.move_forward(r / 2 - 1, dir);
            for (int mu = 0; mu < 3; mu++) {
              if (dir != mu) {
                data_pattern.move_backward(d_max, mu);
                for (int d = -d_max; d < 0; d++) {
                  S2 = schwinger_lines[abs(d) - 1]
                                      [data_pattern.get_lattice_size * 3 + mu];
                  A2 = A1 ^ S2;
                  A2 = A2 * plaket[data_pattern.get_lattice_size * 3 + dir];
                  A2 = A2 * S2;
                  correlator[d + d_max] += A2.multiply_conj_tr(S1);
                  data_pattern.move_forward(1, mu);
                }
                A2 = A1 * plaket[data_pattern.get_lattice_size * 3 + dir];
                correlator[d_max] += A2.multiply_conj_tr(S1);
                for (int d = 1; d <= d_max; d++) {
                  S2 = schwinger_lines[d - 1]
                                      [data_pattern.get_lattice_size * 3 + mu];
                  data_pattern.move_forward(d, mu);
                  A2 = A1 * S2;
                  A2 = A2 * plaket[data_pattern.get_lattice_size * 3 + dir];
                  A2 = A2 ^ S2;
                  correlator[d + d_max] += A2.multiply_conj_tr(S1);
                  data_pattern.move_backward(d, mu);
                }
              }
            }
          }
        }
      }
    }
  }
  std::map<int, double> result;
  for (int i = 0; i < correlator.size(); i++) {
    result[i - d_max] = correlator[i] / (data_pattern.get_lattice_size * 6);
  }
  return result;
}

template <class DataPattern, class MatrixType>
std::map<int, double> wilson_plaket_schwinger_electric_transversal_tr_odd(
    const std::vector<MatrixType> &wilson_loops,
    const std::vector<MatrixType> &wilson_loops_opposite,
    const std::vector<MatrixType> &plaket,
    const std::vector<std::vector<MatrixType>> &schwinger_lines,
    DataPattern &data_pattern, int d_max, int time, int r) {
  std::vector<double> correlator(2 * d_max + 1);
  MatrixType W;
  MatrixType A1;
  MatrixType A2;
  MatrixType S1;
  MatrixType S2;
  int d;
  int place;
#pragma omp parallel for collapse(4) private(W, A1, A2, S1, S2, d, place)      \
    firstprivate(d_max, data_pattern) reduction(vec_double_plus : correlator)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          for (int dir = 0; dir < 3; dir++) {
            data_pattern.lat_coord = {x, y, z, t};
            W = wilson_loops[data_pattern.get_lattice_size * 3 + dir];
            S1 = schwinger_lines[r / 2 - 1]
                                [data_pattern.get_lattice_size * 3 + dir];
            A1 = W * S1;
            data_pattern.move_forward(r / 2, dir);
            for (int mu = 0; mu < 3; mu++) {
              if (dir != mu) {
                data_pattern.move_backward(d_max, mu);
                for (int d = -d_max; d < 0; d++) {
                  S2 = schwinger_lines[abs(d) - 1]
                                      [data_pattern.get_lattice_size * 3 + mu];
                  A2 = A1 ^ S2;
                  A2 = A2 * plaket[data_pattern.get_lattice_size * 3 + dir];
                  A2 = A2 * S2;
                  correlator[d + d_max] += A2.multiply_conj_tr(S1);
                  data_pattern.move_forward(1, mu);
                }
                A2 = A1 * plaket[data_pattern.get_lattice_size * 3 + dir];
                correlator[d_max] += A2.multiply_conj_tr(S1);
                for (int d = 1; d <= d_max; d++) {
                  S2 = schwinger_lines[d - 1]
                                      [data_pattern.get_lattice_size * 3 + mu];
                  data_pattern.move_forward(d, mu);
                  A2 = A1 * S2;
                  A2 = A2 * plaket[data_pattern.get_lattice_size * 3 + dir];
                  A2 = A2 ^ S2;
                  correlator[d + d_max] += A2.multiply_conj_tr(S1);
                  data_pattern.move_backward(d, mu);
                }
              }
            }
            data_pattern.lat_coord = {x, y, z, t};
            data_pattern.move_forward(r, dir);
            W = wilson_loops_opposite[data_pattern.get_lattice_size * 3 + dir];
            data_pattern.move_backward(r / 2 - 1, dir);
            S1 = schwinger_lines[r / 2 - 2]
                                [data_pattern.get_lattice_size * 3 + dir]
                                    .conj();
            A1 = W * S1;
            for (int mu = 0; mu < 3; mu++) {
              if (dir != mu) {
                data_pattern.move_backward(d_max, mu);
                for (int d = -d_max; d < 0; d++) {
                  S2 = schwinger_lines[abs(d) - 1]
                                      [data_pattern.get_lattice_size * 3 + mu];
                  A2 = A1 ^ S2;
                  A2 = A2 * plaket[data_pattern.get_lattice_size * 3 + dir];
                  A2 = A2 * S2;
                  correlator[d + d_max] += A2.multiply_conj_tr(S1);
                  data_pattern.move_forward(1, mu);
                }
                A2 = A1 * plaket[data_pattern.get_lattice_size * 3 + dir];
                correlator[d_max] += A2.multiply_conj_tr(S1);
                for (int d = 1; d <= d_max; d++) {
                  S2 = schwinger_lines[d - 1]
                                      [data_pattern.get_lattice_size * 3 + mu];
                  data_pattern.move_forward(d, mu);
                  A2 = A1 * S2;
                  A2 = A2 * plaket[data_pattern.get_lattice_size * 3 + dir];
                  A2 = A2 ^ S2;
                  correlator[d + d_max] += A2.multiply_conj_tr(S1);
                  data_pattern.move_backward(d, mu);
                }
              }
            }
          }
        }
      }
    }
  }
  std::map<int, double> result;
  for (int i = 0; i < correlator.size(); i++) {
    result[i - d_max] = correlator[i] / (data_pattern.get_lattice_size * 12);
  }
  return result;
}

template <class DataPattern, class MatrixType>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_transversal_tr(
    const Data::LatticeData<DataPattern, MatrixType> &conf,
    const std::vector<std::vector<MatrixType>> &schwinger_lines_short,
    int T_min, int T_max, int R_min, int R_max, int d_max) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> plaket_time =
      calculate_plaket_schwinger_time_tr(conf);
  std::vector<MatrixType> wilson_loops;
  std::vector<MatrixType> wilson_loops_opposite;
  std::map<std::tuple<int, int, int>, double> result;
  std::map<int, double> schwinger_electric;
  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {
      wilson_loops = calculate_wilson_loops_schwinger(conf, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(conf, r, t);
      if (r % 2 == 0) {
        schwinger_electric =
            wilson_plaket_schwinger_electric_transversal_tr_even(
                wilson_loops, plaket_time, schwinger_lines_short, data_pattern,
                d_max, t, r);
      } else {
        schwinger_electric =
            wilson_plaket_schwinger_electric_transversal_tr_odd(
                wilson_loops, wilson_loops_opposite, plaket_time,
                schwinger_lines_short, data_pattern, d_max, t, r);
      }
      for (auto it = schwinger_electric.begin(); it != schwinger_electric.end();
           ++it) {
        result[std::tuple<int, int, double>(t, r, it->first)] = it->second;
      }
    }
  }
  return result;
}

template <class DataPattern, class MatrixType>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_longitudinal_l(
    const Data::LatticeData<DataPattern, MatrixType> &conf,
    const std::vector<std::vector<MatrixType>> &schwinger_lines_short,
    int T_min, int T_max, int R_min, int R_max, int d_ouside) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> plaket_space =
      calculate_plaket_schwinger_space_l(conf);
  std::vector<MatrixType> wilson_loops;
  std::vector<MatrixType> wilson_loops_opposite;
  std::map<std::tuple<int, int, int>, double> result;
  std::map<int, double> schwinger_electric;
  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {
      wilson_loops = calculate_wilson_loops_schwinger(conf, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(conf, r, t);
      schwinger_electric = wilson_plaket_schwinger_longitudinal_tr(
          wilson_loops, wilson_loops_opposite, plaket_space,
          schwinger_lines_short, data_pattern, d_ouside, t, r);
      for (auto it = schwinger_electric.begin(); it != schwinger_electric.end();
           ++it) {
        result[std::tuple<int, int, double>(t, r, it->first)] = it->second;
      }
    }
  }
  return result;
}

template <class DataPattern, class MatrixType>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_longitudinal_tr(
    const Data::LatticeData<DataPattern, MatrixType> &conf,
    const std::vector<std::vector<MatrixType>> &schwinger_lines_short,
    int T_min, int T_max, int R_min, int R_max, int d_ouside) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> plaket_space =
      calculate_plaket_schwinger_space_tr(conf);
  std::vector<MatrixType> wilson_loops;
  std::vector<MatrixType> wilson_loops_opposite;
  std::map<std::tuple<int, int, int>, double> result;
  std::map<int, double> schwinger_electric;
  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {
      wilson_loops = calculate_wilson_loops_schwinger(conf, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(conf, r, t);
      schwinger_electric = wilson_plaket_schwinger_longitudinal_tr(
          wilson_loops, wilson_loops_opposite, plaket_space,
          schwinger_lines_short, data_pattern, d_ouside, t, r);
      for (auto it = schwinger_electric.begin(); it != schwinger_electric.end();
           ++it) {
        result[std::tuple<int, int, double>(t, r, it->first)] = it->second;
      }
    }
  }
  return result;
}

template <class DataPattern, class MatrixType>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_transversal_l(
    const Data::LatticeData<DataPattern, MatrixType> &conf,
    const std::vector<std::vector<MatrixType>> &schwinger_lines_short,
    int T_min, int T_max, int R_min, int R_max, int d_max) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> plaket_space =
      calculate_plaket_schwinger_space_l(conf);
  std::vector<MatrixType> wilson_loops;
  std::vector<MatrixType> wilson_loops_opposite;
  std::map<std::tuple<int, int, int>, double> result;
  std::map<int, double> schwinger_electric;
  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {
      wilson_loops = calculate_wilson_loops_schwinger(conf, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(conf, r, t);
      if (r % 2 == 0) {
        schwinger_electric =
            wilson_plaket_schwinger_electric_transversal_tr_even(
                wilson_loops, plaket_space, schwinger_lines_short, data_pattern,
                d_max, t, r);
      } else {
        schwinger_electric =
            wilson_plaket_schwinger_electric_transversal_tr_odd(
                wilson_loops, wilson_loops_opposite, plaket_space,
                schwinger_lines_short, data_pattern, d_max, t, r);
      }
      for (auto it = schwinger_electric.begin(); it != schwinger_electric.end();
           ++it) {
        result[std::tuple<int, int, double>(t, r, it->first)] = it->second;
      }
    }
  }
  return result;
}

template <class DataPattern, class MatrixType>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_transversal_tr(
    const Data::LatticeData<DataPattern, MatrixType> &conf,
    const std::vector<std::vector<MatrixType>> &schwinger_lines_short,
    int T_min, int T_max, int R_min, int R_max, int d_max) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> plaket_space =
      calculate_plaket_schwinger_space_tr(conf);
  std::vector<MatrixType> wilson_loops;
  std::vector<MatrixType> wilson_loops_opposite;
  std::map<std::tuple<int, int, int>, double> result;
  std::map<int, double> schwinger_electric;
  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {
      wilson_loops = calculate_wilson_loops_schwinger(conf, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(conf, r, t);
      if (r % 2 == 0) {
        schwinger_electric =
            wilson_plaket_schwinger_electric_transversal_tr_even(
                wilson_loops, plaket_space, schwinger_lines_short, data_pattern,
                d_max, t, r);
      } else {
        schwinger_electric =
            wilson_plaket_schwinger_electric_transversal_tr_odd(
                wilson_loops, wilson_loops_opposite, plaket_space,
                schwinger_lines_short, data_pattern, d_max, t, r);
      }
      for (auto it = schwinger_electric.begin(); it != schwinger_electric.end();
           ++it) {
        result[std::tuple<int, int, double>(t, r, it->first)] = it->second;
      }
    }
  }
  return result;
}

template <class DataPattern, class MatrixType>
void flux_schwinger_all(
    const Data::LatticeData<DataPattern, MatrixType> &conf,
    const std::vector<std::vector<MatrixType>> &schwinger_lines_short,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_electric_long_l,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_electric_long_tr,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_electric_trans_l,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_electric_trans_tr,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_magnetic_long_l,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_magnetic_long_tr,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_magnetic_trans_l,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_magnetic_trans_tr,
    int T_min, int T_max, int R_min, int R_max, int d_ouside, int d_max) {

  std::vector<MatrixType> plaket_time_left =
      calculate_plaket_schwinger_time_left(conf);
  std::vector<MatrixType> plaket_time_right =
      calculate_plaket_schwinger_time_right(conf);
  std::vector<MatrixType> plaket_time_tr =
      calculate_plaket_schwinger_time_tr(conf);
  std::vector<MatrixType> plaket_space_l =
      calculate_plaket_schwinger_space_l(conf);
  std::vector<MatrixType> plaket_space_tr =
      calculate_plaket_schwinger_space_tr(conf);

  std::vector<MatrixType> wilson_loops;
  std::vector<MatrixType> wilson_loops_opposite;
  std::map<int, double> schwinger_correlator;

  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {

      wilson_loops = calculate_wilson_loops_schwinger(conf, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(conf, r, t);

      schwinger_correlator = wilson_plaket_schwinger_longitudinal_l(
          wilson_loops, wilson_loops_opposite, plaket_time_left,
          plaket_time_right, schwinger_lines_short, d_ouside, t, r);
      for (auto it = schwinger_correlator.begin();
           it != schwinger_correlator.end(); ++it) {
        flux_tube_schwinger_electric_long_l[std::tuple<int, int, double>(
            t, r, it->first)] = it->second;
      }

      schwinger_correlator = wilson_plaket_schwinger_longitudinal_tr(
          wilson_loops, wilson_loops_opposite, plaket_time_tr,
          schwinger_lines_short, d_ouside, t, r);
      for (auto it = schwinger_correlator.begin();
           it != schwinger_correlator.end(); ++it) {
        flux_tube_schwinger_electric_long_tr[std::tuple<int, int, double>(
            t, r, it->first)] = it->second;
      }

      if (r % 2 == 0) {
        schwinger_correlator =
            wilson_plaket_schwinger_electric_transversal_l_even(
                wilson_loops, wilson_loops_opposite, plaket_time_left,
                plaket_time_right, schwinger_lines_short, d_max, t, r);
      } else {
        schwinger_correlator =
            wilson_plaket_schwinger_electric_transversal_l_odd(
                wilson_loops, plaket_time_left, plaket_time_right,
                schwinger_lines_short, d_max, t, r);
      }
      for (auto it = schwinger_correlator.begin();
           it != schwinger_correlator.end(); ++it) {
        flux_tube_schwinger_electric_trans_l[std::tuple<int, int, double>(
            t, r, it->first)] = it->second;
      }

      if (r % 2 == 0) {
        schwinger_correlator =
            wilson_plaket_schwinger_electric_transversal_tr_even(
                wilson_loops, plaket_time_tr, schwinger_lines_short, d_max, t,
                r);
      } else {
        schwinger_correlator =
            wilson_plaket_schwinger_electric_transversal_tr_odd(
                wilson_loops, wilson_loops_opposite, plaket_time_tr,
                schwinger_lines_short, d_max, t, r);
      }
      for (auto it = schwinger_correlator.begin();
           it != schwinger_correlator.end(); ++it) {
        flux_tube_schwinger_electric_trans_tr[std::tuple<int, int, double>(
            t, r, it->first)] = it->second;
      }

      schwinger_correlator = wilson_plaket_schwinger_longitudinal_tr(
          wilson_loops, wilson_loops_opposite, plaket_space_l,
          schwinger_lines_short, d_ouside, t, r);
      for (auto it = schwinger_correlator.begin();
           it != schwinger_correlator.end(); ++it) {
        flux_tube_schwinger_magnetic_long_l[std::tuple<int, int, double>(
            t, r, it->first)] = it->second;
      }

      schwinger_correlator = wilson_plaket_schwinger_longitudinal_tr(
          wilson_loops, wilson_loops_opposite, plaket_space_tr,
          schwinger_lines_short, d_ouside, t, r);
      for (auto it = schwinger_correlator.begin();
           it != schwinger_correlator.end(); ++it) {
        flux_tube_schwinger_magnetic_long_tr[std::tuple<int, int, double>(
            t, r, it->first)] = it->second;
      }

      if (r % 2 == 0) {
        schwinger_correlator =
            wilson_plaket_schwinger_electric_transversal_tr_even(
                wilson_loops, plaket_space_l, schwinger_lines_short, d_max, t,
                r);
      } else {
        schwinger_correlator =
            wilson_plaket_schwinger_electric_transversal_tr_odd(
                wilson_loops, wilson_loops_opposite, plaket_space_l,
                schwinger_lines_short, d_max, t, r);
      }
      for (auto it = schwinger_correlator.begin();
           it != schwinger_correlator.end(); ++it) {
        flux_tube_schwinger_magnetic_trans_l[std::tuple<int, int, double>(
            t, r, it->first)] = it->second;
      }

      if (r % 2 == 0) {
        schwinger_correlator =
            wilson_plaket_schwinger_electric_transversal_tr_even(
                wilson_loops, plaket_space_tr, schwinger_lines_short, d_max, t,
                r);
      } else {
        schwinger_correlator =
            wilson_plaket_schwinger_electric_transversal_tr_odd(
                wilson_loops, wilson_loops_opposite, plaket_space_tr,
                schwinger_lines_short, d_max, t, r);
      }
      for (auto it = schwinger_correlator.begin();
           it != schwinger_correlator.end(); ++it) {
        flux_tube_schwinger_magnetic_trans_tr[std::tuple<int, int, double>(
            t, r, it->first)] = it->second;
      }
    }
  }
}

// d_left and d_right are non-negative
// d_left + d_right + r + 1 values starting from d_left sites to the negative mu
// directon from wilson loop and ending at d_right sites to the positive mu
// direction from wilson_loop
template <class DataPattern>
void wilson_plaket_correlator_plane_longitudinal(
    std::vector<double> &correlator, const std::vector<double> &wilson_loop_tr,
    const std::vector<double> &plaket_tr, DataPattern &data_pattern, int d_left,
    int d_right, int r, int t, int mu) {
  int index_wilson;
#pragma omp parallel for collapse(4) firstprivate(                             \
        data_pattern, d_left, d_right, r, t, mu) private(index_wilson)         \
    reduction(vec_double_plus : correlator)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index_wilson = data_pattern.get_index_site();
          data_pattern.move_forward(t / 2, 3);
          data_pattern.move_backward(d_left, mu);
          for (int d = 0; d < d_left + d_right + r + 1; d++) {
            correlator[d] += wilson_loop_tr[index_wilson] *
                             plaket_tr[data_pattern.get_index_site()];
            data_pattern.move_forward(1, mu);
          }
        }
      }
    }
  }
}

// d_bot and d_top are non-negative
template <class DataPattern>
void wilson_plaket_correlator_plane_transversal(
    std::vector<double> &correlator, const std::vector<double> &wilson_loop_tr,
    const std::vector<double> &plaket_tr, DataPattern &data_pattern, int d_bot,
    int d_top, int r, int t, int mu) {
  int index_wilson;
#pragma omp parallel for collapse(4)                                           \
    firstprivate(data_pattern, d_bot, d_top, r, t, mu) private(index_wilson)   \
    reduction(vec_double_plus : correlator)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index_wilson = data_pattern.get_index_site();
          data_pattern.move_forward(t / 2, 3);
          data_pattern.move_forward(r / 2, mu);
          for (int nu = 0; nu < 3; nu++) {
            if (nu != mu) {
              data_pattern.move_backward(d_bot, nu);
              for (int d = 0; d < d_bot; d++) {
                correlator[d] += wilson_loop_tr[index_wilson] *
                                 plaket_tr[data_pattern.get_index_site()];
                data_pattern.move_forward(1, nu);
              }
              correlator[d_bot] += wilson_loop_tr[index_wilson] *
                                   plaket_tr[data_pattern.get_index_site()];
              for (int d = d_bot + 1; d < d_bot + d_top + 1; d++) {
                data_pattern.move_forward(1, nu);
                correlator[d] += wilson_loop_tr[index_wilson] *
                                 plaket_tr[data_pattern.get_index_site()];
              }
              data_pattern.move_backward(d_bot, nu);
            }
          }
        }
      }
    }
  }
}

template <class DataPattern, class MatrixType>
std::map<std::tuple<int, int, int>, double> wilson_plaket_correlator(
    const std::vector<double> &plaket_tr,
    const Data::LatticeData<DataPattern, MatrixType> &conf_wilson, int T_min,
    int T_max, int R_min, int R_max, int main_coordinate,
    int transverse_coordinate, std::string direction) {
  DataPattern data_pattern(conf_wilson.lat_dim);
  std::map<std::tuple<int, int, int>, double> flux_tube;
  std::map<int, std::vector<MatrixType>> time_lines;
  std::vector<MatrixType> space_lines;
  for (int t = T_min; t <= T_max; t += 2) {
    time_lines[t] = wilson_lines_indexed(conf_wilson, t, 3);
  }
  int transverse_shift = transverse_coordinate;
  std::vector<double> wilson_tr;
  std::vector<double> correlator;
  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 2) {
      int main_coordinate_min, main_coordinate_max;
      if (direction.compare("longitudinal") == 0) {
        main_coordinate_min = main_coordinate;
        main_coordinate_max = main_coordinate;
        correlator = std::vector<double>(main_coordinate_max +
                                         main_coordinate_min + r + 1);
      } else if (direction.compare("transversal") == 0) {
        main_coordinate_min = main_coordinate;
        main_coordinate_max = main_coordinate;
        transverse_coordinate = r / 2 + transverse_shift;
        correlator =
            std::vector<double>(main_coordinate_max + main_coordinate_min + 1);
      } else {
        std::cout << "wrong direction" << std::endl;
        exit(1);
      }
      for (int mu = 0; mu < 3; mu++) {
        space_lines = wilson_lines_indexed(conf_wilson, r, mu);
        wilson_tr = calculate_wilson_loop_plane_tr(space_lines, time_lines[t],
                                                   data_pattern, mu, 3, r, t);
        if (direction.compare("longitudinal") == 0) {
          wilson_plaket_correlator_plane_longitudinal(
              correlator, wilson_tr, plaket_tr, data_pattern,
              main_coordinate_min, main_coordinate_max, r, t, mu);
        } else if (direction.compare("transversal") == 0) {
          wilson_plaket_correlator_plane_transversal(
              correlator, wilson_tr, plaket_tr, data_pattern,
              main_coordinate_min, main_coordinate_max, r, t, mu);
        }
      }
      if (direction.compare("longitudinal") == 0) {
        for (int D = 0; D <= r / 2 + main_coordinate; D++) {
          flux_tube[std::tuple<int, int, int>(t, r, D)] =
              (correlator[D + main_coordinate + r / 2] +
               correlator[-D + main_coordinate + r / 2]) /
              data_pattern.get_lattice_size() / 3 / 2;
        }
      } else if (direction.compare("transversal") == 0) {
        for (int D = 0; D <= main_coordinate_max; D++) {
          flux_tube[std::tuple<int, int, int>(t, r, D)] =
              (correlator[D - main_coordinate_min] +
               correlator[-D - main_coordinate_min]) /
              data_pattern.get_lattice_size() / 6 / 2;
        }
      }
    }
  }
  return flux_tube;
}

// Wilson_plaket_schwinger_correlator
template <class T>
std::vector<T> calculate_plaket_time_left_down(const std::vector<T> &array);
template <class T>
std::vector<T> calculate_plaket_time_left_up(const std::vector<T> &array);
template <class T>
std::vector<T> calculate_plaket_time_right_down(const std::vector<T> &array);
template <class T>
std::vector<T> calculate_plaket_time_right_up(const std::vector<T> &array);
template <class T>
std::vector<T> calculate_plaket_space(const std::vector<T> &array);
// template <class T>
// std::vector<T> calculate_polyakov_loop(const std::vector<T> &array);
// template <class T>
// std::vector<T> calculate_wilson_loop(const std::vector<T> &array, int r,
//                                      int time);
template <class T>
std::vector<T> calculate_plaket_schwinger_time(const std::vector<T> &array);
template <class T>
std::vector<std::vector<T>>
calculate_plaket_schwinger_space(const std::vector<T> &array);
template <class T>
std::vector<T> calculate_schwinger_lines_short(const std::vector<T> &array,
                                               int d);
template <class T>
std::vector<std::vector<T>>
calculate_schwinger_line(const std::vector<T> &array, int d, int x_trans);
template <class T>
std::vector<T> calculate_wilson_loops_schwinger(const std::vector<T> &array,
                                                int r, int time);
template <class T>
std::vector<T>
calculate_wilson_loops_schwinger_opposite(const std::vector<T> &array, int r,
                                          int time);
template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_longitudinal_l(
    const std::vector<T> &array,
    const std::vector<std::vector<T>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_ouside);
template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_longitudinal_tr(
    const std::vector<T> &array,
    const std::vector<std::vector<T>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_ouside);
template <class T>
double schwinger_electric_long_tr_even_test(const std::vector<T> &conf, int t,
                                            int r);
template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_transversal_l(
    const std::vector<T> &array,
    const std::vector<std::vector<T>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_max);
template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_transversal_tr(
    const std::vector<T> &array,
    const std::vector<std::vector<T>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_max);
template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_longitudinal_l(
    const std::vector<T> &array,
    const std::vector<std::vector<T>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_ouside);
template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_longitudinal_tr(
    const std::vector<T> &array,
    const std::vector<std::vector<T>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_ouside);
template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_transversal_l(
    const std::vector<T> &array,
    const std::vector<std::vector<T>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_max);
template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_transversal_tr(
    const std::vector<T> &array,
    const std::vector<std::vector<T>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_max);
template <class T>
void flux_schwinger_all(
    const std::vector<T> &array,
    const std::vector<std::vector<T>> &schwinger_lines_short,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_electric_long_l,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_electric_long_tr,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_electric_trans_l,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_electric_trans_tr,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_magnetic_long_l,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_magnetic_long_tr,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_magnetic_trans_l,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_magnetic_trans_tr,
    int T_min, int T_max, int R_min, int R_max, int d_ouside, int d_max);

// Wilson_plaket_correlator
template <class T>
std::vector<double> calculate_plaket_time_trace_l(const std::vector<T> &array);
template <class T>
std::vector<double> calculate_plaket_time_trace_tr(const std::vector<T> &array);
template <class T>
std::vector<double> calculate_plaket_space_trace_l(const std::vector<T> array);
template <class T>
std::vector<double>
calculate_plaket_space_trace_tr(const std::vector<T> &array);
template <class T>
std::vector<double> calculate_plaket_space_tr(const std::vector<T> &array);
double plaket4_time(const std::vector<double> &plaket_tr, link1 &link);
double plaket4_space(const std::vector<double> &plaket_tr, link1 &link, int nu);

template <class T>
std::vector<double> calculate_wilson_loop_tr(const std::vector<T> &array, int r,
                                             int time);
template <class T>
std::map<std::tuple<int, int, int>, double>
calculate_wilson_plaket_correlator_electric_longitudinal(
    const std::vector<double> &plaket_tr, const std::vector<T> &conf_wilson,
    int T_min, int T_max, int R_min, int R_max, int d_ouside);
template <class T>
std::map<std::tuple<int, int, int>, double>
calculate_wilson_plaket_correlator_electric_transversal(
    const std::vector<double> &plaket_tr, const std::vector<T> &conf_wilson,
    int T_min, int T_max, int R_min, int R_max, int d_ouside);

std::map<int, double>
wilson_plaket_correlator_electric(const std::vector<double> &wilson_loop_tr,
                                  const std::vector<double> &plaket_tr, int r,
                                  int time, int x_trans, int d_min, int d_max);
std::map<int, double>
wilson_plaket_correlator_electric_x(const std::vector<double> &wilson_loop_tr,
                                    const std::vector<double> &plaket_tr, int r,
                                    int time, int x_trans_max, int d);
std::map<int, double>
wilson_plaket_correlator_magnetic(const std::vector<double> &wilson_loop_tr,
                                  const std::vector<double> &plaket_tr, int r,
                                  int time, int x_trans, int d_min, int d_max);
std::map<int, double>
wilson_plaket_correlator_magnetic_x(const std::vector<double> &wilson_loop_tr,
                                    const std::vector<double> &plaket_tr, int r,
                                    int time, int x_trans_max, int d);

template <class T>
std::vector<double> wilson_plane_tr(const std::vector<T> &wilson_lines_mu,
                                    const std::vector<T> &wilson_lines_nu,
                                    int size_mu1, int size_mu2, int size_nu1,
                                    int size_nu2, int length_mu, int length_nu);

template <class T>
std::vector<double> plaket_plane_tr(const std::vector<T> &conf_mu,
                                    const std::vector<T> &conf_nu, int size_mu1,
                                    int size_mu2, int size_nu1, int size_nu2);

void plaket_plane_aver(std::vector<double> &plaket_aver_tr,
                       const std::vector<double> &plaket_tr, int size_mu1,
                       int size_mu2, int size_nu1, int size_nu2);

template <class T>
std::vector<double>
plaket_aver_tr_time(const std::vector<std::vector<T>> &conf);

template <class T>
std::vector<double>
plaket_aver_tr_space(const std::vector<std::vector<T>> &conf);

void wilson_plaket_correlator_plane_longitudinal(
    std::vector<double> &correlator, const std::vector<double> &wilson_loop_tr,
    const std::vector<double> &plaket_tr, int size_mu1, int size_mu2,
    int x_trans, int d_min, int d_max);

void wilson_plaket_correlator_plane_transversal(
    std::vector<double> &correlator, const std::vector<double> &wilson_loop_tr,
    const std::vector<double> &plaket_tr, int size_mu1, int size_mu2,
    int size_nu1, int size_nu2, int d, int x_trans_min, int x_trans_max);

template <class T>
std::vector<double> plaket_tr_single_dir_test(const std::vector<T> &conf,
                                              int mu);

template <class T>
std::map<std::tuple<int, int, int>, double>
wilson_plaket_correlator(const std::vector<double> &plaket_tr,
                         const std::vector<std::vector<T>> &conf_wilson,
                         int T_min, int T_max, int R_min, int R_max,
                         int main_coordinate, int transverse_coordinate,
                         std::string direction);