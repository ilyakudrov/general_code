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
  std::vector<MatrixType> vec(data_pattern.get_lattice_size() * 3);
#pragma omp parallel for collapse(4) firstprivate(data_pattern)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          for (int mu = 0; mu < 3; mu++) {
            vec[data_pattern.get_index_site() * 3 + mu] =
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
  std::vector<MatrixType> vec(data_pattern.get_lattice_size() * 3);
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
  std::vector<MatrixType> vec(data_pattern.get_lattice_size() * 3);
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
                wilson_loop_schwinger_opposite(conf, data_pattern, r, time, mu);
          }
        }
      }
    }
  }
  return vec;
}

template <class DataPattern>
std::map<int, double> wilson_plaket_correlator_electric_longitudinal(
    DataPattern &data_pattern, const std::vector<double> &wilson_loop_tr,
    const std::vector<double> &plaket_tr, int r, int time, int d_outside) {
  std::vector<double> correlator(2 * d_outside + r);
  double W;
#pragma omp parallel for collapse(4) private(W)                                \
    firstprivate(d_outside, data_pattern)                                      \
    reduction(vec_double_plus : correlator)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          for (int dir = 0; dir < 3; dir++) {
            data_pattern.lat_coord = {x, y, z, t};
            W = wilson_loop_tr[data_pattern.get_index_site() * 3 + dir];
            data_pattern.move_forward(time / 2, 3);
            data_pattern.move_backward(d_outside, dir);
            for (int d = 0; d < 2 * d_outside + r; d++) {
              correlator[d] +=
                  W * plaket_tr[data_pattern.get_index_site() * 3 + dir];
              data_pattern.move_forward(1, dir);
            }
          }
        }
      }
    }
  }
  std::map<int, double> result;
  for (int i = 0; i < correlator.size(); i++) {
    result[i] = correlator[i] / (data_pattern.get_lattice_size() * 3);
  }
  return result;
}

template <class DataPattern, class MatrixType>
std::map<std::tuple<int, int, int>, double>
calculate_wilson_plaket_correlator_electric_longitudinal(
    const std::vector<double> &plaket_tr,
    const Data::LatticeData<DataPattern, MatrixType> &conf, int T_min,
    int T_max, int R_min, int R_max, int d_ouside) {
  DataPattern data_pattern(conf.lat_dim);
  std::map<std::tuple<int, int, int>, double> flux_tube;
  std::vector<double> wilson_loop_tr;
  std::map<int, double> correlator;
  for (int time = T_min; time <= T_max; time += 2) {
    for (int r = R_min; r <= R_max; r += 1) {
      wilson_loop_tr = calculate_wilson_loop_time_tr(conf, r, time);
      correlator = wilson_plaket_correlator_electric_longitudinal(
          data_pattern, wilson_loop_tr, plaket_tr, r, time, d_ouside);
      for (auto i = correlator.begin(); i != correlator.end(); i++) {
        flux_tube[std::tuple<int, int, int>(time, r, i->first)] = i->second;
      }
    }
  }

  return flux_tube;
}

template <class DataPattern, class MatrixType>
std::map<int, double> wilson_plaket_schwinger_l_dir_l_plaket(
    const std::vector<MatrixType> &wilson_loops,
    const std::vector<MatrixType> &wilson_loops_opposite,
    const std::vector<MatrixType> &plaket_left,
    const std::vector<MatrixType> &plaket_right,
    const std::vector<std::vector<MatrixType>> &schwinger_lines,
    DataPattern &data_pattern, int d_outside, int time, int r) {
  int correlator_size = 2 * d_outside + r;
  std::vector<double> correlator(correlator_size);
  MatrixType W;
  MatrixType A;
  MatrixType S;
  int d; // place of plaket in correlator vector
  int index_schwinger;
  int place;
#pragma omp parallel for collapse(4) private(W, A, S, d, place,                \
                                                 index_schwinger)              \
    firstprivate(d_outside, correlator_size, data_pattern)                     \
    reduction(vec_double_plus : correlator)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          for (int dir = 0; dir < 3; dir++) {
            data_pattern.lat_coord = {x, y, z, t};
            W = wilson_loops[data_pattern.get_index_site() * 3 + dir];
            d = 0;
            index_schwinger = d_outside - 2;
            data_pattern.move_backward(d_outside - 1, dir);
            while (d < d_outside - 1) {
              place = data_pattern.get_index_site() * 3 + dir;
              S = schwinger_lines[index_schwinger][place];
              A = W ^ S;
              A = A ^ plaket_right[place];
              correlator[d] += A.multiply_tr(S);
              data_pattern.move_forward(1, dir);
              d++;
              index_schwinger--;
            }
            place = data_pattern.get_index_site() * 3;
            correlator[d] += W.multiply_conj_tr(plaket_right[place + dir]);
            d++;
            correlator[d] += W.multiply_tr(plaket_left[place + dir]);
            d++;
            index_schwinger = 0;
            while (d < d_outside + (r + 1) / 2) {
              S = schwinger_lines[index_schwinger]
                                 [data_pattern.get_index_site() * 3 + dir];
              A = W * S;
              data_pattern.move_forward(index_schwinger + 1, dir);
              A = A * plaket_left[data_pattern.get_index_site() * 3 + dir];
              correlator[d] += A.multiply_conj_tr(S);
              data_pattern.move_backward(index_schwinger + 1, dir);
              d++;
              index_schwinger++;
            }
            data_pattern.move_forward(r, dir);
            W = wilson_loops_opposite[data_pattern.get_index_site() * 3 + dir];
            data_pattern.move_backward(r / 2 - 1, dir);
            index_schwinger = r / 2 - 2;
            while (d < d_outside + r - 1) {
              place = data_pattern.get_index_site() * 3 + dir;
              S = schwinger_lines[index_schwinger][place];
              A = W ^ S;
              A = A * plaket_right[place];
              correlator[d] += A.multiply_tr(S);
              data_pattern.move_forward(1, dir);
              d++;
              index_schwinger--;
            }
            place = data_pattern.get_index_site() * 3 + dir;
            correlator[d] += W.multiply_tr(plaket_right[place]);
            d++;
            correlator[d] += W.multiply_conj_tr(plaket_left[place]);
            d++;
            index_schwinger = 0;
            while (d < r + 2 * d_outside) {
              place = data_pattern.get_index_site() * 3 + dir;
              S = schwinger_lines[index_schwinger]
                                 [data_pattern.get_index_site() * 3 + dir];
              A = W * S;
              data_pattern.move_forward(index_schwinger + 1, dir);
              A = A ^ plaket_left[data_pattern.get_index_site() * 3 + dir];
              correlator[d] += A.multiply_conj_tr(S);
              data_pattern.move_backward(index_schwinger + 1, dir);
              d++;
              index_schwinger++;
            }
          }
        }
      }
    }
  }
  std::map<int, double> result;
  for (int i = 0; i < correlator.size(); i++) {
    result[i] = correlator[i] / (data_pattern.get_lattice_size() * 3);
  }
  return result;
}

template <class DataPattern, class MatrixType>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_l_dir_l_plaket(
    const Data::LatticeData<DataPattern, MatrixType> &conf_plaket,
    const Data::LatticeData<DataPattern, MatrixType> &conf_wilson,
    const std::vector<std::vector<MatrixType>> &schwinger_lines_short,
    int T_min, int T_max, int R_min, int R_max, int d_ouside) {
  DataPattern data_pattern(conf_wilson.lat_dim);
  std::vector<MatrixType> plaket_time_left =
      calculate_plaket_schwinger_time_left(conf_plaket);
  std::vector<MatrixType> plaket_time_right =
      calculate_plaket_schwinger_time_right(conf_plaket);
  std::vector<MatrixType> wilson_loops;
  std::vector<MatrixType> wilson_loops_opposite;
  std::map<std::tuple<int, int, int>, double> result;
  std::map<int, double> schwinger_electric;
  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {
      wilson_loops = calculate_wilson_loops_schwinger(conf_wilson, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(conf_wilson, r, t);
      schwinger_electric = wilson_plaket_schwinger_l_dir_l_plaket(
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
std::map<int, double> wilson_plaket_schwinger_l_dir_tr_plaket(
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
flux_schwinger_electric_l_dir_tr_plaket(
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
std::map<int, double> wilson_plaket_schwinger_electric_tr_dir_l_plaket_even(
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
std::map<int, double> wilson_plaket_schwinger_electric_tr_dir_l_plaket_odd(
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
flux_schwinger_electric_tr_dir_l_plaket(
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
std::map<int, double> wilson_plaket_schwinger_electric_tr_dir_tr_plaket_even(
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
std::map<int, double> wilson_plaket_schwinger_electric_tr_dir_tr_plaket_odd(
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
flux_schwinger_electric_tr_dir_tr_plaket(
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
flux_schwinger_magnetic_l_dir_l_plaket(
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
flux_schwinger_magnetic_l_dir_tr_plaket(
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
flux_schwinger_magnetic_tr_dir_l_plaket(
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
flux_schwinger_magnetic_tr_dir_tr_plaket(
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

// for both electric and magnetic correlators
// only for case if R is even number
// d_outside - number of sites outside of wilson loop, is non-negative
// d_outside + R / 2 + 1 values, corresponding to sites along spatial direction
// of wilson loop, where index 0 corresponds to the values in the middle of a
// string
template <class DataPattern>
std::tuple<std::vector<double>, std::vector<double>>
wilson_plaket_correlator_longitudinal_even(
    const std::vector<double> &wilson_loop_tr,
    const std::vector<double> &plaket_time_tr,
    const std::vector<double> &plaket_space_tr, DataPattern &data_pattern,
    int d_outside, int T, int R) {
  int index_wilson;
  int correlator_size = d_outside + R / 2 + 1;
  std::vector<double> correlator_electric(correlator_size);
  std::vector<double> correlator_magnetic(correlator_size);
#pragma omp parallel for collapse(4)                                           \
    firstprivate(data_pattern, d_outside, R, T) private(index_wilson)          \
    reduction(vec_double_plus : correlator_electric)                           \
    reduction(vec_double_plus : correlator_magnetic)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index_wilson = data_pattern.get_index_site() * 3;
          data_pattern.move_forward(T / 2, 3);
          for (int mu = 0; mu < 3; mu++) {
            data_pattern.move_backward(d_outside, mu);
            for (int d = correlator_size - 1; d >= 0; d--) {
              correlator_electric[d] +=
                  wilson_loop_tr[index_wilson + mu] *
                  plaket_time_tr[data_pattern.get_index_site()];
              correlator_magnetic[d] +=
                  wilson_loop_tr[index_wilson + mu] *
                  plaket_space_tr[data_pattern.get_index_site()];
              data_pattern.move_forward(1, mu);
            }
            for (int d = 1; d < correlator_size; d++) {
              correlator_electric[d] +=
                  wilson_loop_tr[index_wilson + mu] *
                  plaket_time_tr[data_pattern.get_index_site()];
              correlator_magnetic[d] +=
                  wilson_loop_tr[index_wilson + mu] *
                  plaket_space_tr[data_pattern.get_index_site()];
              data_pattern.move_forward(1, mu);
            }
            data_pattern.move_backward(R + d_outside + 1, mu);
          }
        }
      }
    }
  }
  correlator_electric[0] /= data_pattern.get_lattice_size() * 3;
  for (int i = 1; i < correlator_electric.size(); i++) {
    correlator_electric[i] /= data_pattern.get_lattice_size() * 6;
  }
  correlator_magnetic[0] /= data_pattern.get_lattice_size() * 3;
  for (int i = 1; i < correlator_magnetic.size(); i++) {
    correlator_magnetic[i] /= data_pattern.get_lattice_size() * 6;
  }
  return {correlator_electric, correlator_magnetic};
}

// for both electric and magnetic correlators
// only for case if R is odd number
// d_outside - number of sites outside of wilson loop, is non-negative
// d_outside + (R + 1) / 2 values, corresponding to sites along spatial
// direction of wilson loop, where index 0 corresponds to the values closest to
// the middle of a string
template <class DataPattern>
std::tuple<std::vector<double>, std::vector<double>>
wilson_plaket_correlator_longitudinal_odd(
    const std::vector<double> &wilson_loop_tr,
    const std::vector<double> &plaket_time_tr,
    const std::vector<double> &plaket_space_tr, DataPattern &data_pattern,
    int d_outside, int T, int R) {
  int index_wilson;
  int correlator_size = d_outside + (R + 1) / 2;
  std::vector<double> correlator_electric(correlator_size);
  std::vector<double> correlator_magnetic(correlator_size);
#pragma omp parallel for collapse(4)                                           \
    firstprivate(data_pattern, d_outside, R, T) private(index_wilson)          \
    reduction(vec_double_plus : correlator_electric)                           \
    reduction(vec_double_plus : correlator_magnetic)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index_wilson = data_pattern.get_index_site() * 3;
          data_pattern.move_forward(T / 2, 3);
          for (int mu = 0; mu < 3; mu++) {
            data_pattern.move_backward(d_outside, mu);
            for (int d = correlator_size - 1; d >= 0; d--) {
              correlator_electric[d] +=
                  wilson_loop_tr[index_wilson + mu] *
                  plaket_time_tr[data_pattern.get_index_site()];
              correlator_magnetic[d] +=
                  wilson_loop_tr[index_wilson + mu] *
                  plaket_space_tr[data_pattern.get_index_site()];
              data_pattern.move_forward(1, mu);
            }
            for (int d = 0; d < correlator_size; d++) {
              correlator_electric[d] +=
                  wilson_loop_tr[index_wilson + mu] *
                  plaket_time_tr[data_pattern.get_index_site()];
              correlator_magnetic[d] +=
                  wilson_loop_tr[index_wilson + mu] *
                  plaket_space_tr[data_pattern.get_index_site()];
              data_pattern.move_forward(1, mu);
            }
            data_pattern.move_backward(R + d_outside + 1, mu);
          }
        }
      }
    }
  }
  for (int i = 0; i < correlator_electric.size(); i++) {
    correlator_electric[i] /= data_pattern.get_lattice_size() * 6;
  }
  for (int i = 0; i < correlator_magnetic.size(); i++) {
    correlator_magnetic[i] /= data_pattern.get_lattice_size() * 6;
  }
  return {correlator_electric, correlator_magnetic};
}

// for both electric and magnetic correlators
// only for case if R is even number
// d_outside - number of sites outside of wilson loop, is non-negative
// d_outside + 1 values, corresponding to sites along transversal
// direction in the center of wilson loop, where index 0 corresponds to the
// values in the middle of a string
template <class DataPattern>
std::tuple<std::vector<double>, std::vector<double>>
wilson_plaket_correlator_transversal_even(
    const std::vector<double> &wilson_loop_tr,
    const std::vector<double> &plaket_time_tr,
    const std::vector<double> &plaket_space_tr, DataPattern &data_pattern,
    int d_outside, int T, int R) {
  int index_wilson;
  int correlator_size = d_outside + 1;
  std::vector<double> correlator_electric(correlator_size);
  std::vector<double> correlator_magnetic(correlator_size);
#pragma omp parallel for collapse(4)                                           \
    firstprivate(data_pattern, d_outside, R, T) private(index_wilson)          \
    reduction(vec_double_plus : correlator_electric)                           \
    reduction(vec_double_plus : correlator_magnetic)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index_wilson = data_pattern.get_index_site() * 3;
          data_pattern.move_forward(T / 2, 3);
          for (int mu = 0; mu < 3; mu++) {
            data_pattern.move_forward(R / 2, mu);
            for (int nu = 0; nu < 3; nu++) {
              if (nu != mu) {
                data_pattern.move_backward(d_outside, nu);
                for (int d = correlator_size - 1; d >= 0; d--) {
                  correlator_electric[d] +=
                      wilson_loop_tr[index_wilson + mu] *
                      plaket_time_tr[data_pattern.get_index_site()];
                  correlator_magnetic[d] +=
                      wilson_loop_tr[index_wilson + mu] *
                      plaket_space_tr[data_pattern.get_index_site()];
                  data_pattern.move_forward(1, nu);
                }
                for (int d = 1; d < correlator_size; d++) {
                  correlator_electric[d] +=
                      wilson_loop_tr[index_wilson + mu] *
                      plaket_time_tr[data_pattern.get_index_site()];
                  correlator_magnetic[d] +=
                      wilson_loop_tr[index_wilson + mu] *
                      plaket_space_tr[data_pattern.get_index_site()];
                  data_pattern.move_forward(1, nu);
                }
                data_pattern.move_backward(d_outside + 1, nu);
              }
            }
            data_pattern.move_backward(R / 2, mu);
          }
        }
      }
    }
  }
  correlator_electric[0] /= data_pattern.get_lattice_size() * 6;
  for (int i = 1; i < correlator_electric.size(); i++) {
    correlator_electric[i] /= data_pattern.get_lattice_size() * 12;
  }
  correlator_magnetic[0] /= data_pattern.get_lattice_size() * 6;
  for (int i = 1; i < correlator_magnetic.size(); i++) {
    correlator_magnetic[i] /= data_pattern.get_lattice_size() * 12;
  }
  return {correlator_electric, correlator_magnetic};
}

// for both electric and magnetic correlators
// only for case if R is odd number
// d_outside - number of sites outside of wilson loop, is non-negative
// d_outside + 1 values, corresponding to sites along transversal
// direction in the center of wilson loop, where index 0 corresponds to the
// values in the middle of a string
template <class DataPattern>
std::tuple<std::vector<double>, std::vector<double>>
wilson_plaket_correlator_transversal_odd(
    const std::vector<double> &wilson_loop_tr,
    const std::vector<double> &plaket_time_tr,
    const std::vector<double> &plaket_space_tr, DataPattern &data_pattern,
    int d_outside, int T, int R) {
  int index_wilson;
  int correlator_size = d_outside + 1;
  std::vector<double> correlator_electric(correlator_size);
  std::vector<double> correlator_magnetic(correlator_size);
#pragma omp parallel for collapse(4)                                           \
    firstprivate(data_pattern, d_outside, R, T) private(index_wilson)          \
    reduction(vec_double_plus : correlator_electric)                           \
    reduction(vec_double_plus : correlator_magnetic)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index_wilson = data_pattern.get_index_site() * 3;
          data_pattern.move_forward(T / 2, 3);
          for (int mu = 0; mu < 3; mu++) {
            data_pattern.move_forward(R / 2, mu);
            for (int nu = 0; nu < 3; nu++) {
              if (nu != mu) {
                data_pattern.move_backward(d_outside, nu);
                for (int d = d_outside; d >= 0; d--) {
                  correlator_electric[d] +=
                      wilson_loop_tr[index_wilson + mu] *
                      plaket_time_tr[data_pattern.get_index_site()];
                  correlator_magnetic[d] +=
                      wilson_loop_tr[index_wilson + mu] *
                      plaket_space_tr[data_pattern.get_index_site()];
                  data_pattern.move_forward(1, mu);
                  correlator_electric[d] +=
                      wilson_loop_tr[index_wilson + mu] *
                      plaket_time_tr[data_pattern.get_index_site()];
                  correlator_magnetic[d] +=
                      wilson_loop_tr[index_wilson + mu] *
                      plaket_space_tr[data_pattern.get_index_site()];
                  data_pattern.move_backward(1, mu);
                  data_pattern.move_forward(1, nu);
                }
                for (int d = 1; d < correlator_size; d++) {
                  correlator_electric[d] +=
                      wilson_loop_tr[index_wilson + mu] *
                      plaket_time_tr[data_pattern.get_index_site()];
                  correlator_magnetic[d] +=
                      wilson_loop_tr[index_wilson + mu] *
                      plaket_space_tr[data_pattern.get_index_site()];
                  data_pattern.move_forward(1, mu);
                  correlator_electric[d] +=
                      wilson_loop_tr[index_wilson + mu] *
                      plaket_time_tr[data_pattern.get_index_site()];
                  correlator_magnetic[d] +=
                      wilson_loop_tr[index_wilson + mu] *
                      plaket_space_tr[data_pattern.get_index_site()];
                  data_pattern.move_backward(1, mu);
                  data_pattern.move_forward(1, nu);
                }
                data_pattern.move_backward(d_outside + 1, nu);
              }
            }
            data_pattern.move_backward(R / 2, mu);
          }
        }
      }
    }
  }
  correlator_electric[0] /= data_pattern.get_lattice_size() * 12;
  for (int i = 1; i < correlator_electric.size(); i++) {
    correlator_electric[i] /= data_pattern.get_lattice_size() * 24;
  }
  correlator_magnetic[0] /= data_pattern.get_lattice_size() * 12;
  for (int i = 1; i < correlator_magnetic.size(); i++) {
    correlator_magnetic[i] /= data_pattern.get_lattice_size() * 24;
  }
  return {correlator_electric, correlator_magnetic};
}

// template <class DataPattern, class MatrixType>
// void wilson_plaket_correlator_all(
//     std::map<std::tuple<int, int, int>, double>
//         &flux_tube_electric_longitudinal,
//     std::map<std::tuple<int, int, int>, double>
//     &flux_tube_electric_transversal, std::map<std::tuple<int, int, int>,
//     double>
//         &flux_tube_magnetic_longitudinal,
//     std::map<std::tuple<int, int, int>, double>
//     &flux_tube_magnetic_transversal, const std::vector<double>
//     &plaket_time_tr, const std::vector<double> &plaket_space_tr, const
//     Data::LatticeData<DataPattern, MatrixType> &conf, int T_min, int T_max,
//     int R_min, int R_max, int d_outside) {
//   DataPattern data_pattern(conf.lat_dim);
//   std::vector<MatrixType> time_lines;
//   std::array<std::vector<MatrixType>, 3> space_lines;
//   std::vector<double> wilson_loop_tr;
//   std::tuple<std::vector<double>, std::vector<double>> correlator;
//   time_lines = wilson_lines_indexed(conf, T_min, 3);
//   for (int t = T_min; t <= T_max; t += 2) {
//     for (int mu = 0; mu < 3; mu++) {
//       space_lines[mu] = wilson_lines_indexed(conf, R_min, mu);
//     }
//     for (int r = R_min; r <= R_max; r++) {
//       wilson_loop_tr = calculate_wilson_loop_time_tr(time_lines, space_lines,
//                                                      data_pattern, t, r);
//       if (r % 2 == 0) {
//         correlator = wilson_plaket_correlator_longitudinal_even(
//             wilson_loop_tr, plaket_time_tr, plaket_space_tr, data_pattern,
//             d_outside, t, r);
//         for (int d = 0; d < std::get<0>(correlator).size(); d++) {
//           flux_tube_electric_longitudinal[{t, r, d}] =
//               std::get<0>(correlator)[d];
//         }
//         for (int d = 0; d < std::get<1>(correlator).size(); d++) {
//           flux_tube_magnetic_longitudinal[{t, r, d}] =
//               std::get<1>(correlator)[d];
//         }
//         correlator = wilson_plaket_correlator_transversal_even(
//             wilson_loop_tr, plaket_time_tr, plaket_space_tr, data_pattern,
//             d_outside, t, r);
//         for (int d = 0; d < std::get<0>(correlator).size(); d++) {
//           flux_tube_electric_transversal[{t, r, d}] =
//               std::get<0>(correlator)[d];
//         }
//         for (int d = 0; d < std::get<1>(correlator).size(); d++) {
//           flux_tube_magnetic_transversal[{t, r, d}] =
//               std::get<1>(correlator)[d];
//         }
//       } else {
//         correlator = wilson_plaket_correlator_longitudinal_odd(
//             wilson_loop_tr, plaket_time_tr, plaket_space_tr, data_pattern,
//             d_outside, t, r);
//         for (int d = 0; d < std::get<0>(correlator).size(); d++) {
//           flux_tube_electric_longitudinal[{t, r, d}] =
//               std::get<0>(correlator)[d];
//         }
//         for (int d = 0; d < std::get<1>(correlator).size(); d++) {
//           flux_tube_magnetic_longitudinal[{t, r, d}] =
//               std::get<1>(correlator)[d];
//         }
//         correlator = wilson_plaket_correlator_transversal_odd(
//             wilson_loop_tr, plaket_time_tr, plaket_space_tr, data_pattern,
//             d_outside, t, r);
//         for (int d = 0; d < std::get<0>(correlator).size(); d++) {
//           flux_tube_electric_transversal[{t, r, d}] =
//               std::get<0>(correlator)[d];
//         }
//         for (int d = 0; d < std::get<1>(correlator).size(); d++) {
//           flux_tube_magnetic_transversal[{t, r, d}] =
//               std::get<1>(correlator)[d];
//         }
//       }
//       for (int mu = 0; mu < 3; mu++) {
//         wilson_lines_prolong(conf, space_lines[mu], r, mu);
//       }
//     }
//     wilson_lines_prolong(conf, time_lines, t, 3);
//     wilson_lines_prolong(conf, time_lines, t + 1, 3);
//   }
// }

template <class DataPattern, class MatrixType>
void wilson_plaket_correlator_all(
    std::map<std::tuple<int, int>, double> &wilson_loops,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_electric_longitudinal,
    std::map<std::tuple<int, int, int>, double> &flux_tube_electric_transversal,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_magnetic_longitudinal,
    std::map<std::tuple<int, int, int>, double> &flux_tube_magnetic_transversal,
    const std::vector<double> &plaket_time_tr,
    const std::vector<double> &plaket_space_tr,
    const Data::LatticeData<DataPattern, MatrixType> &conf, int T_min,
    int T_max, int R_min, int R_max, int d_outside) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<std::vector<MatrixType>> time_lines((T_max - T_min) / 2 + 1);
  std::array<std::vector<MatrixType>, 3> space_lines;
  std::vector<double> wilson_loop_tr;
  double wilson_aver;
  std::tuple<std::vector<double>, std::vector<double>> correlator;
  time_lines[0] = wilson_lines_indexed(conf, T_min, 3);
  for (int T = T_min + 2; T <= T_max; T += 2) {
    time_lines[(T - T_min) / 2] = time_lines[(T - T_min - 2) / 2];
    wilson_lines_prolong(conf, time_lines[(T - T_min) / 2], T - 2, 3);
    wilson_lines_prolong(conf, time_lines[(T - T_min) / 2], T - 1, 3);
  }
  for (int mu = 0; mu < 3; mu++) {
    space_lines[mu] = wilson_lines_indexed(conf, R_min, mu);
  }
  for (int R = R_min; R <= R_max; R++) {
    for (int T = T_min; T <= T_max; T += 2) {
      wilson_loop_tr = calculate_wilson_loop_time_tr(
          time_lines[(T - T_min) / 2], space_lines, data_pattern, T, R);
      wilson_aver = 0;
#pragma omp parallel for reduction(+ : wilson_aver)
      for (int i = 0; i < wilson_loop_tr.size(); i++) {
        wilson_aver += wilson_loop_tr[i];
      }
      wilson_loops[{T, R}] = wilson_aver / wilson_loop_tr.size();
      if (R % 2 == 0) {
        correlator = wilson_plaket_correlator_longitudinal_even(
            wilson_loop_tr, plaket_time_tr, plaket_space_tr, data_pattern,
            d_outside, T, R);
        for (int d = 0; d < std::get<0>(correlator).size(); d++) {
          flux_tube_electric_longitudinal[{T, R, d}] =
              std::get<0>(correlator)[d];
        }
        for (int d = 0; d < std::get<1>(correlator).size(); d++) {
          flux_tube_magnetic_longitudinal[{T, R, d}] =
              std::get<1>(correlator)[d];
        }
        correlator = wilson_plaket_correlator_transversal_even(
            wilson_loop_tr, plaket_time_tr, plaket_space_tr, data_pattern,
            d_outside, T, R);
        for (int d = 0; d < std::get<0>(correlator).size(); d++) {
          flux_tube_electric_transversal[{T, R, d}] =
              std::get<0>(correlator)[d];
        }
        for (int d = 0; d < std::get<1>(correlator).size(); d++) {
          flux_tube_magnetic_transversal[{T, R, d}] =
              std::get<1>(correlator)[d];
        }
      } else {
        correlator = wilson_plaket_correlator_longitudinal_odd(
            wilson_loop_tr, plaket_time_tr, plaket_space_tr, data_pattern,
            d_outside, T, R);
        for (int d = 0; d < std::get<0>(correlator).size(); d++) {
          flux_tube_electric_longitudinal[{T, R, d}] =
              std::get<0>(correlator)[d];
        }
        for (int d = 0; d < std::get<1>(correlator).size(); d++) {
          flux_tube_magnetic_longitudinal[{T, R, d}] =
              std::get<1>(correlator)[d];
        }
        correlator = wilson_plaket_correlator_transversal_odd(
            wilson_loop_tr, plaket_time_tr, plaket_space_tr, data_pattern,
            d_outside, T, R);
        for (int d = 0; d < std::get<0>(correlator).size(); d++) {
          flux_tube_electric_transversal[{T, R, d}] =
              std::get<0>(correlator)[d];
        }
        for (int d = 0; d < std::get<1>(correlator).size(); d++) {
          flux_tube_magnetic_transversal[{T, R, d}] =
              std::get<1>(correlator)[d];
        }
      }
    }
    for (int mu = 0; mu < 3; mu++) {
      wilson_lines_prolong(conf, space_lines[mu], R, mu);
    }
  }
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