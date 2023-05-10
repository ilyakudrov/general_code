#define data_size 4 * x_size *y_size *z_size *t_size
#define PLACE3_DIR                                                             \
  (t) * 3 * x_size *y_size *z_size + (z)*3 * x_size *y_size + (y)*3 * x_size + \
      (x)*3 + dir
#define PLACE3_LINK_NODIR                                                      \
  (link.coordinate[3]) * 3 * x_size *y_size *z_size +                          \
      (link.coordinate[2]) * 3 * x_size *y_size +                              \
      (link.coordinate[1]) * 3 * x_size + (link.coordinate[0]) * 3
#define PLACE1_NODIR                                                           \
  (t) * x_size *y_size *z_size + (z)*x_size *y_size + (y)*x_size + (x)
#define PLACE_PLAKET_TIME                                                      \
  (link.coordinate[3]) * 3 * x_size *y_size *z_size +                          \
      (link.coordinate[2]) * 3 * x_size *y_size +                              \
      (link.coordinate[1]) * 3 * x_size + (link.coordinate[0]) * 3 +           \
      link.direction
#define PLACE_PLAKET_SPACE                                                     \
  (link.coordinate[3]) * 3 * x_size *y_size *z_size +                          \
      (link.coordinate[2]) * 3 * x_size *y_size +                              \
      (link.coordinate[1]) * 3 * x_size + (link.coordinate[0]) * 3
#define PLACE1_LINK_NODIR                                                      \
  (link.coordinate[3]) * x_size *y_size *z_size +                              \
      (link.coordinate[2]) * x_size *y_size + (link.coordinate[1]) * x_size +  \
      (link.coordinate[0])

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

#define ITER_START_ZYX                                                         \
  for (int z = 0; z < z_size; z++) {                                           \
    for (int y = 0; y < y_size; y++) {                                         \
      for (int x = 0; x < x_size; x++) {                                       \
        link.go(x, y, z, t);                                                   \
        link.update(0);                                                        \
        link.update(1);                                                        \
        link.update(2);

#define ITER_END_3                                                             \
  }                                                                            \
  }                                                                            \
  }

#include "../include/flux_tube.h"
#include "../include/basic_observables.h"
#include "../include/link.h"
#include <omp.h>

#include <algorithm>

template <class T>
std::vector<T> calculate_plaket_time_left_down(const std::vector<T> array) {
  std::vector<T> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
#pragma omp parallel for collapse(4) private(link)
  SPACE_ITER_START;
  for (int dir = 0; dir < 3; dir++) {
    link.move_dir(dir);
    vec[link.place / 4 * 3 + dir] = link.plaket_left_down(array, 3);
  }
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<T> calculate_plaket_time_left_up(const std::vector<T> array) {
  std::vector<T> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
#pragma omp parallel for collapse(4) private(link)
  SPACE_ITER_START;
  for (int dir = 0; dir < 3; dir++) {
    link.move_dir(dir);
    vec[link.place / 4 * 3 + dir] = link.plaket_left_up(array, 3);
  }
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<T> calculate_plaket_time_right_down(const std::vector<T> array) {
  std::vector<T> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
#pragma omp parallel for collapse(4) private(link)
  SPACE_ITER_START;
  for (int dir = 0; dir < 3; dir++) {
    link.move_dir(dir);
    vec[link.place / 4 * 3 + dir] = link.plaket_right_down(array, 3);
  }
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<T> calculate_plaket_time_right_up(const std::vector<T> array) {
  std::vector<T> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
#pragma omp parallel for collapse(4) private(link)
  SPACE_ITER_START;
  for (int dir = 0; dir < 3; dir++) {
    link.move_dir(dir);
    vec[link.place / 4 * 3 + dir] = link.plaket_right_up(array, 3);
  }
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<T>
calculate_plaket_schwinger_time_left(const std::vector<T> array) {
  std::vector<T> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
#pragma omp parallel for collapse(4) private(link)
  SPACE_ITER_START;
  for (int dir = 0; dir < 3; dir++) {
    link.move_dir(dir);
    vec[link.place / 4 * 3 + dir] =
        (link.plaket_left_down(array, 3) + link.plaket_left_up(array, 3)) * 0.5;
  }
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<T>
calculate_plaket_schwinger_time_right(const std::vector<T> array) {
  std::vector<T> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
#pragma omp parallel for collapse(4) private(link)
  SPACE_ITER_START;
  for (int dir = 0; dir < 3; dir++) {
    link.move_dir(dir);
    vec[link.place / 4 * 3 + dir] =
        (link.plaket_right_down(array, 3) + link.plaket_right_up(array, 3)) *
        0.5;
  }
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<T> calculate_plaket_schwinger_time_tr(const std::vector<T> array) {
  std::vector<T> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
#pragma omp parallel for collapse(4) private(link)
  SPACE_ITER_START;
  for (int dir = 0; dir < 3; dir++) {
    for (int mu = 0; mu < 3; mu++) {
      if (dir != mu) {
        link.move_dir(mu);
        vec[link.place / 4 * 3 + dir] =
            vec[link.place / 4 * 3 + dir] +
            (link.plaket_right_down(array, 3).conj() +
             link.plaket_right_up(array, 3).conj() +
             link.plaket_left_down(array, 3) + link.plaket_left_up(array, 3)) *
                0.125;
      }
    }
  }
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<T> calculate_plaket_schwinger_space_l(const std::vector<T> &array) {
  std::vector<T> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
#pragma omp parallel for collapse(4) private(link)
  SPACE_ITER_START;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i != j) {
        link.move_dir(j);
        for (int k = 0; k < 3; k++) {
          if (k != i && k != j) {
            vec[link.place / 4 * 3 + i] = vec[link.place / 4 * 3 + i] +
                                          (link.plaket_right_down(array, k) +
                                           link.plaket_right_up(array, k) +
                                           link.plaket_left_down(array, k) +
                                           link.plaket_left_up(array, k)) *
                                              0.25;
          }
        }
      }
    }
  }
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<T>
calculate_plaket_schwinger_space_tr(const std::vector<T> &array) {
  std::vector<T> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
#pragma omp parallel for collapse(4) private(link)
  SPACE_ITER_START;
  for (int dir = 0; dir < 3; dir++) {
    for (int i = 0; i < 3; i++) {
      if (i != dir) {
        for (int j = 0; j < 3; j++) {
          if (j != i) {
            link.move_dir(j);
            for (int k = 0; k < 3; k++) {
              if (k != i && k != j) {
                vec[link.place / 4 * 3 + dir] =
                    vec[link.place / 4 * 3 + dir] +
                    (link.plaket_right_down(array, k) +
                     link.plaket_right_up(array, k) +
                     link.plaket_left_down(array, k) +
                     link.plaket_left_up(array, k)) *
                        0.125;
              }
            }
          }
        }
      }
    }
  }
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<T> calculate_polyakov_loop(const std::vector<T> &array) {
  std::vector<T> vec((data_size) / 4);
  link1 link(x_size, y_size, z_size, t_size);
  link.move_dir(3);
  SPACE_ITER_START;
  vec[PLACE1_NODIR] = link.polyakov_loop(array);
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<T> calculate_wilson_loop(const std::vector<T> &array, int r,
                                     int time) {
  std::vector<T> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
  for (int dir = 0; dir < 3; dir++) {
    link.move_dir(dir);
    SPACE_ITER_START;
    vec[PLACE3_DIR] = link.wilson_loop(array, r, time);
    SPACE_ITER_END;
  }
  return vec;
}

template <class T>
std::vector<T> calculate_schwinger_lines_short(const std::vector<T> &array,
                                               int d) {
  std::vector<T> vec(data_size * 3 / 4);
  link1 link(x_size, y_size, z_size, t_size);
#pragma omp parallel for collapse(4) private(link)
  SPACE_ITER_START
  for (int mu = 0; mu < 3; mu++) {
    link.move_dir(mu);
    vec[link.place * 3 / 4 + mu] = link.wilson_line(array, d);
  }
  SPACE_ITER_END
  return vec;
}

template <class T>
std::vector<T> calculate_wilson_loops_schwinger(const std::vector<T> &array,
                                                int r, int time) {
  std::vector<T> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
#pragma omp parallel for collapse(4) private(link)
  SPACE_ITER_START;
  for (int mu = 0; mu < 3; mu++) {
    link.move_dir(mu);
    vec[link.place / 4 * 3 + mu] = link.wilson_loop_schwinger(array, r, time);
  }
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<T>
calculate_wilson_loops_schwinger_opposite(const std::vector<T> &array, int r,
                                          int time) {
  std::vector<T> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
#pragma omp parallel for collapse(4) private(link)
  SPACE_ITER_START;
  for (int mu = 0; mu < 3; mu++) {
    link.move_dir(mu);
    vec[link.place / 4 * 3 + mu] =
        link.wilson_loop_schwinger_opposite(array, r, time);
  }
  SPACE_ITER_END;
  return vec;
}

#pragma omp declare reduction(vec_double_plus : std::vector<double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

template <class T>
std::map<int, double> wilson_plaket_schwinger_longitudinal_l(
    const std::vector<T> &wilson_loops,
    const std::vector<T> &wilson_loops_opposite,
    const std::vector<T> &plaket_left, const std::vector<T> &plaket_right,
    const std::vector<std::vector<T>> &schwinger_lines, int d_ouside, int time,
    int r) {
  link1 link(x_size, y_size, z_size, t_size);

  std::vector<double> correlator(2 * d_ouside + r);
  T W;
  T A;
  T S;
  int d;
  int place;

#pragma omp parallel for collapse(4) private(W, A, S, link, d, place)          \
    firstprivate(d_ouside) reduction(vec_double_plus                           \
                                     : correlator)
  SPACE_ITER_START

  for (int dir = 0; dir < 3; dir++) {
    link.go_update(x, y, z, t);
    W = wilson_loops[link.place / 4 * 3 + dir];
    d = -d_ouside;
    link.move(dir, 1 - d_ouside);
    while (d < -1) {
      place = link.place / 4 * 3;
      S = schwinger_lines[abs(d) - 2][place + dir];
      A = W ^ S;
      A = A ^ plaket_right[place + dir];
      correlator[d + d_ouside] += A.multiply_tr(S);
      link.move(dir, 1);
      d++;
    }
    correlator[d + d_ouside] +=
        W.multiply_conj_tr(plaket_right[link.place / 4 * 3 + dir]);
    d++;
    correlator[d + d_ouside] +=
        W.multiply_tr(plaket_left[link.place / 4 * 3 + dir]);
    d++;
    while (d < (r + 1) / 2) {
      S = schwinger_lines[d - 1][link.place / 4 * 3 + dir];
      A = W * S;
      link.move(dir, d);
      A = A * plaket_left[link.place / 4 * 3 + dir];
      correlator[d + d_ouside] += A.multiply_conj_tr(S);
      link.move(dir, -d);
      d++;
    }
    link.move(dir, r);
    W = wilson_loops_opposite[link.place / 4 * 3 + dir];
    link.move(dir, -(r / 2 - 1));
    while (d < r - 1) {
      S = schwinger_lines[r - d - 2][link.place / 4 * 3 + dir];
      A = W ^ S;
      A = A * plaket_right[link.place / 4 * 3 + dir];
      correlator[d + d_ouside] += A.multiply_tr(S);
      link.move(dir, 1);
      d++;
    }
    correlator[d + d_ouside] +=
        W.multiply_tr(plaket_right[link.place / 4 * 3 + dir]);
    d++;
    correlator[d + d_ouside] +=
        W.multiply_conj_tr(plaket_left[link.place / 4 * 3 + dir]);
    d++;
    while (d < r + d_ouside) {
      S = schwinger_lines[d - r - 1][link.place / 4 * 3 + dir];
      A = W * S;
      link.move(dir, d - r);
      A = A ^ plaket_left[link.place / 4 * 3 + dir];
      correlator[d + d_ouside] += A.multiply_conj_tr(S);
      link.move(dir, -(d - r));
      d++;
    }
  }
  SPACE_ITER_END

  std::map<int, double> result;
  for (int i = 0; i < correlator.size(); i++) {
    result[i - d_ouside - r / 2] =
        correlator[i] / (x_size * y_size * z_size * t_size * 3);
  }
  return result;
}

template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_longitudinal_l(
    const std::vector<T> &array,
    std::vector<std::vector<T>> &schwinger_lines_short, int T_min, int T_max,
    int R_min, int R_max, int d_ouside) {

  std::vector<T> plaket_time_left = calculate_plaket_schwinger_time_left(array);

  std::vector<T> plaket_time_right =
      calculate_plaket_schwinger_time_right(array);

  std::vector<T> wilson_loops;
  std::vector<T> wilson_loops_opposite;
  std::map<std::tuple<int, int, int>, double> result;
  std::map<int, double> schwinger_electric;

  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {

      wilson_loops = calculate_wilson_loops_schwinger(array, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(array, r, t);

      schwinger_electric = wilson_plaket_schwinger_longitudinal_l(
          wilson_loops, wilson_loops_opposite, plaket_time_left,
          plaket_time_right, schwinger_lines_short, d_ouside, t, r);
      for (auto it = schwinger_electric.begin(); it != schwinger_electric.end();
           ++it) {
        result[std::tuple<int, int, double>(t, r, it->first)] = it->second;
      }
    }
  }

  return result;
}

template <class T>
std::map<int, double> wilson_plaket_schwinger_longitudinal_tr(
    std::vector<T> &wilson_loops, std::vector<T> &wilson_loops_opposite,
    const std::vector<T> &plaket,
    const std::vector<std::vector<T>> &schwinger_lines, int d_ouside, int time,
    int r) {
  link1 link(x_size, y_size, z_size, t_size);

  std::vector<double> correlator(2 * d_ouside + r + 1);
  T W;
  T A;
  T S;
  int d;
  int place;

#pragma omp parallel for collapse(4) private(W, A, S, link, d, place)          \
    firstprivate(d_ouside) reduction(vec_double_plus                           \
                                     : correlator)
  SPACE_ITER_START

  for (int dir = 0; dir < 3; dir++) {
    link.go_update(x, y, z, t);
    W = wilson_loops[link.place / 4 * 3 + dir];
    d = -d_ouside;
    link.move(dir, -d_ouside);
    while (d < 0) {
      place = link.place / 4 * 3;
      S = schwinger_lines[abs(d) - 1][place + dir];
      A = W ^ S;
      A = A * plaket[place + dir];
      correlator[d + d_ouside] += A.multiply_tr(S);
      link.move(dir, 1);
      d++;
    }
    correlator[d + d_ouside] += W.multiply_tr(plaket[link.place / 4 * 3 + dir]);
    d++;
    while (d <= r / 2) {
      S = schwinger_lines[d - 1][link.place / 4 * 3 + dir];
      A = W * S;
      link.move(dir, d);
      A = A * plaket[link.place / 4 * 3 + dir];
      correlator[d + d_ouside] += A.multiply_conj_tr(S);
      link.move(dir, -d);
      d++;
    }
    link.move(dir, r);
    W = wilson_loops_opposite[link.place / 4 * 3 + dir];
    link.move(dir, 1 - (r + 1) / 2);
    while (d < r) {
      S = schwinger_lines[r - d - 1][link.place / 4 * 3 + dir];
      A = W ^ S;
      A = A ^ plaket[link.place / 4 * 3 + dir];
      correlator[d + d_ouside] += A.multiply_tr(S);
      link.move(dir, 1);
      d++;
    }
    correlator[d + d_ouside] +=
        W.multiply_conj_tr(plaket[link.place / 4 * 3 + dir]);
    d++;
    while (d <= r + d_ouside) {
      S = schwinger_lines[d - r - 1][link.place / 4 * 3 + dir];
      A = W * S;
      link.move(dir, d - r);
      A = A ^ plaket[link.place / 4 * 3 + dir];
      correlator[d + d_ouside] += A.multiply_conj_tr(S);
      link.move(dir, -(d - r));
      d++;
    }
  }
  SPACE_ITER_END

  std::map<int, double> result;
  for (int i = 0; i < correlator.size(); i++) {
    result[i - d_ouside - r / 2] =
        correlator[i] / (x_size * y_size * z_size * t_size * 3);
  }
  return result;
}

template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_longitudinal_tr(
    std::vector<T> &array, std::vector<std::vector<T>> &schwinger_lines_short,
    int T_min, int T_max, int R_min, int R_max, int d_ouside) {

  std::vector<T> plaket_time = calculate_plaket_schwinger_time_tr(array);

  std::vector<T> wilson_loops;
  std::vector<T> wilson_loops_opposite;
  std::map<std::tuple<int, int, int>, double> result;
  std::map<int, double> schwinger_electric;

  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {

      wilson_loops = calculate_wilson_loops_schwinger(array, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(array, r, t);

      schwinger_electric = wilson_plaket_schwinger_longitudinal_tr(
          wilson_loops, wilson_loops_opposite, plaket_time,
          schwinger_lines_short, d_ouside, t, r);
      for (auto it = schwinger_electric.begin(); it != schwinger_electric.end();
           ++it) {
        result[std::tuple<int, int, double>(t, r, it->first)] = it->second;
      }
    }
  }

  return result;
}

template <class T>
std::map<int, double> wilson_plaket_schwinger_electric_transversal_l_even(
    std::vector<T> &wilson_loops, std::vector<T> &wilson_loops_opposite,
    const std::vector<T> &plaket_left, const std::vector<T> &plaket_right,
    const std::vector<std::vector<T>> &schwinger_lines, int d_max, int time,
    int r) {
  link1 link(x_size, y_size, z_size, t_size);

  std::vector<double> correlator(2 * d_max + 1);
  T W;
  T A1;
  T A2;
  T S1;
  T S2;
  int d;
  int place;

#pragma omp parallel for collapse(4) private(W, A1, A2, S1, S2, link, d,       \
                                             place) firstprivate(d_max)        \
    reduction(vec_double_plus                                                  \
              : correlator)
  SPACE_ITER_START

  for (int dir = 0; dir < 3; dir++) {
    link.go_update(x, y, z, t);
    W = wilson_loops[link.place / 4 * 3 + dir];
    S1 = schwinger_lines[r / 2 - 2][link.place / 4 * 3 + dir];
    A1 = W * S1;
    link.move(dir, r / 2 - 1);
    for (int mu = 0; mu < 3; mu++) {
      if (dir != mu) {
        link.move(mu, -d_max);
        for (int d = -d_max; d < 0; d++) {
          S2 = schwinger_lines[abs(d) - 1][link.place / 4 * 3 + mu];
          A2 = A1 ^ S2;
          A2 = A2 * plaket_left[link.place / 4 * 3 + dir];
          A2 = A2 * S2;
          correlator[d + d_max] += A2.multiply_conj_tr(S1);
          link.move(mu, 1);
        }
        A2 = A1 * plaket_left[link.place / 4 * 3 + dir];
        correlator[d_max] += A2.multiply_conj_tr(S1);
        for (int d = 1; d <= d_max; d++) {
          S2 = schwinger_lines[d - 1][link.place / 4 * 3 + mu];
          link.move(mu, d);
          A2 = A1 * S2;
          A2 = A2 * plaket_left[link.place / 4 * 3 + dir];
          A2 = A2 ^ S2;
          correlator[d + d_max] += A2.multiply_conj_tr(S1);
          link.move(mu, -d);
        }
      }
    }
    link.go_update(x, y, z, t);
    link.move(dir, r);
    W = wilson_loops_opposite[link.place / 4 * 3 + dir];
    link.move(dir, 1 - r / 2);
    S1 = schwinger_lines[r / 2 - 2][link.place / 4 * 3 + dir].conj();
    A1 = W * S1;
    for (int mu = 0; mu < 3; mu++) {
      if (dir != mu) {
        link.move(mu, -d_max);
        for (int d = -d_max; d < 0; d++) {
          S2 = schwinger_lines[abs(d) - 1][link.place / 4 * 3 + mu];
          A2 = A1 ^ S2;
          A2 = A2 ^ plaket_right[link.place / 4 * 3 + dir];
          A2 = A2 * S2;
          correlator[d + d_max] += A2.multiply_conj_tr(S1);
          link.move(mu, 1);
        }
        A2 = A1 * plaket_right[link.place / 4 * 3 + dir];
        correlator[d_max] += A2.multiply_conj_tr(S1);
        for (int d = 1; d <= d_max; d++) {
          S2 = schwinger_lines[d - 1][link.place / 4 * 3 + mu];
          link.move(mu, d);
          A2 = A1 * S2;
          A2 = A2 ^ plaket_right[link.place / 4 * 3 + dir];
          A2 = A2 ^ S2;
          correlator[d + d_max] += A2.multiply_conj_tr(S1);
          link.move(mu, -d);
        }
      }
    }
  }
  SPACE_ITER_END

  std::map<int, double> result;
  for (int i = 0; i < correlator.size(); i++) {
    result[i - d_max] =
        correlator[i] / (x_size * y_size * z_size * t_size * 12);
  }
  return result;
}

template <class T>
std::map<int, double> wilson_plaket_schwinger_electric_transversal_l_odd(
    std::vector<T> &wilson_loops, const std::vector<T> &plaket_left,
    const std::vector<T> &plaket_right,
    const std::vector<std::vector<T>> &schwinger_lines, int d_max, int time,
    int r) {
  link1 link(x_size, y_size, z_size, t_size);

  std::vector<double> correlator(2 * d_max + 1);
  T W;
  T A1;
  T A2;
  T S1;
  T S2;
  int d;
  int place;

#pragma omp parallel for collapse(4) private(W, A1, A2, S1, S2, link, d,       \
                                             place) firstprivate(d_max)        \
    reduction(vec_double_plus                                                  \
              : correlator)
  SPACE_ITER_START

  for (int dir = 0; dir < 3; dir++) {
    link.go_update(x, y, z, t);
    W = wilson_loops[link.place / 4 * 3 + dir];
    S1 = schwinger_lines[r / 2 - 2][link.place / 4 * 3 + dir];
    A1 = W * S1;
    link.move(dir, r / 2);
    for (int mu = 0; mu < 3; mu++) {
      if (dir != mu) {
        link.move(mu, -d_max);
        for (int d = -d_max; d < 0; d++) {
          S2 = schwinger_lines[abs(d) - 1][link.place / 4 * 3 + mu];
          A2 = A1 ^ S2;
          A2 = A2 * plaket_left[link.place / 4 * 3 + dir];
          A2 = A2 * S2;
          correlator[d + d_max] += A2.multiply_conj_tr(S1);
          link.move(mu, 1);
        }
        A2 = A1 * plaket_left[link.place / 4 * 3 + dir];
        correlator[d_max] += A2.multiply_conj_tr(S1);
        for (int d = 1; d <= d_max; d++) {
          S2 = schwinger_lines[d - 1][link.place / 4 * 3 + mu];
          link.move(mu, d);
          A2 = A1 * S2;
          A2 = A2 * plaket_left[link.place / 4 * 3 + dir];
          A2 = A2 ^ S2;
          correlator[d + d_max] += A2.multiply_conj_tr(S1);
          link.move(mu, -d);
        }
      }
    }
  }
  SPACE_ITER_END

  std::map<int, double> result;
  for (int i = 0; i < correlator.size(); i++) {
    result[i - d_max] = correlator[i] / (x_size * y_size * z_size * t_size * 6);
  }
  return result;
}

template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_transversal_l(
    std::vector<T> &array, std::vector<std::vector<T>> &schwinger_lines_short,
    int T_min, int T_max, int R_min, int R_max, int d_max) {

  std::vector<T> plaket_time_left = calculate_plaket_schwinger_time_left(array);

  std::vector<T> plaket_time_right =
      calculate_plaket_schwinger_time_right(array);

  std::vector<T> wilson_loops;
  std::vector<T> wilson_loops_opposite;
  std::map<std::tuple<int, int, int>, double> result;
  std::map<int, double> schwinger_electric;

  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {

      wilson_loops = calculate_wilson_loops_schwinger(array, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(array, r, t);

      if (r % 2 == 0) {
        schwinger_electric =
            wilson_plaket_schwinger_electric_transversal_l_even(
                wilson_loops, wilson_loops_opposite, plaket_time_left,
                plaket_time_right, schwinger_lines_short, d_max, t, r);
      } else {
        schwinger_electric = wilson_plaket_schwinger_electric_transversal_l_odd(
            wilson_loops, plaket_time_left, plaket_time_right,
            schwinger_lines_short, d_max, t, r);
      }
      for (auto it = schwinger_electric.begin(); it != schwinger_electric.end();
           ++it) {
        result[std::tuple<int, int, double>(t, r, it->first)] = it->second;
      }
    }
  }

  return result;
}

template <class T>
std::map<int, double> wilson_plaket_schwinger_electric_transversal_tr_even(
    std::vector<T> &wilson_loops, const std::vector<T> &plaket,
    const std::vector<std::vector<T>> &schwinger_lines, int d_max, int time,
    int r) {
  link1 link(x_size, y_size, z_size, t_size);

  std::vector<double> correlator(2 * d_max + 1);
  T W;
  T A1;
  T A2;
  T S1;
  T S2;
  int d;
  int place;

#pragma omp parallel for collapse(4) private(W, A1, A2, S1, S2, link, d,       \
                                             place) firstprivate(d_max)        \
    reduction(vec_double_plus                                                  \
              : correlator)
  SPACE_ITER_START

  for (int dir = 0; dir < 3; dir++) {
    link.go_update(x, y, z, t);
    W = wilson_loops[link.place / 4 * 3 + dir];
    S1 = schwinger_lines[r / 2 - 2][link.place / 4 * 3 + dir];
    A1 = W * S1;
    link.move(dir, r / 2 - 1);
    for (int mu = 0; mu < 3; mu++) {
      if (dir != mu) {
        link.move(mu, -d_max);
        for (int d = -d_max; d < 0; d++) {
          S2 = schwinger_lines[abs(d) - 1][link.place / 4 * 3 + mu];
          A2 = A1 ^ S2;
          A2 = A2 * plaket[link.place / 4 * 3 + dir];
          A2 = A2 * S2;
          correlator[d + d_max] += A2.multiply_conj_tr(S1);
          link.move(mu, 1);
        }
        A2 = A1 * plaket[link.place / 4 * 3 + dir];
        correlator[d_max] += A2.multiply_conj_tr(S1);
        for (int d = 1; d <= d_max; d++) {
          S2 = schwinger_lines[d - 1][link.place / 4 * 3 + mu];
          link.move(mu, d);
          A2 = A1 * S2;
          A2 = A2 * plaket[link.place / 4 * 3 + dir];
          A2 = A2 ^ S2;
          correlator[d + d_max] += A2.multiply_conj_tr(S1);
          link.move(mu, -d);
        }
      }
    }
  }
  SPACE_ITER_END

  std::map<int, double> result;
  for (int i = 0; i < correlator.size(); i++) {
    result[i - d_max] = correlator[i] / (x_size * y_size * z_size * t_size * 6);
  }
  return result;
}

template <class T>
std::map<int, double> wilson_plaket_schwinger_electric_transversal_tr_odd(
    std::vector<T> &wilson_loops, std::vector<T> &wilson_loops_opposite,
    const std::vector<T> &plaket,
    const std::vector<std::vector<T>> &schwinger_lines, int d_max, int time,
    int r) {
  link1 link(x_size, y_size, z_size, t_size);

  std::vector<double> correlator(2 * d_max + 1);
  T W;
  T A1;
  T A2;
  T S1;
  T S2;
  int d;
  int place;

#pragma omp parallel for collapse(4) private(W, A1, A2, S1, S2, link, d,       \
                                             place) firstprivate(d_max)        \
    reduction(vec_double_plus                                                  \
              : correlator)
  SPACE_ITER_START

  for (int dir = 0; dir < 3; dir++) {
    link.go_update(x, y, z, t);
    W = wilson_loops[link.place / 4 * 3 + dir];
    S1 = schwinger_lines[r / 2 - 2][link.place / 4 * 3 + dir];
    A1 = W * S1;
    link.move(dir, r / 2);
    for (int mu = 0; mu < 3; mu++) {
      if (dir != mu) {
        link.move(mu, -d_max);
        for (int d = -d_max; d < 0; d++) {
          S2 = schwinger_lines[abs(d) - 1][link.place / 4 * 3 + mu];
          A2 = A1 ^ S2;
          A2 = A2 * plaket[link.place / 4 * 3 + dir];
          A2 = A2 * S2;
          correlator[d + d_max] += A2.multiply_conj_tr(S1);
          link.move(mu, 1);
        }
        A2 = A1 * plaket[link.place / 4 * 3 + dir];
        correlator[d_max] += A2.multiply_conj_tr(S1);
        for (int d = 1; d <= d_max; d++) {
          S2 = schwinger_lines[d - 1][link.place / 4 * 3 + mu];
          link.move(mu, d);
          A2 = A1 * S2;
          A2 = A2 * plaket[link.place / 4 * 3 + dir];
          A2 = A2 ^ S2;
          correlator[d + d_max] += A2.multiply_conj_tr(S1);
          link.move(mu, -d);
        }
      }
    }
    link.go_update(x, y, z, t);
    link.move(dir, r);
    W = wilson_loops_opposite[link.place / 4 * 3 + dir];
    link.move(dir, 1 - r / 2);
    S1 = schwinger_lines[r / 2 - 2][link.place / 4 * 3 + dir].conj();
    A1 = W * S1;
    for (int mu = 0; mu < 3; mu++) {
      if (dir != mu) {
        link.move(mu, -d_max);
        for (int d = -d_max; d < 0; d++) {
          S2 = schwinger_lines[abs(d) - 1][link.place / 4 * 3 + mu];
          A2 = A1 ^ S2;
          A2 = A2 ^ plaket[link.place / 4 * 3 + dir];
          A2 = A2 * S2;
          correlator[d + d_max] += A2.multiply_conj_tr(S1);
          link.move(mu, 1);
        }
        A2 = A1 * plaket[link.place / 4 * 3 + dir];
        correlator[d_max] += A2.multiply_conj_tr(S1);
        for (int d = 1; d <= d_max; d++) {
          S2 = schwinger_lines[d - 1][link.place / 4 * 3 + mu];
          link.move(mu, d);
          A2 = A1 * S2;
          A2 = A2 ^ plaket[link.place / 4 * 3 + dir];
          A2 = A2 ^ S2;
          correlator[d + d_max] += A2.multiply_conj_tr(S1);
          link.move(mu, -d);
        }
      }
    }
  }
  SPACE_ITER_END

  std::map<int, double> result;
  for (int i = 0; i < correlator.size(); i++) {
    result[i - d_max] =
        correlator[i] / (x_size * y_size * z_size * t_size * 12);
  }
  return result;
}

template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_transversal_tr(
    std::vector<T> &array, std::vector<std::vector<T>> &schwinger_lines_short,
    int T_min, int T_max, int R_min, int R_max, int d_max) {

  std::vector<T> plaket_time = calculate_plaket_schwinger_time_tr(array);

  std::vector<T> wilson_loops;
  std::vector<T> wilson_loops_opposite;
  std::map<std::tuple<int, int, int>, double> result;
  std::map<int, double> schwinger_electric;

  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {

      wilson_loops = calculate_wilson_loops_schwinger(array, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(array, r, t);

      if (r % 2 == 0) {
        schwinger_electric =
            wilson_plaket_schwinger_electric_transversal_tr_even(
                wilson_loops, plaket_time, schwinger_lines_short, d_max, t, r);
      } else {
        schwinger_electric =
            wilson_plaket_schwinger_electric_transversal_tr_odd(
                wilson_loops, wilson_loops_opposite, plaket_time,
                schwinger_lines_short, d_max, t, r);
      }
      for (auto it = schwinger_electric.begin(); it != schwinger_electric.end();
           ++it) {
        result[std::tuple<int, int, double>(t, r, it->first)] = it->second;
      }
    }
  }

  return result;
}

template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_longitudinal_l(
    std::vector<T> &array, std::vector<std::vector<T>> &schwinger_lines_short,
    int T_min, int T_max, int R_min, int R_max, int d_ouside) {

  std::vector<T> plaket_space = calculate_plaket_schwinger_space_l(array);

  std::vector<T> wilson_loops;
  std::vector<T> wilson_loops_opposite;
  std::map<std::tuple<int, int, int>, double> result;
  std::map<int, double> schwinger_electric;

  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {

      wilson_loops = calculate_wilson_loops_schwinger(array, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(array, r, t);

      schwinger_electric = wilson_plaket_schwinger_longitudinal_tr(
          wilson_loops, wilson_loops_opposite, plaket_space,
          schwinger_lines_short, d_ouside, t, r);
      for (auto it = schwinger_electric.begin(); it != schwinger_electric.end();
           ++it) {
        result[std::tuple<int, int, double>(t, r, it->first)] = it->second;
      }
    }
  }

  return result;
}

template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_longitudinal_tr(
    std::vector<T> &array, std::vector<std::vector<T>> &schwinger_lines_short,
    int T_min, int T_max, int R_min, int R_max, int d_ouside) {

  std::vector<T> plaket_space = calculate_plaket_schwinger_space_tr(array);

  std::vector<T> wilson_loops;
  std::vector<T> wilson_loops_opposite;
  std::map<std::tuple<int, int, int>, double> result;
  std::map<int, double> schwinger_electric;

  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {

      wilson_loops = calculate_wilson_loops_schwinger(array, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(array, r, t);

      schwinger_electric = wilson_plaket_schwinger_longitudinal_tr(
          wilson_loops, wilson_loops_opposite, plaket_space,
          schwinger_lines_short, d_ouside, t, r);
      for (auto it = schwinger_electric.begin(); it != schwinger_electric.end();
           ++it) {
        result[std::tuple<int, int, double>(t, r, it->first)] = it->second;
      }
    }
  }

  return result;
}

template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_transversal_l(
    std::vector<T> &array, std::vector<std::vector<T>> &schwinger_lines_short,
    int T_min, int T_max, int R_min, int R_max, int d_max) {

  std::vector<T> plaket_space = calculate_plaket_schwinger_space_l(array);

  std::vector<T> wilson_loops;
  std::vector<T> wilson_loops_opposite;
  std::map<std::tuple<int, int, int>, double> result;
  std::map<int, double> schwinger_electric;

  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {

      wilson_loops = calculate_wilson_loops_schwinger(array, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(array, r, t);

      if (r % 2 == 0) {
        schwinger_electric =
            wilson_plaket_schwinger_electric_transversal_tr_even(
                wilson_loops, plaket_space, schwinger_lines_short, d_max, t, r);
      } else {
        schwinger_electric =
            wilson_plaket_schwinger_electric_transversal_tr_odd(
                wilson_loops, wilson_loops_opposite, plaket_space,
                schwinger_lines_short, d_max, t, r);
      }
      for (auto it = schwinger_electric.begin(); it != schwinger_electric.end();
           ++it) {
        result[std::tuple<int, int, double>(t, r, it->first)] = it->second;
      }
    }
  }

  return result;
}

template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_transversal_tr(
    std::vector<T> &array, std::vector<std::vector<T>> &schwinger_lines_short,
    int T_min, int T_max, int R_min, int R_max, int d_max) {

  std::vector<T> plaket_space = calculate_plaket_schwinger_space_tr(array);

  std::vector<T> wilson_loops;
  std::vector<T> wilson_loops_opposite;
  std::map<std::tuple<int, int, int>, double> result;
  std::map<int, double> schwinger_electric;

  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {

      wilson_loops = calculate_wilson_loops_schwinger(array, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(array, r, t);

      if (r % 2 == 0) {
        schwinger_electric =
            wilson_plaket_schwinger_electric_transversal_tr_even(
                wilson_loops, plaket_space, schwinger_lines_short, d_max, t, r);
      } else {
        schwinger_electric =
            wilson_plaket_schwinger_electric_transversal_tr_odd(
                wilson_loops, wilson_loops_opposite, plaket_space,
                schwinger_lines_short, d_max, t, r);
      }
      for (auto it = schwinger_electric.begin(); it != schwinger_electric.end();
           ++it) {
        result[std::tuple<int, int, double>(t, r, it->first)] = it->second;
      }
    }
  }

  return result;
}

template <class T>
void flux_schwinger_all(std::vector<T> &array,
                        std::vector<std::vector<T>> &schwinger_lines_short,
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
                        int T_min, int T_max, int R_min, int R_max,
                        int d_ouside, int d_max) {

  std::vector<T> plaket_time_left = calculate_plaket_schwinger_time_left(array);
  std::vector<T> plaket_time_right =
      calculate_plaket_schwinger_time_right(array);
  std::vector<T> plaket_time_tr = calculate_plaket_schwinger_time_tr(array);
  std::vector<T> plaket_space_l = calculate_plaket_schwinger_space_l(array);
  std::vector<T> plaket_space_tr = calculate_plaket_schwinger_space_tr(array);

  std::vector<T> wilson_loops;
  std::vector<T> wilson_loops_opposite;
  std::map<int, double> schwinger_correlator;

  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 1) {

      wilson_loops = calculate_wilson_loops_schwinger(array, r, t);
      wilson_loops_opposite =
          calculate_wilson_loops_schwinger_opposite(array, r, t);

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

template <class T>
double wilson_plaket_electric_longitudinal(const std::vector<T> &array,
                                           const std::vector<T> &plaket,
                                           int time, int r) {
  link1 link(x_size, y_size, z_size, t_size);

  std::vector<T> wilson_loops =
      calculate_wilson_loops_schwinger(array, r, time);
  int place = 0;
  double correlator = 0;

  SPACE_ITER_START

  place = link.place / 4 * 3;

  for (int dir = 0; dir < 3; dir++) {
    correlator += wilson_loops[place + dir].tr() * plaket[place + dir].tr();
  }

  SPACE_ITER_END

  return correlator;
}

template <class T>
std::vector<double> calculate_plaket_time_trace_l(const std::vector<T> &array) {
  std::vector<double> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
#pragma omp parallel for collapse(4) private(link)
  SPACE_ITER_START;
  for (int dir = 0; dir < 3; dir++) {
    link.move_dir(dir);
    vec[link.place / 4 * 3 + dir] = (link.plaket_left_down(array, 3).tr() +
                                     link.plaket_left_up(array, 3).tr()) /
                                    2;
  }
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<double>
calculate_plaket_time_trace_tr(const std::vector<T> &array) {
  std::vector<double> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
#pragma omp parallel for collapse(4) private(link)
  SPACE_ITER_START;
  for (int dir = 0; dir < 3; dir++) {
    for (int mu = 0; mu < 3; mu++) {
      if (dir != mu) {
        link.move_dir(mu);
        vec[link.place / 4 * 3 + dir] +=
            (link.plaket_left_down(array, 3).tr() +
             link.plaket_left_up(array, 3).tr() +
             link.plaket_right_down(array, 3).tr() +
             link.plaket_right_up(array, 3).tr()) /
            8;
      }
    }
  }
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<double> calculate_plaket_space_trace_l(const std::vector<T> array) {
  std::vector<double> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
#pragma omp parallel for collapse(4) private(link)
  SPACE_ITER_START;
  for (int dir = 0; dir < 3; dir++) {
    for (int mu = 0; mu < 3; mu++) {
      if (dir != mu) {
        link.move_dir(mu);
        vec[link.place / 4 * 3 + dir] = vec[link.place / 4 * 3 + dir] +
                                        (link.plaket_right_down(array, 3).tr() +
                                         link.plaket_right_up(array, 3).tr() +
                                         link.plaket_left_down(array, 3).tr() +
                                         link.plaket_left_up(array, 3).tr()) *
                                            0.125;
      }
    }
  }
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<double>
calculate_plaket_space_trace_tr(const std::vector<T> &array) {
  std::vector<double> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
#pragma omp parallel for collapse(4) private(link)
  SPACE_ITER_START;
  for (int i = 0; i < 3; i++) {
    for (int j = i + 1; j < 3; j++) {
      if (i != j) {
        link.move_dir(j);
        for (int k = 0; k < 3; k++) {
          if (k != i && k != j) {
            vec[link.place / 4 * 3 + i] =
                vec[link.place / 4 * 3 + i] +
                (link.plaket_right_down(array, k).tr() +
                 link.plaket_right_up(array, k).tr() +
                 link.plaket_left_down(array, k).tr() +
                 link.plaket_left_up(array, k).tr()) *
                    0.25;
          }
        }
      }
    }
  }
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<double> calculate_plaket_space_tr(const std::vector<T> &array) {
  std::vector<double> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
#pragma omp parallel for collapse(4) private(link)
  SPACE_ITER_START;
  for (int dir = 0; dir < 3; dir++) {
    for (int j = dir + 1; j < 3; j++) {
      link.move_dir(dir);
      vec[link.place / 4 * 3 + dir + j] = link.plaket_left_down(array, j).tr();
    }
  }
  SPACE_ITER_END;
  return vec;
}

double plaket4_time(const std::vector<double> &plaket_tr, link1 &link) {
  double a = plaket_tr[link.place / 4 * 3 + link.direction];
  link.move(link.direction, -1);
  a += plaket_tr[link.place / 4 * 3 + link.direction];
  link.move(3, -1);
  a += plaket_tr[link.place / 4 * 3 + link.direction];
  link.move(link.direction, 1);
  a += plaket_tr[link.place / 4 * 3 + link.direction];
  link.move(3, 1);
  return a / 4;
}

double plaket4_space(const std::vector<double> &plaket_tr, link1 &link,
                     int nu) {
  double a = plaket_tr[link.place / 4 * 3 + link.direction + nu];
  link.move(link.direction, -1);
  a += plaket_tr[link.place / 4 * 3 + link.direction + nu];
  link.move(nu, -1);
  a += plaket_tr[link.place / 4 * 3 + link.direction + nu];
  link.move(link.direction, 1);
  a += plaket_tr[link.place / 4 * 3 + link.direction + nu];
  link.move(nu, 1);
  return a / 4;
}

template <class T>
std::vector<double> calculate_wilson_loop_tr(const std::vector<T> &array, int r,
                                             int time) {
  std::vector<double> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
#pragma omp parallel for collapse(4) private(link)
  SPACE_ITER_START;
  for (int dir = 0; dir < 3; dir++) {
    link.move_dir(dir);
    vec[link.place / 4 * 3 + dir] = link.wilson_loop(array, r, time).tr();
  }
  SPACE_ITER_END;
  return vec;
}

std::map<int, double> wilson_plaket_correlator_electric_longitudinal(
    const std::vector<double> &wilson_loop_tr,
    const std::vector<double> &plaket_tr, int r, int time, int d_ouside) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<double> correlator(2 * d_ouside + r + 1);
  double a;
  SPACE_ITER_START
  for (int dir = 0; dir < 3; dir++) {
    link.go_update(x, y, z, t);
    a = wilson_loop_tr[link.place / 4 * 3 + dir];
    link.move(3, time / 2);
    link.move(dir, -d_ouside);
    for (int d = -d_ouside; d <= d_ouside + r; d++) {
      correlator[d + d_ouside] += a * plaket_tr[link.place / 4 * 3 + dir];
      link.move(dir, 1);
    }
  }
  SPACE_ITER_END
  int count;
  count = data_size / 4 * 3;
  std::map<int, double> result;
  for (int i = 0; i < correlator.size(); i++) {
    result[i - d_ouside - r / 2] = correlator[i] / count;
  }
  return result;
}

std::map<int, double> wilson_plaket_correlator_electric_transversal(
    const std::vector<double> &wilson_loop_tr,
    const std::vector<double> &plaket_tr, int r, int time, int d_max) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<double> correlator(2 * d_max + 1);
  double a;
  SPACE_ITER_START
  for (int dir = 0; dir < 3; dir++) {
    link.go_update(x, y, z, t);
    a = wilson_loop_tr[link.place / 4 * 3 + dir];
    link.move(3, time / 2);
    link.move(dir, r / 2 - 1);
    for (int mu = 0; mu < 3; mu++) {
      if (dir != mu) {
        link.move(mu, -d_max);
        for (int d = -d_max; d <= d_max; d++) {
          correlator[d + d_max] += a * plaket_tr[link.place / 4 * 3 + dir];
          link.move(dir, 1);
        }
        link.move(mu, -(d_max + 1));
      }
    }
  }
  SPACE_ITER_END
  int count;
  count = data_size / 4 * 6;
  std::map<int, double> result;
  for (int i = 0; i < correlator.size(); i++) {
    result[i - d_max] = correlator[i] / count;
  }
  return result;
}

template <class T>
std::map<std::tuple<int, int, int>, double>
calculate_wilson_plaket_correlator_electric_longitudinal(
    const std::vector<double> &plaket_tr, const std::vector<T> &conf_wilson,
    int T_min, int T_max, int R_min, int R_max, int d_ouside) {

  std::map<std::tuple<int, int, int>, double> flux_tube;

  std::vector<double> wilson_loop_tr;
  std::map<int, double> correlator;

  for (int time = T_min; time <= T_max; time += 2) {
    for (int r = R_min; r <= R_max; r += 1) {

      wilson_loop_tr = calculate_wilson_loop_tr(conf_wilson, r, time);

      correlator = wilson_plaket_correlator_electric_longitudinal(
          wilson_loop_tr, plaket_tr, r, time, d_ouside);

      for (auto i = correlator.begin(); i != correlator.end(); i++) {
        flux_tube[std::tuple<int, int, int>(time, r, i->first)] = i->second;
      }
    }
  }

  return flux_tube;
}

template <class T>
std::map<std::tuple<int, int, int>, double>
calculate_wilson_plaket_correlator_electric_transversal(
    const std::vector<double> &plaket_tr, const std::vector<T> &conf_wilson,
    int T_min, int T_max, int R_min, int R_max, int d_ouside) {

  std::map<std::tuple<int, int, int>, double> flux_tube;

  std::vector<double> wilson_loop_tr;
  std::map<int, double> correlator;

  for (int time = T_min; time <= T_max; time += 2) {
    for (int r = R_min; r <= R_max; r += 1) {

      wilson_loop_tr = calculate_wilson_loop_tr(conf_wilson, r, time);

      correlator = wilson_plaket_correlator_electric_transversal(
          wilson_loop_tr, plaket_tr, r, time, d_ouside);

      for (auto i = correlator.begin(); i != correlator.end(); i++) {
        flux_tube[std::tuple<int, int, int>(time, r, i->first)] = i->second;
      }
    }
  }

  return flux_tube;
}

std::map<int, double>
wilson_plaket_correlator_electric(const std::vector<double> &wilson_loop_tr,
                                  const std::vector<double> &plaket_tr, int r,
                                  int time, int x_trans, int d_min, int d_max) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<double> correlator(d_max - d_min + 1, 0.0);
  double a;
  bool if_test = false;
  for (int dir = 0; dir < 3; dir++) {
    SPACE_ITER_START
    a = wilson_loop_tr[link.place / 4 * 3 + dir];
    link.move(3, time / 2);
    if_test = false;
    link.move(dir, d_min);
    for (int d = d_min; d <= d_max; d++) {
      if (x_trans == 0) {
        for (int mu = 0; mu < 3; mu++) {
          link.move_dir(mu);
          correlator[d - d_min] += a * plaket4_time(plaket_tr, link);
        }
      } else {
        for (int nu = 0; nu < 3; nu++) {
          if (nu != dir) {
            link.move(nu, x_trans);
            for (int mu = 0; mu < 3; mu++) {
              link.move_dir(mu);
              correlator[d - d_min] += a * plaket4_time(plaket_tr, link);
            }
            link.move(nu, -2 * x_trans);
            for (int mu = 0; mu < 3; mu++) {
              link.move_dir(mu);
              correlator[d - d_min] += a * plaket4_time(plaket_tr, link);
            }
            link.move(nu, x_trans);
          }
        }
      }
      link.move(dir, 1);
    }
    SPACE_ITER_END
  }
  int count;
  if (x_trans == 0)
    count = data_size / 4 * 9;
  else
    count = data_size / 4 * 36;
  std::map<int, double> result;
  for (int i = 0; i < correlator.size(); i++) {
    result[i + d_min] = correlator[i] / count;
  }
  return result;
}

std::map<int, double>
wilson_plaket_correlator_electric_x(const std::vector<double> &wilson_loop_tr,
                                    const std::vector<double> &plaket_tr, int r,
                                    int time, int x_trans_max, int d) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<double> correlator(x_trans_max + 1, 0.0);
  double a;
  for (int dir = 0; dir < 3; dir++) {
    SPACE_ITER_START
    a = wilson_loop_tr[link.place / 4 * 3 + dir];
    link.move(3, time / 2);
    link.move(dir, d);
    for (int nu = 0; nu < 3; nu++) {
      if (nu != dir) {
        link.move(nu, -x_trans_max);
        for (int x_trans = x_trans_max; x_trans > 0; --x_trans) {
          for (int mu = 0; mu < 3; mu++) {
            link.move_dir(mu);
            correlator[x_trans] += a * plaket4_time(plaket_tr, link);
          }
          link.move(nu, 1);
        }
        for (int mu = 0; mu < 3; mu++) {
          link.move_dir(mu);
          correlator[0] += a * plaket4_time(plaket_tr, link);
        }
        for (int x_trans = 1; x_trans <= x_trans_max; ++x_trans) {
          link.move(nu, 1);
          for (int mu = 0; mu < 3; mu++) {
            link.move_dir(mu);
            correlator[x_trans] += a * plaket4_time(plaket_tr, link);
          }
        }
        link.move(nu, -x_trans_max);
      }
    }
    SPACE_ITER_END
  }
  int count;
  count = x_size * y_size * z_size * t_size * 36;
  std::map<int, double> result;
  result[0] = correlator[0] / (x_size * y_size * z_size * t_size * 18);
  for (int i = 1; i < correlator.size(); i++) {
    result[i] = correlator[i] / count;
  }
  return result;
}

std::map<int, double>
wilson_plaket_correlator_magnetic(const std::vector<double> &wilson_loop_tr,
                                  const std::vector<double> &plaket_tr, int r,
                                  int time, int x_trans, int d_min, int d_max) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<double> correlator(d_max - d_min + 1, 0.0);
  double a;
  for (int dir = 0; dir < 3; dir++) {
    SPACE_ITER_START
    link.move_dir(dir);
    a = wilson_loop_tr[link.place / 4 * 3 + dir];
    link.move(3, time / 2);
    link.move(dir, d_min);
    for (int d = d_min; d <= d_max; d++) {
      if (x_trans == 0) {
        for (int mu = 0; mu < 3; mu++) {
          for (int j = mu + 1; j < 3; j++) {
            link.move_dir(mu);
            correlator[d - d_min] += a * plaket4_space(plaket_tr, link, j);
          }
        }
      } else {
        for (int nu = 0; nu < 3; nu++) {
          if (nu != dir) {
            link.move(nu, x_trans);
            for (int mu = 0; mu < 3; mu++) {
              for (int j = mu + 1; j < 3; j++) {
                link.move_dir(mu);
                correlator[d - d_min] += a * plaket4_space(plaket_tr, link, j);
              }
            }
            link.move(nu, -2 * x_trans);
            for (int mu = 0; mu < 3; mu++) {
              for (int j = mu + 1; j < 3; j++) {
                link.move_dir(mu);
                correlator[d - d_min] += a * plaket4_space(plaket_tr, link, j);
              }
            }
            link.move(nu, x_trans);
          }
        }
      }
      link.move(dir, 1);
    }
    SPACE_ITER_END
  }
  int count;
  if (x_trans == 0)
    count = data_size / 4 * 9;
  else
    count = data_size / 4 * 36;
  std::map<int, double> result;
  for (int i = 0; i < correlator.size(); i++) {
    result[i + d_min] = correlator[i] / count;
  }
  return result;
}

std::map<int, double>
wilson_plaket_correlator_magnetic_x(const std::vector<double> &wilson_loop_tr,
                                    const std::vector<double> &plaket_tr, int R,
                                    int T, int x_trans_max, int d) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<double> correlator(x_trans_max + 1, 0.0);
  double a;
  for (int dir = 0; dir < 3; dir++) {
    SPACE_ITER_START
    a = wilson_loop_tr[link.place / 4 * 3 + dir];
    link.move(3, T / 2);
    link.move(dir, d);
    for (int nu = 0; nu < 3; nu++) {
      if (nu != dir) {
        link.move(nu, -x_trans_max);
        for (int x_trans = x_trans_max; x_trans > 0; --x_trans) {
          for (int mu = 0; mu < 3; mu++) {
            for (int j = mu + 1; j < 3; j++) {
              link.move_dir(mu);
              correlator[x_trans] += a * plaket4_space(plaket_tr, link, j);
            }
          }
          link.move(nu, 1);
        }
        for (int mu = 0; mu < 3; mu++) {
          for (int j = mu + 1; j < 3; j++) {
            link.move_dir(mu);
            correlator[0] += a * plaket4_space(plaket_tr, link, j);
          }
        }
        for (int x_trans = 1; x_trans <= x_trans_max; ++x_trans) {
          link.move(nu, 1);
          for (int mu = 0; mu < 3; mu++) {
            for (int j = mu + 1; j < 3; j++) {
              link.move_dir(mu);
              correlator[x_trans] += a * plaket4_space(plaket_tr, link, j);
            }
          }
        }
        link.move(nu, -x_trans_max);
      }
    }
    SPACE_ITER_END
  }
  int count;
  count = x_size * y_size * z_size * t_size * 36;
  std::map<int, double> result;
  result[0] = correlator[0] / (x_size * y_size * z_size * t_size * 18);
  for (int i = 1; i < correlator.size(); i++) {
    result[i] = correlator[i] / count;
  }
  return result;
}

template <class T>
std::vector<double>
wilson_plane_tr(std::vector<T> &wilson_lines_mu,
                std::vector<T> &wilson_lines_nu, int size_mu1, int size_mu2,
                int size_nu1, int size_nu2, int length_mu, int length_nu) {
  int data_size1 = x_size * y_size * z_size * t_size;

  std::vector<double> wilson_loops_tr(x_size * y_size * z_size * t_size, 0.0);

  T loops;

#pragma omp parallel for collapse(3) private(loops)
  for (int k = 0; k < data_size1; k += size_nu2) {
    for (int i = 0; i < size_nu2; i += size_mu2) {
      for (int j = 0; j < size_mu2; j++) {
        if (j < size_mu2 - length_mu * size_mu1)
          loops = wilson_lines_mu[i + k + j] *
                  wilson_lines_nu[i + k + j + length_mu * size_mu1];
        else
          loops = wilson_lines_mu[i + k + j] *
                  wilson_lines_nu[i + k + j - size_mu2 + length_mu * size_mu1];
        if (i + j < size_nu2 - length_nu * size_nu1)
          loops = loops ^ wilson_lines_mu[i + k + j + length_nu * size_nu1];
        else
          loops = loops ^
                  wilson_lines_mu[i + k + j - size_nu2 + length_nu * size_nu1];

        if (k + i + length_nu / 2 * size_nu1 >= size_nu2)
          wilson_loops_tr[i + k + j + length_nu / 2 * size_nu1 - size_nu2] =
              loops.multiply_conj_tr(wilson_lines_nu[i + k + j]);
        else
          wilson_loops_tr[i + k + j + length_nu / 2 * size_nu1] =
              loops.multiply_conj_tr(wilson_lines_nu[i + k + j]);
      }
    }
  }

  return wilson_loops_tr;
}

template <class T>
std::vector<double> plaket_plane_tr(std::vector<T> &conf_mu,
                                    std::vector<T> &conf_nu, int size_mu1,
                                    int size_mu2, int size_nu1, int size_nu2) {
  int data_size1 = x_size * y_size * z_size * t_size;

  std::vector<double> wilson_loops_tr(x_size * y_size * z_size * t_size);

  T loops;

  // #pragma omp parallel for collapse(3) private(loops) num_threads(1)
  for (int k = 0; k < data_size1; k += size_nu2) {
    for (int i = 0; i < size_nu2; i += size_mu2) {
      for (int j = 0; j < size_mu2; j++) {
        if (j < size_mu2 - size_mu1)
          loops = conf_mu[i + k + j] * conf_nu[i + k + j + size_mu1];
        else
          loops = conf_mu[i + k + j] * conf_nu[i + k + j - size_mu2 + size_mu1];
        if (i + j < size_nu2 - size_nu1)
          loops = loops ^ conf_mu[i + k + j + size_nu1];
        else
          loops = loops ^ conf_mu[i + k + j - size_nu2 + size_nu1];

        wilson_loops_tr[i + k + j] = loops.multiply_conj_tr(conf_nu[i + k + j]);
      }
    }
  }

  return wilson_loops_tr;
}

void plaket_plane_aver(std::vector<double> &plaket_aver_tr,
                       std::vector<double> &plaket_tr, int size_mu1,
                       int size_mu2, int size_nu1, int size_nu2) {

  int data_size1 = x_size * y_size * z_size * t_size;

  // #pragma omp parallel for collapse(3) num_threads(1)
  for (int k = 0; k < data_size1; k += size_nu2) {
    for (int i = 0; i < size_nu2; i += size_mu2) {
      for (int j = 0; j < size_mu2; j++) {

        if (i >= size_nu1) {
          plaket_aver_tr[i + k + j] +=
              plaket_tr[i + k + j] + plaket_tr[i + k + j - size_nu1];

          if (j >= size_mu1) {
            plaket_aver_tr[i + k + j] +=
                plaket_tr[i + k + j - size_mu1] +
                plaket_tr[i + k + j - size_mu1 - size_nu1];
          }

          else {
            plaket_aver_tr[i + k + j] +=
                plaket_tr[i + k + j + size_mu2 - size_mu1] +
                plaket_tr[i + k + j + size_mu2 - size_mu1 - size_nu1];
          }
        } else {
          plaket_aver_tr[i + k + j] +=
              plaket_tr[i + k + j] + plaket_tr[i + k + j + size_nu2 - size_nu1];

          if (j >= size_mu1) {
            plaket_aver_tr[i + k + j] +=
                plaket_tr[i + k + j - size_mu1] +
                plaket_tr[i + k + j - size_mu1 + size_nu2 - size_nu1];
          } else {
            plaket_aver_tr[i + k + j] +=
                plaket_tr[i + k + j + size_mu2 - size_mu1] +
                plaket_tr[i + k + j + size_mu2 - size_mu1 + size_nu2 -
                          size_nu1];
          }
        }
      }
    }
  }
}

template <class T>
std::vector<double> plaket_aver_tr_time(std::vector<std::vector<T>> conf) {

  int data_size1 = x_size * y_size * z_size * t_size;

  std::vector<int> steps = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  std::vector<double> plaket_result(data_size1, 0.0);
  std::vector<double> plaket_plane;

  for (int mu = 0; mu < 3; mu++) {

    plaket_plane = plaket_plane_tr(conf[mu], conf[3], steps[mu], steps[mu + 1],
                                   steps[3], steps[4]);

    plaket_plane_aver(plaket_result, plaket_plane, steps[mu], steps[mu + 1],
                      steps[3], steps[4]);
  }

  for (auto it = plaket_result.begin(); it != plaket_result.end(); it++) {
    *it = *it / 12;
  }

  return plaket_result;
}

template <class T>
std::vector<double> plaket_aver_tr_space(std::vector<std::vector<T>> conf) {

  int data_size1 = x_size * y_size * z_size * t_size;

  std::vector<int> steps = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  std::vector<double> plaket_result(data_size1, 0.0);
  std::vector<double> plaket_plane;

  for (int mu = 0; mu < 3; mu++) {
    for (int nu = mu + 1; nu < 3; nu++) {

      plaket_plane = plaket_plane_tr(conf[mu], conf[nu], steps[mu],
                                     steps[mu + 1], steps[nu], steps[nu + 1]);

      plaket_plane_aver(plaket_result, plaket_plane, steps[mu], steps[mu + 1],
                        steps[nu], steps[nu + 1]);
    }
  }

  for (auto it = plaket_result.begin(); it != plaket_result.end(); it++) {
    *it = *it / 12;
  }

  return plaket_result;
}

void wilson_plaket_correlator_plane_longitudinal(
    std::vector<double> &correlator, const std::vector<double> &wilson_loop_tr,
    const std::vector<double> &plaket_tr, int size_mu1, int size_mu2,
    int x_trans, int d_min, int d_max) {
  int data_size1 = x_size * y_size * z_size * t_size;

  for (int i = 0; i < data_size1; i += size_mu2) {
    for (int j = 0; j < size_mu2; j++) {
      for (int d = d_min; d <= d_max; d++) {
        if (j + d * size_mu1 < 0) {
          correlator[d - d_min] += wilson_loop_tr[i + j] *
                                   plaket_tr[i + j + d * size_mu1 + size_mu2];
        } else if (j + d * size_mu1 < size_mu2) {
          correlator[d - d_min] +=
              wilson_loop_tr[i + j] * plaket_tr[i + j + d * size_mu1];
        } else {
          correlator[d - d_min] += wilson_loop_tr[i + j] *
                                   plaket_tr[i + j + d * size_mu1 - size_mu2];
        }
      }
    }
  }
}

void wilson_plaket_correlator_plane_transversal(
    std::vector<double> &correlator, const std::vector<double> &wilson_loop_tr,
    const std::vector<double> &plaket_tr, int size_mu1, int size_mu2,
    int size_nu1, int size_nu2, int d, int x_trans_min, int x_trans_max) {
  int data_size1 = x_size * y_size * z_size * t_size;

  int mu_shift;
  int nu_shift;

  if (size_mu2 < size_nu2) {

    // #pragma omp parallel for collapse(3) private(mu_shift, nu_shift)
    // num_threads(4)
    for (int i = 0; i < data_size1; i += size_nu2) {
      for (int j = 0; j < size_nu2; j += size_mu2) {
        for (int k = 0; k < size_mu2; k++) {

          if (k + d * size_mu1 < 0)
            mu_shift = d * size_mu1 + size_mu2;
          else if (k + d * size_mu1 < size_mu2)
            mu_shift = d * size_mu1;
          else
            mu_shift = d * size_mu1 - size_mu2;

          for (int x = x_trans_min; x <= x_trans_max; x++) {

            if (j + x * size_nu1 < 0)
              nu_shift = x * size_nu1 + size_nu2;
            else if (j + x * size_nu1 < size_nu2)
              nu_shift = x * size_nu1;
            else
              nu_shift = x * size_nu1 - size_nu2;

            correlator[x - x_trans_min] +=
                wilson_loop_tr[i + j + k] *
                plaket_tr[i + j + k + nu_shift + mu_shift];
          }
        }
      }
    }

  } else {

    // #pragma omp parallel for collapse(3) private(mu_shift, nu_shift)
    // num_threads(4)
    for (int i = 0; i < data_size1; i += size_mu2) {
      for (int j = 0; j < size_mu2; j += size_nu2) {
        for (int k = 0; k < size_nu2; k++) {

          if (j + d * size_mu1 < 0)
            mu_shift = d * size_mu1 + size_mu2;
          else if (j + d * size_mu1 < size_mu2)
            mu_shift = d * size_mu1;
          else
            mu_shift = d * size_mu1 - size_mu2;

          for (int x = x_trans_min; x <= x_trans_max; x++) {

            if (k + x * size_nu1 < 0)
              nu_shift = x * size_nu1 + size_nu2;
            else if (k + x * size_nu1 < size_nu2)
              nu_shift = x * size_nu1;
            else
              nu_shift = x * size_nu1 - size_nu2;

            correlator[x - x_trans_min] +=
                wilson_loop_tr[i + j + k] *
                plaket_tr[i + j + k + nu_shift + mu_shift];
          }
        }
      }
    }
  }
}

template <class T>
std::map<std::tuple<int, int, int>, double>
wilson_plaket_correlator(std::vector<double> &plaket_tr,
                         std::vector<std::vector<T>> &conf_wilson, int T_min,
                         int T_max, int R_min, int R_max, int main_coordinate,
                         int transverse_coordinate, std::string direction) {
  int data_size1 = x_size * y_size * z_size * t_size;

  std::vector<int> steps = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  std::map<std::tuple<int, int, int>, double> flux_tube;

  std::map<int, std::vector<T>> time_lines;
  std::vector<T> space_lines;

  for (int t = T_min; t <= T_max; t += 2) {
    time_lines[t] = wilson_lines(conf_wilson[3], t, steps[3], steps[4]);
  }

  int transverse_shift = transverse_coordinate;

  std::vector<double> wilson_tr;

  for (int t = T_min; t <= T_max; t += 2) {
    for (int r = R_min; r <= R_max; r += 2) {

      int main_coordinate_min, main_coordinate_max;

      if (direction.compare("longitudinal") == 0) {
        main_coordinate_min = -main_coordinate;
        main_coordinate_max = r + main_coordinate;
      } else if (direction.compare("transversal") == 0) {
        main_coordinate_min = -main_coordinate;
        main_coordinate_max = main_coordinate;
        transverse_coordinate = r / 2 + transverse_shift;
      } else {
        std::cout << "wrong direction" << std::endl;
        exit(1);
      }

      std::vector<double> correlator(
          main_coordinate_max - main_coordinate_min + 1, 0.0);

      for (int mu = 0; mu < 3; mu++) {

        space_lines =
            wilson_lines(conf_wilson[mu], r, steps[mu], steps[mu + 1]);

        wilson_tr = wilson_plane_tr(space_lines, time_lines[t], steps[mu],
                                    steps[mu + 1], steps[3], steps[4], r, t);

        if (direction.compare("longitudinal") == 0) {
          wilson_plaket_correlator_plane_longitudinal(
              correlator, wilson_tr, plaket_tr, steps[mu], steps[mu + 1],
              transverse_coordinate, main_coordinate_min, main_coordinate_max);
        } else if (direction.compare("transversal") == 0) {
          for (int nu = 0; nu < 3; nu++) {
            if (mu != nu) {
              wilson_plaket_correlator_plane_transversal(
                  correlator, wilson_tr, plaket_tr, steps[mu], steps[mu + 1],
                  steps[nu], steps[nu + 1], transverse_coordinate,
                  main_coordinate_min, main_coordinate_max);
            }
          }
        }
      }

      int norm;
      if (direction.compare("longitudinal") == 0)
        norm = 3;

      else if (direction.compare("transversal") == 0)
        norm = 6;

      if (direction.compare("longitudinal") == 0) {
        for (int D = 0; D <= r / 2 + main_coordinate; D++) {
          flux_tube[std::tuple<int, int, int>(t, r, D)] =
              (correlator[D + main_coordinate + r / 2] +
               correlator[-D + main_coordinate + r / 2]) /
              data_size1 / norm / 2;
        }
      }

      else if (direction.compare("transversal") == 0) {
        for (int D = 0; D <= main_coordinate_max; D++) {
          flux_tube[std::tuple<int, int, int>(t, r, D)] =
              (correlator[D - main_coordinate_min] +
               correlator[-D - main_coordinate_min]) /
              data_size1 / norm / 2;
        }
      }
    }
  }

  return flux_tube;
}

// su2
template std::vector<su2>
calculate_plaket_time_left_down(const std::vector<su2> array);
template std::vector<su2>
calculate_plaket_time_left_up(const std::vector<su2> array);
template std::vector<su2>
calculate_plaket_time_right_down(const std::vector<su2> array);
template std::vector<su2>
calculate_plaket_time_right_up(const std::vector<su2> array);
template std::vector<su2>
calculate_wilson_loops_schwinger(const std::vector<su2> &array, int r,
                                 int time);
template std::vector<su2>
calculate_wilson_loops_schwinger_opposite(const std::vector<su2> &array, int r,
                                          int time);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_longitudinal_l(
    const std::vector<su2> &array,
    std::vector<std::vector<su2>> &schwinger_lines_short, int T_min, int T_max,
    int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_longitudinal_tr(
    std::vector<su2> &array,
    std::vector<std::vector<su2>> &schwinger_lines_short, int T_min, int T_max,
    int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_transversal_l(
    std::vector<su2> &array,
    std::vector<std::vector<su2>> &schwinger_lines_short, int T_min, int T_max,
    int R_min, int R_max, int d_max);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_transversal_tr(
    std::vector<su2> &array,
    std::vector<std::vector<su2>> &schwinger_lines_short, int T_min, int T_max,
    int R_min, int R_max, int d_max);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_longitudinal_l(
    std::vector<su2> &array,
    std::vector<std::vector<su2>> &schwinger_lines_short, int T_min, int T_max,
    int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_longitudinal_tr(
    std::vector<su2> &array,
    std::vector<std::vector<su2>> &schwinger_lines_short, int T_min, int T_max,
    int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_transversal_l(
    std::vector<su2> &array,
    std::vector<std::vector<su2>> &schwinger_lines_short, int T_min, int T_max,
    int R_min, int R_max, int d_max);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_transversal_tr(
    std::vector<su2> &array,
    std::vector<std::vector<su2>> &schwinger_lines_short, int T_min, int T_max,
    int R_min, int R_max, int d_max);
template void
flux_schwinger_all(std::vector<su2> &array,
                   std::vector<std::vector<su2>> &schwinger_lines_short,
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
                   int T_min, int T_max, int R_min, int R_max, int d_ouside,
                   int d_max);
template std::vector<su2>
calculate_schwinger_lines_short(const std::vector<su2> &array, int d);
template std::vector<double>
calculate_plaket_time_trace_l(const std::vector<su2> &array);
template std::vector<double>
calculate_plaket_time_trace_tr(const std::vector<su2> &array);
template std::vector<double>
calculate_plaket_space_trace_l(const std::vector<su2> array);
template std::vector<double>
calculate_plaket_space_trace_tr(const std::vector<su2> &array);
template std::vector<double>
calculate_plaket_space_tr(const std::vector<su2> &array);
template std::vector<su2>
calculate_polyakov_loop(const std::vector<su2> &array);
template std::vector<su2> calculate_wilson_loop(const std::vector<su2> &array,
                                                int r, int time);
template std::vector<double>
calculate_wilson_loop_tr(const std::vector<su2> &array, int r, int time);
template std::map<std::tuple<int, int, int>, double>
calculate_wilson_plaket_correlator_electric_longitudinal(
    const std::vector<double> &plaket_tr, const std::vector<su2> &conf_wilson,
    int T_min, int T_max, int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
calculate_wilson_plaket_correlator_electric_transversal(
    const std::vector<double> &plaket_tr, const std::vector<su2> &conf_wilson,
    int T_min, int T_max, int R_min, int R_max, int d_ouside);

template std::vector<double> wilson_plane_tr(std::vector<su2> &wilson_lines_mu,
                                             std::vector<su2> &wilson_lines_nu,
                                             int size_mu1, int size_mu2,
                                             int size_nu1, int size_nu2,
                                             int length_mu, int length_nu);

template std::vector<double> plaket_plane_tr(std::vector<su2> &conf_mu,
                                             std::vector<su2> &conf_nu,
                                             int size_mu1, int size_mu2,
                                             int size_nu1, int size_nu2);

template std::vector<double>
plaket_aver_tr_time(std::vector<std::vector<su2>> conf);
template std::vector<double>
plaket_aver_tr_space(std::vector<std::vector<su2>> conf);

template std::map<std::tuple<int, int, int>, double>
wilson_plaket_correlator(std::vector<double> &plaket_tr,
                         std::vector<std::vector<su2>> &conf_wilson, int T_min,
                         int T_max, int R_min, int R_max, int main_coordinate,
                         int transverse_coordinate, std::string direction);

// abelian

template std::vector<abelian>
calculate_plaket_time_left_down(const std::vector<abelian> array);
template std::vector<abelian>
calculate_plaket_time_left_up(const std::vector<abelian> array);
template std::vector<abelian>
calculate_plaket_time_right_down(const std::vector<abelian> array);
template std::vector<abelian>
calculate_plaket_time_right_up(const std::vector<abelian> array);
template std::vector<abelian>
calculate_wilson_loops_schwinger(const std::vector<abelian> &array, int r,
                                 int time);
template std::vector<abelian>
calculate_wilson_loops_schwinger_opposite(const std::vector<abelian> &array,
                                          int r, int time);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_longitudinal_l(
    const std::vector<abelian> &array,
    std::vector<std::vector<abelian>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_longitudinal_tr(
    std::vector<abelian> &array,
    std::vector<std::vector<abelian>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_transversal_l(
    std::vector<abelian> &array,
    std::vector<std::vector<abelian>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_max);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_transversal_tr(
    std::vector<abelian> &array,
    std::vector<std::vector<abelian>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_max);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_longitudinal_l(
    std::vector<abelian> &array,
    std::vector<std::vector<abelian>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_longitudinal_tr(
    std::vector<abelian> &array,
    std::vector<std::vector<abelian>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_transversal_l(
    std::vector<abelian> &array,
    std::vector<std::vector<abelian>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_max);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_transversal_tr(
    std::vector<abelian> &array,
    std::vector<std::vector<abelian>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_max);
template void
flux_schwinger_all(std::vector<abelian> &array,
                   std::vector<std::vector<abelian>> &schwinger_lines_short,
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
                   int T_min, int T_max, int R_min, int R_max, int d_ouside,
                   int d_max);
template std::vector<abelian>
calculate_schwinger_lines_short(const std::vector<abelian> &array, int d);
template std::vector<double>
calculate_plaket_time_trace_l(const std::vector<abelian> &array);
template std::vector<double>
calculate_plaket_time_trace_tr(const std::vector<abelian> &array);
template std::vector<double>
calculate_plaket_space_trace_l(const std::vector<abelian> array);
template std::vector<double>
calculate_plaket_space_trace_tr(const std::vector<abelian> &array);
template std::vector<double>
calculate_plaket_space_tr(const std::vector<abelian> &array);
template std::vector<abelian>
calculate_polyakov_loop(const std::vector<abelian> &array);
template std::vector<abelian>
calculate_wilson_loop(const std::vector<abelian> &array, int r, int time);
template std::vector<double>
calculate_wilson_loop_tr(const std::vector<abelian> &array, int r, int time);
template std::map<std::tuple<int, int, int>, double>
calculate_wilson_plaket_correlator_electric_longitudinal(
    const std::vector<double> &plaket_tr,
    const std::vector<abelian> &conf_wilson, int T_min, int T_max, int R_min,
    int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
calculate_wilson_plaket_correlator_electric_transversal(
    const std::vector<double> &plaket_tr,
    const std::vector<abelian> &conf_wilson, int T_min, int T_max, int R_min,
    int R_max, int d_ouside);

template std::vector<double>
wilson_plane_tr(std::vector<abelian> &wilson_lines_mu,
                std::vector<abelian> &wilson_lines_nu, int size_mu1,
                int size_mu2, int size_nu1, int size_nu2, int length_mu,
                int length_nu);

template std::vector<double> plaket_plane_tr(std::vector<abelian> &conf_mu,
                                             std::vector<abelian> &conf_nu,
                                             int size_mu1, int size_mu2,
                                             int size_nu1, int size_nu2);

template std::vector<double>
plaket_aver_tr_time(std::vector<std::vector<abelian>> conf);
template std::vector<double>
plaket_aver_tr_space(std::vector<std::vector<abelian>> conf);

template std::map<std::tuple<int, int, int>, double>
wilson_plaket_correlator(std::vector<double> &plaket_tr,
                         std::vector<std::vector<abelian>> &conf_wilson,
                         int T_min, int T_max, int R_min, int R_max,
                         int main_coordinate, int transverse_coordinate,
                         std::string direction);

// su3

template std::vector<su3>
calculate_plaket_time_left_down(const std::vector<su3> array);
template std::vector<su3>
calculate_plaket_time_left_up(const std::vector<su3> array);
template std::vector<su3>
calculate_plaket_time_right_down(const std::vector<su3> array);
template std::vector<su3>
calculate_plaket_time_right_up(const std::vector<su3> array);
template std::vector<su3>
calculate_wilson_loops_schwinger(const std::vector<su3> &array, int r,
                                 int time);
template std::vector<su3>
calculate_wilson_loops_schwinger_opposite(const std::vector<su3> &array, int r,
                                          int time);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_longitudinal_l(
    const std::vector<su3> &array,
    std::vector<std::vector<su3>> &schwinger_lines_short, int T_min, int T_max,
    int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_longitudinal_tr(
    std::vector<su3> &array,
    std::vector<std::vector<su3>> &schwinger_lines_short, int T_min, int T_max,
    int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_transversal_l(
    std::vector<su3> &array,
    std::vector<std::vector<su3>> &schwinger_lines_short, int T_min, int T_max,
    int R_min, int R_max, int d_max);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_transversal_tr(
    std::vector<su3> &array,
    std::vector<std::vector<su3>> &schwinger_lines_short, int T_min, int T_max,
    int R_min, int R_max, int d_max);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_longitudinal_l(
    std::vector<su3> &array,
    std::vector<std::vector<su3>> &schwinger_lines_short, int T_min, int T_max,
    int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_longitudinal_tr(
    std::vector<su3> &array,
    std::vector<std::vector<su3>> &schwinger_lines_short, int T_min, int T_max,
    int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_transversal_l(
    std::vector<su3> &array,
    std::vector<std::vector<su3>> &schwinger_lines_short, int T_min, int T_max,
    int R_min, int R_max, int d_max);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_transversal_tr(
    std::vector<su3> &array,
    std::vector<std::vector<su3>> &schwinger_lines_short, int T_min, int T_max,
    int R_min, int R_max, int d_max);
template void
flux_schwinger_all(std::vector<su3> &array,
                   std::vector<std::vector<su3>> &schwinger_lines_short,
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
                   int T_min, int T_max, int R_min, int R_max, int d_ouside,
                   int d_max);
template std::vector<su3>
calculate_schwinger_lines_short(const std::vector<su3> &array, int d);
template std::vector<double>
calculate_plaket_time_trace_l(const std::vector<su3> &array);
template std::vector<double>
calculate_plaket_time_trace_tr(const std::vector<su3> &array);
template std::vector<double>
calculate_plaket_space_trace_l(const std::vector<su3> array);
template std::vector<double>
calculate_plaket_space_trace_tr(const std::vector<su3> &array);
template std::vector<double>
calculate_plaket_space_tr(const std::vector<su3> &array);
template std::vector<su3>
calculate_polyakov_loop(const std::vector<su3> &array);
template std::vector<su3> calculate_wilson_loop(const std::vector<su3> &array,
                                                int r, int time);
template std::vector<double>
calculate_wilson_loop_tr(const std::vector<su3> &array, int r, int time);
template std::map<std::tuple<int, int, int>, double>
calculate_wilson_plaket_correlator_electric_longitudinal(
    const std::vector<double> &plaket_tr, const std::vector<su3> &conf_wilson,
    int T_min, int T_max, int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
calculate_wilson_plaket_correlator_electric_transversal(
    const std::vector<double> &plaket_tr, const std::vector<su3> &conf_wilson,
    int T_min, int T_max, int R_min, int R_max, int d_ouside);

template std::vector<double> wilson_plane_tr(std::vector<su3> &wilson_lines_mu,
                                             std::vector<su3> &wilson_lines_nu,
                                             int size_mu1, int size_mu2,
                                             int size_nu1, int size_nu2,
                                             int length_mu, int length_nu);

template std::vector<double> plaket_plane_tr(std::vector<su3> &conf_mu,
                                             std::vector<su3> &conf_nu,
                                             int size_mu1, int size_mu2,
                                             int size_nu1, int size_nu2);

template std::vector<double>
plaket_aver_tr_time(std::vector<std::vector<su3>> conf);
template std::vector<double>
plaket_aver_tr_space(std::vector<std::vector<su3>> conf);

template std::map<std::tuple<int, int, int>, double>
wilson_plaket_correlator(std::vector<double> &plaket_tr,
                         std::vector<std::vector<su3>> &conf_wilson, int T_min,
                         int T_max, int R_min, int R_max, int main_coordinate,
                         int transverse_coordinate, std::string direction);

// su3_abelian

template std::vector<su3_abelian>
calculate_plaket_time_left_down(const std::vector<su3_abelian> array);
template std::vector<su3_abelian>
calculate_plaket_time_left_up(const std::vector<su3_abelian> array);
template std::vector<su3_abelian>
calculate_plaket_time_right_down(const std::vector<su3_abelian> array);
template std::vector<su3_abelian>
calculate_plaket_time_right_up(const std::vector<su3_abelian> array);
template std::vector<su3_abelian>
calculate_wilson_loops_schwinger(const std::vector<su3_abelian> &array, int r,
                                 int time);
template std::vector<su3_abelian>
calculate_wilson_loops_schwinger_opposite(const std::vector<su3_abelian> &array,
                                          int r, int time);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_longitudinal_l(
    const std::vector<su3_abelian> &array,
    std::vector<std::vector<su3_abelian>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_longitudinal_tr(
    std::vector<su3_abelian> &array,
    std::vector<std::vector<su3_abelian>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_transversal_l(
    std::vector<su3_abelian> &array,
    std::vector<std::vector<su3_abelian>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_max);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_transversal_tr(
    std::vector<su3_abelian> &array,
    std::vector<std::vector<su3_abelian>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_max);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_longitudinal_l(
    std::vector<su3_abelian> &array,
    std::vector<std::vector<su3_abelian>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_longitudinal_tr(
    std::vector<su3_abelian> &array,
    std::vector<std::vector<su3_abelian>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_transversal_l(
    std::vector<su3_abelian> &array,
    std::vector<std::vector<su3_abelian>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_max);
template std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_transversal_tr(
    std::vector<su3_abelian> &array,
    std::vector<std::vector<su3_abelian>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_max);
template void
flux_schwinger_all(std::vector<su3_abelian> &array,
                   std::vector<std::vector<su3_abelian>> &schwinger_lines_short,
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
                   int T_min, int T_max, int R_min, int R_max, int d_ouside,
                   int d_max);
template std::vector<su3_abelian>
calculate_schwinger_lines_short(const std::vector<su3_abelian> &array, int d);
template std::vector<double>
calculate_plaket_time_trace_l(const std::vector<su3_abelian> &array);
template std::vector<double>
calculate_plaket_time_trace_tr(const std::vector<su3_abelian> &array);
template std::vector<double>
calculate_plaket_space_trace_l(const std::vector<su3_abelian> array);
template std::vector<double>
calculate_plaket_space_trace_tr(const std::vector<su3_abelian> &array);
template std::vector<double>
calculate_plaket_space_tr(const std::vector<su3_abelian> &array);
template std::vector<su3_abelian>
calculate_polyakov_loop(const std::vector<su3_abelian> &array);
template std::vector<su3_abelian>
calculate_wilson_loop(const std::vector<su3_abelian> &array, int r, int time);
template std::vector<double>
calculate_wilson_loop_tr(const std::vector<su3_abelian> &array, int r,
                         int time);
template std::map<std::tuple<int, int, int>, double>
calculate_wilson_plaket_correlator_electric_longitudinal(
    const std::vector<double> &plaket_tr,
    const std::vector<su3_abelian> &conf_wilson, int T_min, int T_max,
    int R_min, int R_max, int d_ouside);
template std::map<std::tuple<int, int, int>, double>
calculate_wilson_plaket_correlator_electric_transversal(
    const std::vector<double> &plaket_tr,
    const std::vector<su3_abelian> &conf_wilson, int T_min, int T_max,
    int R_min, int R_max, int d_ouside);

template std::vector<double>
wilson_plane_tr(std::vector<su3_abelian> &wilson_lines_mu,
                std::vector<su3_abelian> &wilson_lines_nu, int size_mu1,
                int size_mu2, int size_nu1, int size_nu2, int length_mu,
                int length_nu);

template std::vector<double> plaket_plane_tr(std::vector<su3_abelian> &conf_mu,
                                             std::vector<su3_abelian> &conf_nu,
                                             int size_mu1, int size_mu2,
                                             int size_nu1, int size_nu2);

template std::vector<double>
plaket_aver_tr_time(std::vector<std::vector<su3_abelian>> conf);
template std::vector<double>
plaket_aver_tr_space(std::vector<std::vector<su3_abelian>> conf);
template std::map<std::tuple<int, int, int>, double>
wilson_plaket_correlator(std::vector<double> &plaket_tr,
                         std::vector<std::vector<su3_abelian>> &conf_wilson,
                         int T_min, int T_max, int R_min, int R_max,
                         int main_coordinate, int transverse_coordinate,
                         std::string direction);
