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
#include "../include/link.h"

template <class T>
std::vector<T> calculate_plaket_time(const std::vector<T> array) {
  std::vector<T> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
  SPACE_ITER_START;
  for (int dir = 1; dir < 3; dir++) {
    link.move_dir(dir);
    vec.push_back(link.plaket_mu(array, 3));
  }
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<T> calculate_plaket_space(const std::vector<T> &array) {
  std::vector<T> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
  SPACE_ITER_START;
  for (int dir = 0; dir < 3; dir++) {
    for (int j = dir + 1; j < 3; j++) {
      link.move_dir(dir);
      vec[PLACE1_NODIR + dir + j - 2] = link.plaket_mu(array, j);
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
std::vector<T> calculate_plaket_schwinger_time(const std::vector<T> &array) {
  std::vector<T> vec(data_size / 4);
  link1 link(x_size, y_size, z_size, t_size);
  link.move_dir(3);
  SPACE_ITER_START;
  vec[link.place / 4] = (link.plaket_schwinger_average(array, 0) +
                         link.plaket_schwinger_average(array, 1) +
                         link.plaket_schwinger_average(array, 2)) *
                        (1. / 3);
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<std::vector<T>>
calculate_plaket_schwinger_space(const std::vector<T> &array) {
  std::vector<std::vector<T>> vec(3, std::vector<T>(data_size / 4));
  link1 link(x_size, y_size, z_size, t_size);
  int count = 0;
  for (int mu = 0; mu < 3; mu++) {
    for (int nu = mu + 1; nu < 3; nu++) {
      link.move_dir(mu);
      SPACE_ITER_START;
      vec[count][link.place / 4 * 3] = link.plaket_schwinger_average(array, nu);
      SPACE_ITER_END;
      count++;
    }
  }
  return vec;
}

template <class T>
std::vector<T> calculate_schwinger_lines_short(const std::vector<T> &array,
                                               int d) {
  std::vector<T> vec(data_size * 3 / 4);
  link1 link(x_size, y_size, z_size, t_size);
  for (int mu = 0; mu < 3; mu++) {
    link.move_dir(mu);
    SPACE_ITER_START;
    vec[link.place * 3 / 4 + mu] = link.wilson_line(array, d);
    SPACE_ITER_END;
  }
  return vec;
}

template <class T>
std::vector<std::vector<T>>
calculate_schwinger_line(const std::vector<T> &array, int d, int x_trans) {
  std::vector<std::vector<T>> vec(3, std::vector<T>(data_size));
  link1 link(x_size, y_size, z_size, t_size);
  for (int dir = 0; dir < 3; dir++) {
    link.move_dir(dir);
    SPACE_ITER_START;
    for (int nu = 0; nu < 3; nu++) {
      if (nu != dir)
        vec[nu - 1][PLACE3_DIR] = link.schwinger_line(array, d, nu, x_trans);
    }
    SPACE_ITER_END;
  }
  return vec;
}

template <class T>
std::vector<std::vector<T>>
calculate_wilson_loops_schwinger(const std::vector<T> &array, int r, int t) {
  std::vector<std::vector<T>> vec(3, std::vector<T>(data_size / 4));
  link1 link(x_size, y_size, z_size, t_size);
  for (int mu = 0; mu < 3; mu++) {
    link.move_dir(mu);
    SPACE_ITER_START;
    vec[mu][link.place / 4] = link.wilson_loop_schwinger(array, r, t);
    SPACE_ITER_END;
  }
  return vec;
}

template <class T>
std::map<int, double>
wilson_plaket_schwinger_electric(const std::vector<T> &array,
                                 const std::vector<T> &plaket, int d_min,
                                 int d_max, int t, int r) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<T> schwinger_lines;

  std::vector<std::vector<T>> wilson_loops =
      calculate_wilson_loops_schwinger(array, r, t);
  std::map<int, double> result;
  T A;
  for (int d = d_min; d <= d_max; d++) {

    schwinger_lines = calculate_schwinger_lines_short(array, abs(d));

    SPACE_ITER_START;

    for (int dir = 0; dir < 3; dir++) {
      A = wilson_loops[dir][link.place / 4];

      if (d < 0) {

        link.move(dir, d);
        A = (A ^ schwinger_lines[link.place * 3 / 4 + dir]);
        A = A * plaket[link.place / 4];
        A = A * schwinger_lines[link.place * 3 / 4 + dir];
        result[d] += A.tr();
        link.move(dir, -d);

      } else if (d > 0) {

        A = (A * schwinger_lines[link.place * 3 / 4 + dir]);
        link.move(dir, d);
        A = A * plaket[link.place / 4];
        link.move(dir, -d);
        A = A ^ schwinger_lines[link.place * 3 / 4 + dir];
        result[d] += A.tr();

      } else {

        A = A * plaket[link.place / 4];
        result[d] += A.tr();
      }
    }

    SPACE_ITER_END;
  }

  for (auto it = result.begin(); it != result.end(); ++it) {
    it->second = it->second / (x_size * y_size * z_size * t_size * 3);
  }

  return result;
}

template <class T>
std::vector<double> calculate_plaket_time_tr(const std::vector<T> &array) {
  std::vector<double> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
  SPACE_ITER_START;
  for (int dir = 0; dir < 3; dir++) {
    link.move_dir(dir);
    vec[link.place / 4 * 3 + dir] = link.plaket_mu(array, 3).tr();
  }
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<double> calculate_plaket_space_tr(const std::vector<T> &array) {
  std::vector<double> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
  int place_dir;
  SPACE_ITER_START;
  for (int dir = 0; dir < 3; dir++) {
    for (int j = dir + 1; j < 3; j++) {
      link.move_dir(dir);
      vec[link.place / 4 * 3 + dir + j] = link.plaket_mu(array, j).tr();
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
  for (int dir = 0; dir < 3; dir++) {
    link.move_dir(dir);
    SPACE_ITER_START;
    vec[link.place / 4 * 3 + dir] = link.wilson_loop(array, r, time).tr();
    SPACE_ITER_END;
  }
  return vec;
}

// calculate std::vector of wilson loop traces with particular time size
// template <class T>
// std::vector<std::vector<double>> calculate_wilson_loop_tr(const
// std::vector<T> &array,
//                                                std::vector<int>
//                                                space_sizes, int time_size)
//                                                {
//   link1 link(x_size, y_size, z_size, t_size);
//   std::vector<std::vector<double>> wilson(space_sizes.size());
//   std::vector<std::vector<T>> time_lines(time_max - time_min + 1);
//   std::vector<T> space_lines;
//   for (int i = time_min; i <= time_max; i++) {
//     time_lines[i - time_min] = wilson_lines(array, 3, i);
//   }
//   T A;
//   for (int dir = 0; dir < 3; dir++) {
//     for (int r = r_min; r <= r_max; r++) {
//       if (r == r_min)
//         space_lines = wilson_lines(array, dir, r);
//       else
//         space_lines = wilson_line_increase(array, space_lines, dir, r - 1);
//       for (int time = time_min; time <= time_max; time++) {

//         SPACE_ITER_START

//         A = time_lines[time - time_min][link.place / 4];
//         link.move(3, time);
//         A = A * space_lines[link.place / 4];
//         link.move(3, -time);
//         link.move(dir, r);
//         A = A * time_lines[time - time_min][link.place / 4].conj();
//         link.move(dir, -r);
//         A = A * space_lines[link.place / 4].conj();

//         wilson[(r - r_min) + (time - time_min) * (r_max - r_min + 1)] +=
//         A.tr();

//         SPACE_ITER_END
//       }
//     }
//   }
//   for (int i = 0; i < (time_max - time_min + 1) * (r_max - r_min + 1); i++)
//   {
//     wilson[i] = wilson[i] / (DATA_SIZE / 4 * 3);
//   }
//   return wilson;
// }

std::map<int, double>
wilson_plaket_correlator_electric(const std::vector<double> &wilson_loop_tr,
                                  const std::vector<double> &plaket_tr, int r,
                                  int time, int x_trans, int d_min, int d_max) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<double> correlator(d_max - d_min + 1, 0.0);
  double a;
  for (int dir = 0; dir < 3; dir++) {
    SPACE_ITER_START
    a = wilson_loop_tr[link.place / 4 * 3 + dir];
    link.move(3, time / 2);
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

// su2

template std::map<int, double>
wilson_plaket_schwinger_electric(const std::vector<su2> &array,
                                 const std::vector<su2> &plaket, int d_min,
                                 int d_max, int t, int r);
template std::vector<su2>
calculate_plaket_schwinger_time(const std::vector<su2> &array);
template std::vector<std::vector<su2>>
calculate_plaket_schwinger_space(const std::vector<su2> &array);
template std::vector<su2>
calculate_schwinger_lines_short(const std::vector<su2> &array, int d);
template std::vector<std::vector<su2>>
calculate_schwinger_line(const std::vector<su2> &array, int d, int x_trans);
template std::vector<double>
calculate_plaket_time_tr(const std::vector<su2> &array);
template std::vector<double>
calculate_plaket_space_tr(const std::vector<su2> &array);
template std::vector<su2>
calculate_polyakov_loop(const std::vector<su2> &array);
template std::vector<su2> calculate_wilson_loop(const std::vector<su2> &array,
                                                int r, int time);
template std::vector<double>
calculate_wilson_loop_tr(const std::vector<su2> &array, int r, int time);

// abelian

template std::map<int, double>
wilson_plaket_schwinger_electric(const std::vector<abelian> &array,
                                 const std::vector<abelian> &plaket, int d_min,
                                 int d_max, int t, int r);
template std::vector<abelian>
calculate_plaket_schwinger_time(const std::vector<abelian> &array);
template std::vector<std::vector<abelian>>
calculate_plaket_schwinger_space(const std::vector<abelian> &array);
template std::vector<std::vector<abelian>>
calculate_schwinger_line(const std::vector<abelian> &array, int d, int x_trans);
template std::vector<abelian>
calculate_schwinger_lines_short(const std::vector<abelian> &array, int d);
template std::vector<double>
calculate_plaket_time_tr(const std::vector<abelian> &array);
template std::vector<double>
calculate_plaket_space_tr(const std::vector<abelian> &array);
template std::vector<abelian>
calculate_polyakov_loop(const std::vector<abelian> &array);
template std::vector<abelian>
calculate_wilson_loop(const std::vector<abelian> &array, int r, int time);
template std::vector<double>
calculate_wilson_loop_tr(const std::vector<abelian> &array, int r, int time);
