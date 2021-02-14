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

template <class T> vector<T> calculate_plaket_time(const vector<T> array) {
  vector<T> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
  SPACE_ITER_START;
  for (int dir = 1; dir < 3; dir++) {
    link.move_dir(dir);
    vec.push_back(link.plaket_mu(array, 3));
  }
  SPACE_ITER_END;
  return vec;
}

template <class T> vector<T> calculate_plaket_space(const vector<T> &array) {
  vector<T> vec(data_size / 4 * 3);
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

template <class T> vector<T> calculate_polyakov_loop(const vector<T> &array) {
  vector<T> vec((data_size) / 4);
  link1 link(x_size, y_size, z_size, t_size);
  link.move_dir(3);
  SPACE_ITER_START;
  vec[PLACE1_NODIR] = link.polyakov_loop(array);
  SPACE_ITER_END;
  return vec;
}

template <class T>
vector<T> calculate_wilson_loop(const vector<T> &array, int r, int time) {
  vector<T> vec(data_size / 4 * 3);
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
vector<vector<T>> calculate_schwinger_line(const vector<T> &array, int d,
                                           int x_trans) {
  vector<vector<T>> vec(3, vector<T>(data_size));
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
void wilson_plaket_schwinger_electric(const vector<vector<T>> &schwinger_line,
                                      const vector<T> &plaket,
                                      const vector<T> &polyakov_loop,
                                      vector<vector<result>> &field1,
                                      vector<vector<result>> &field2,
                                      vector<result> &field3, int d, int D,
                                      int x_trans) {
  link1 link(x_size, y_size, z_size, t_size);
  for (int dir = 0; dir < 3; dir++) {
    T A, B;
    SPACE_ITER_START;
    A = polyakov_loop[PLACE1_LINK_NODIR];
    // B = schwinger_line;
    link.move(dir, d);
    if (x_trans == 0) {
      link.move(dir, d);
      A = A * B * plaket[PLACE3_DIR] * B.conj();
    }
    // for(int mu = 1;mu < 4;mu++){
    //   if(mu != dir){
    //     link.move()
    //   }
    // }
    SPACE_ITER_END;
  }
}

template <class T>
vector<FLOAT> calculate_plaket_time_tr(const vector<T> &array) {
  vector<FLOAT> vec(data_size / 4 * 3);
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
vector<FLOAT> calculate_plaket_space_tr(const vector<T> &array) {
  vector<FLOAT> vec(data_size / 4 * 3);
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

FLOAT plaket4_time(const vector<FLOAT> &plaket_tr, link1 &link) {
  FLOAT a = plaket_tr[link.place / 4 * 3 + link.direction];
  link.move(link.direction, -1);
  a += plaket_tr[link.place / 4 * 3 + link.direction];
  link.move(3, -1);
  a += plaket_tr[link.place / 4 * 3 + link.direction];
  link.move(link.direction, 1);
  a += plaket_tr[link.place / 4 * 3 + link.direction];
  link.move(3, 1);
  return a / 4;
}

FLOAT plaket4_space(const vector<FLOAT> &plaket_tr, link1 &link, int nu) {
  FLOAT a = plaket_tr[link.place / 4 * 3 + link.direction + nu];
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
vector<FLOAT> calculate_wilson_loop_tr(const vector<T> &array, int r,
                                       int time) {
  vector<FLOAT> vec(data_size / 4 * 3);
  link1 link(x_size, y_size, z_size, t_size);
  for (int dir = 0; dir < 3; dir++) {
    link.move_dir(dir);
    SPACE_ITER_START;
    vec[link.place / 4 * 3 + dir] = link.wilson_loop(array, r, time).tr();
    SPACE_ITER_END;
  }
  return vec;
}

result wilson_plaket_correlator_electric(const vector<FLOAT> &wilson_loop_tr,
                                         const vector<FLOAT> &plaket_tr, int r,
                                         int time, int x_trans, int d_min,
                                         int d_max) {
  link1 link(x_size, y_size, z_size, t_size);
  FLOAT vec[d_max - d_min + 1];
  result final(0);
  FLOAT aver[2];
  FLOAT a;
  for (int dir = 0; dir < 3; dir++) {
    SPACE_ITER_START
    a = wilson_loop_tr[link.place / 4 * 3 + dir];
    link.move(3, time / 2);
    link.move(dir, d_min);
    for (int d = d_min; d <= d_max; d++) {
      if (x_trans == 0) {
        for (int mu = 0; mu < 3; mu++) {
          link.move_dir(mu);
          vec[d - d_min] += a * plaket4_time(plaket_tr, link);
        }
      } else {
        for (int nu = 0; nu < 3; nu++) {
          if (nu != dir) {
            link.move(nu, x_trans);
            for (int mu = 0; mu < 3; mu++) {
              link.move_dir(mu);
              vec[d - d_min] += a * plaket4_time(plaket_tr, link);
            }
            link.move(nu, -2 * x_trans);
            for (int mu = 0; mu < 3; mu++) {
              link.move_dir(mu);
              vec[d - d_min] += a * plaket4_time(plaket_tr, link);
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
  for (int d = d_min; d <= d_max; d++) {
    final.array.push_back(vec[d - d_min] / count);
  }
  return final;
}

result wilson_plaket_correlator_electric_x(const vector<FLOAT> &wilson_loop_tr,
                                           const vector<FLOAT> &plaket_tr,
                                           int r, int time, int x_trans_min,
                                           int x_trans_max, int d) {
  link1 link(x_size, y_size, z_size, t_size);
  FLOAT vec[x_trans_max - x_trans_min + 1];
  result final(0);
  FLOAT aver[2];
  FLOAT a;
  for (int dir = 0; dir < 3; dir++) {
    SPACE_ITER_START
    a = wilson_loop_tr[link.place / 4 * 3 + dir];
    link.move(3, time / 2);
    link.move(dir, d);
    for (int nu = 0; nu < 3; nu++) {
      if (nu != dir) {
        link.move(nu, x_trans_min - 1);
        for (int x_trans = x_trans_min; x_trans <= x_trans_max; x_trans++) {
          link.move(nu, 1);
          for (int mu = 0; mu < 3; mu++) {
            link.move_dir(mu);
            vec[x_trans - x_trans_min] += a * plaket4_time(plaket_tr, link);
          }
        }
        link.move(nu, -x_trans_max);
      }
    }
    SPACE_ITER_END
  }
  int count;
  count = data_size / 4 * 18;
  for (int x_trans = x_trans_min; x_trans <= x_trans_max; x_trans++) {
    final.array.push_back(vec[x_trans - x_trans_min] / count);
  }
  return final;
}

result wilson_plaket_correlator_magnetic(const vector<FLOAT> &wilson_loop_tr,
                                         const vector<FLOAT> &plaket_tr, int r,
                                         int time, int x_trans, int d_min,
                                         int d_max) {
  link1 link(x_size, y_size, z_size, t_size);
  FLOAT vec[d_max - d_min + 1];
  result final(0);
  FLOAT aver[2];
  FLOAT a;
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
            vec[d - d_min] += a * plaket4_space(plaket_tr, link, j);
          }
        }
      } else {
        for (int nu = 0; nu < 3; nu++) {
          if (nu != dir) {
            link.move(nu, x_trans);
            for (int mu = 0; mu < 3; mu++) {
              for (int j = mu + 1; j < 3; j++) {
                link.move_dir(mu);
                vec[d - d_min] += a * plaket4_space(plaket_tr, link, j);
              }
            }
            link.move(nu, -2 * x_trans);
            for (int mu = 0; mu < 3; mu++) {
              for (int j = mu + 1; j < 3; j++) {
                link.move_dir(mu);
                vec[d - d_min] += a * plaket4_space(plaket_tr, link, j);
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
  for (int d = d_min; d <= d_max; d++) {
    final.array.push_back(vec[d - d_min] / count);
  }
  return final;
}

result wilson_plaket_correlator_magnetic_x(const vector<FLOAT> &wilson_loop_tr,
                                           const vector<FLOAT> &plaket_tr,
                                           int R, int T, int x_trans_min,
                                           int x_trans_max, int d) {
  link1 link(x_size, y_size, z_size, t_size);
  FLOAT vec[x_trans_max - x_trans_min + 1];
  result final(0);
  FLOAT aver[2];
  FLOAT a;
  for (int dir = 0; dir < 3; dir++) {
    SPACE_ITER_START
    a = wilson_loop_tr[link.place / 4 * 3 + dir];
    link.move(3, T / 2);
    link.move(dir, d);
    for (int nu = 0; nu < 3; nu++) {
      if (nu != dir) {
        link.move(nu, x_trans_min - 1);
        for (int x_trans = x_trans_min; x_trans <= x_trans_max; x_trans++) {
          link.move(nu, 1);
          for (int mu = 0; mu < 3; mu++) {
            for (int j = mu + 1; j < 3; j++) {
              link.move_dir(mu);
              vec[x_trans - x_trans_min] +=
                  a * plaket4_space(plaket_tr, link, j);
            }
          }
        }
        link.move(nu, -x_trans_max);
      }
    }
    SPACE_ITER_END
  }
  int count;
  count = data_size / 4 * 18;
  for (int x_trans = x_trans_min; x_trans <= x_trans_max; x_trans++) {
    final.array.push_back(vec[x_trans - x_trans_min] / count);
  }
  return final;
}

// su2

template void wilson_plaket_schwinger_electric(
    const vector<vector<su2>> &schwinger_line, const vector<su2> &plaket,
    const vector<su2> &polyakov_loop, vector<vector<result>> &field1,
    vector<vector<result>> &field2, vector<result> &field3, int d, int D,
    int x_trans);
template vector<vector<su2>> calculate_schwinger_line(const vector<su2> &array,
                                                      int d, int x_trans);
template vector<FLOAT> calculate_plaket_time_tr(const vector<su2> &array);
template vector<FLOAT> calculate_plaket_space_tr(const vector<su2> &array);
template vector<su2> calculate_polyakov_loop(const vector<su2> &array);
template vector<su2> calculate_wilson_loop(const vector<su2> &array, int r,
                                           int time);
template vector<FLOAT> calculate_wilson_loop_tr(const vector<su2> &array, int r,
                                                int time);

// abelian

template void wilson_plaket_schwinger_electric(
    const vector<vector<abelian>> &schwinger_line,
    const vector<abelian> &plaket, const vector<abelian> &polyakov_loop,
    vector<vector<result>> &field1, vector<vector<result>> &field2,
    vector<result> &field3, int d, int D, int x_trans);
template vector<vector<abelian>>
calculate_schwinger_line(const vector<abelian> &array, int d, int x_trans);
template vector<FLOAT> calculate_plaket_time_tr(const vector<abelian> &array);
template vector<FLOAT> calculate_plaket_space_tr(const vector<abelian> &array);
template vector<abelian> calculate_polyakov_loop(const vector<abelian> &array);
template vector<abelian> calculate_wilson_loop(const vector<abelian> &array,
                                               int r, int time);
template vector<FLOAT> calculate_wilson_loop_tr(const vector<abelian> &array,
                                                int r, int time);