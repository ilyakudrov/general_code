#include "../include/link.h"
#include "../include/data.h"
#include <cmath>

#define Pi 3.141592653589793238462643383279502884
#define data_size                                                              \
  4 * lattice_size[0] * lattice_size[1] * lattice_size[2] * lattice_size[3]
#define PLACE3_LINK_DIR                                                        \
  (coordinate[3]) * 3 * lattice_size[0] * lattice_size[1] * lattice_size[2] +  \
      (coordinate[2]) * 3 * lattice_size[0] * lattice_size[1] +                \
      (coordinate[1]) * 3 * lattice_size[0] + (coordinate[0]) * 3 +            \
      abs(direction) - 1
#define PLACE4_LINK_DIR                                                        \
  (coordinate[3]) * 4 * lattice_size[0] * lattice_size[1] * lattice_size[2] +  \
      (coordinate[2]) * 4 * lattice_size[0] * lattice_size[1] +                \
      (coordinate[1]) * 4 * lattice_size[0] + (coordinate[0]) * 4 +            \
      abs(direction) - 1
#define PLACE1_LINK_NODIR                                                      \
  (coordinate[3]) * lattice_size[0] * lattice_size[1] * lattice_size[2] +      \
      (coordinate[2]) * lattice_size[0] * lattice_size[1] +                    \
      (coordinate[1]) * lattice_size[0] + (coordinate[0])

template <class T>
link1<T>::link1(int lattice_size_x, int lattice_size_y, int lattice_size_z,
                int lattice_size_t) {
  lattice_size[0] = lattice_size_x;
  lattice_size[1] = lattice_size_y;
  lattice_size[2] = lattice_size_z;
  lattice_size[3] = lattice_size_t;
  direction = 1;
  for (int i = 0; i < 4; i++) {
    coordinate[i] = 0;
  }
}

template <class T> ostream &operator<<(ostream &os, const link1<T> &link) {
  os << "x: " << link.coordinate[0] << " y: " << link.coordinate[1]
     << " z: " << link.coordinate[2] << " t: " << link.coordinate[3]
     << " dir: " << link.direction << endl;
  return os;
}

template <class T> void link1<T>::move(int dir, int step) {
  coordinate[abs(dir) - 1] +=
      (step * (dir / abs(dir)) + lattice_size[abs(dir) - 1]);
  coordinate[abs(dir) - 1] =
      coordinate[abs(dir) - 1] % lattice_size[abs(dir) - 1];
}
template <class T> void link1<T>::go(int x, int y, int z, int t) {
  coordinate[0] = x;
  coordinate[1] = y;
  coordinate[2] = z;
  coordinate[3] = t;
}
template <class T> void link1<T>::move_dir(int dir) { direction = dir; }
template <class T> int link1<T>::get_place() { return PLACE4_LINK_DIR; }
template <class T> int link1<T>::get_place1() { return PLACE1_LINK_NODIR; }

template <class T> T link1<T>::get_matrix(const vector<T> &vec) {
  T A;
  if (direction / abs(direction) == 1) {
    return vec[PLACE4_LINK_DIR];
  }
  if (direction / abs(direction) == -1) {
    move(direction, 1);
    A = vec[PLACE4_LINK_DIR];
    move(direction, -1);
    return A.conj();
  } else {
    cout << "error in get_matrix" << endl;
    return A;
  }
}

template <class T> FLOAT link1<T>::border_sign(int mu) {
  if (abs(mu) == 4) {
    if (0 > (coordinate[3] + mu / abs(mu)) ||
        (coordinate[3] + mu / abs(mu)) > lattice_size[3] - 1)
      return -1.;
    else
      return 1.;
  } else
    return 1.;
}

template <> FLOAT link1<su2>::get_angle_abelian(const vector<su2> &vec) {
  FLOAT angle;
  if (direction / abs(direction) == 1) {
    angle = atan2(vec[PLACE4_LINK_DIR].a3, vec[PLACE4_LINK_DIR].a0);
    return angle;
  }
  if (direction / abs(direction) == -1) {
    move(direction, 1);
    angle = atan2(vec[PLACE4_LINK_DIR].a3, vec[PLACE4_LINK_DIR].a0);
    move(direction, -1);
    return -angle;
  } else {
    cout << "error in get_matrix" << endl;
    return -1;
  }
}

template <class T>
T link1<T>::schwinger_line(const vector<T> &array, int d, int dir, int x) {
  int dir1 = direction;
  T A;
  for (int i = 0; i < d; i++) {
    A = A * get_matrix(array);
    move(dir1, 1);
  }
  move_dir(dir);
  for (int i = 0; i < x; i++) {
    A = A * get_matrix(array);
    move(dir, 1);
  }
  move_dir(dir1);
  move(-dir, x);
  move(-dir1, d);
  return A;
}

template <class T> T link1<T>::plaket_mu(const vector<T> &array, int mu) {
  int dir = direction;
  T A = get_matrix(array);
  move(dir, 1);
  move_dir(-mu); // the same temporal direction as polyakov loop, connected to
                 // schwinger line, has
  A = A * get_matrix(array);
  move(-mu, 1);
  move_dir(-dir);
  A = A * get_matrix(array);
  move(-dir, 1);
  move_dir(mu);
  A = A * get_matrix(array);
  move(mu, 1);
  move_dir(dir);
  return A;
}

template <class T> T link1<T>::plaket_average4(const vector<T> &array, int mu) {
  T A = plaket_mu(array, mu);
  move(mu, 1);
  A = A + plaket_mu(array, mu);
  move(-direction, 1);
  A = A + plaket_mu(array, mu);
  move(mu, -1);
  A = A + plaket_mu(array, mu);
  move(direction, 1);
  return A;
}

template <class T> T link1<T>::polyakov_loop(const vector<T> &array) {
  T A;
  for (int i = 0; i < lattice_size[3]; i++) {
    A = A * get_matrix(array);
    move(direction, 1);
  }
  return A;
}

template <class T>
T link1<T>::wilson_loop(const vector<T> &array, int r, int t) {
  int dir = direction;
  T A;
  for (int i = 0; i < r; i++) {
    A = A * get_matrix(array);
    move(dir, 1);
  }
  move_dir(4);
  for (int i = 0; i < t; i++) {
    A = A * get_matrix(array);
    move(4, 1);
  }
  move_dir(-dir);
  for (int i = 0; i < r; i++) {
    A = A * get_matrix(array);
    move(-dir, 1);
  }
  move_dir(-4);
  for (int i = 0; i < t; i++) {
    A = A * get_matrix(array);
    move(-4, 1);
  }
  move_dir(dir);
  return A;
}

template <class T> T link1<T>::wilson_line(const vector<T> &array, int length) {
  int dir = direction;
  T A;
  for (int i = 0; i < length; i++) {
    A = A * get_matrix(array);
    move(dir, 1);
  }
  return A;
}

template <class T>
FLOAT link1<T>::field1(const vector<vector<T>> &schwinger_line,
                       const vector<T> &plaket, const vector<T> &polyakov_loop,
                       int d, int D, int dir, int x) {
  int dir1 = direction;
  T C = schwinger_line[dir - 1][PLACE3_LINK_DIR];
  move(dir1, d);
  move(dir, x);
  T B = plaket[PLACE3_LINK_DIR];
  move(-dir, x);
  move(-dir1, d);
  move_dir(-4);
  T A = polyakov_loop[PLACE1_LINK_NODIR];
  A = A.conj() * C * B;
  A = A * C.conj();
  move(dir1, D);
  move_dir(4);
  B = polyakov_loop[PLACE1_LINK_NODIR];
  move(-dir1, D);
  move_dir(dir1);
  return A.tr() * B.tr();
}

template <class T>
FLOAT link1<T>::field2(const vector<T> &plaket, const vector<T> &polyakov_loop,
                       int d, int D, int dir, int x) {
  int dir1 = direction;
  move_dir(-4);
  T A = polyakov_loop[PLACE1_LINK_NODIR];
  move(dir1, d);
  move(dir, x);
  move_dir(dir1);
  T B = plaket[PLACE3_LINK_DIR];
  move(-dir, x);
  move(dir1, D - d);
  move_dir(4);
  T C = polyakov_loop[PLACE1_LINK_NODIR];
  move(-dir1, D);
  move_dir(dir1);
  return B.tr() * C.tr() * A.conj().tr();
}

template <class T>
FLOAT link1<T>::field3(const vector<T> &polyakov_loop, int D, int x) {
  int dir1 = direction;
  move_dir(-4);
  T A = polyakov_loop[PLACE1_LINK_NODIR];
  move(dir1, D);
  move_dir(4);
  T B = polyakov_loop[PLACE1_LINK_NODIR];
  move(-dir1, D);
  move_dir(dir1);
  return A.conj().tr() * B.tr();
}

template <class T> T link1<T>::staples_first(const vector<T> &vec, int eta) {
  T A;
  T B;
  int dir = direction;
  move_dir(eta);
  A = get_matrix(vec);
  move(eta, 1);
  move_dir(dir);
  A = A * get_matrix(vec);
  move(dir, 1);
  move_dir(-eta);
  A = A * get_matrix(vec);
  move(-eta, 1);
  move(-dir, 1);
  B = get_matrix(vec);
  move(-eta, 1);
  move_dir(dir);
  B = B * get_matrix(vec);
  move(dir, 1);
  move_dir(eta);
  B = B * get_matrix(vec);
  move(eta, 1);
  move(-dir, 1);
  move_dir(dir);
  return (A + B);
}

template <class T>
T link1<T>::staples_second(const vector<vector<T>> &smearing_first,
                           unordered_map<int, int> &indexes, int rho, int mu,
                           int nu) {
  T A;
  T B;
  int a = 0;
  A = smearing_first[indexes[rho * 100 + nu * 10 + mu]][PLACE1_LINK_NODIR];
  move(rho, 1);
  A = A * smearing_first[indexes[mu * 100 + rho * 10 + nu]][PLACE1_LINK_NODIR];
  move(direction, 1);
  move(-rho, 1);
  A = A * smearing_first[indexes[rho * 100 + nu * 10 + mu]][PLACE1_LINK_NODIR]
              .conj();
  move(-direction, 1);
  move(-rho, 1);
  B = smearing_first[indexes[rho * 100 + nu * 10 + mu]][PLACE1_LINK_NODIR]
          .conj();
  B = B * smearing_first[indexes[mu * 100 + rho * 10 + nu]][PLACE1_LINK_NODIR];
  move(direction, 1);
  B = B * smearing_first[indexes[rho * 100 + nu * 10 + mu]][PLACE1_LINK_NODIR];
  move(rho, 1);
  move(-direction, 1);
  return (A + B);
}

template <class T>
T link1<T>::staples_second_refresh(const vector<T> &vec, int eta, int nu,
                                   FLOAT alpha3) {
  T A;
  T B;
  int dir = direction;
  move_dir(eta);
  A = smearing_first_refresh(vec, nu, dir, alpha3);
  move(eta, 1);
  move_dir(dir);
  A = A * smearing_first_refresh(vec, nu, eta, alpha3);
  move(dir, 1);
  move_dir(-eta);
  A = A * smearing_first_refresh(vec, nu, dir, alpha3);
  move(-eta, 1);
  move(-dir, 1);
  B = smearing_first_refresh(vec, nu, dir, alpha3);
  move(-eta, 1);
  move_dir(dir);
  B = B * smearing_first_refresh(vec, nu, eta, alpha3);
  move(dir, 1);
  move_dir(eta);
  B = B * smearing_first_refresh(vec, nu, dir, alpha3);
  move(eta, 1);
  move(-dir, 1);
  move_dir(dir);
  return (A + B);
}

template <class T>
T link1<T>::staples_third(const vector<vector<T>> &smearing_second,
                          unordered_map<int, int> indexes, int nu, int mu) {
  T A;
  T B;
  A = smearing_second[indexes[nu * 10 + mu]][PLACE1_LINK_NODIR];
  move(nu, 1);
  A = A * smearing_second[indexes[mu * 10 + nu]][PLACE1_LINK_NODIR];
  move(direction, 1);
  move(-nu, 1);
  A = A * smearing_second[indexes[nu * 10 + mu]][PLACE1_LINK_NODIR].conj();
  move(-direction, 1);
  move(-nu, 1);
  B = smearing_second[indexes[nu * 10 + mu]][PLACE1_LINK_NODIR].conj();
  B = B * smearing_second[indexes[mu * 10 + nu]][PLACE1_LINK_NODIR];
  move(direction, 1);
  B = B * smearing_second[indexes[nu * 10 + mu]][PLACE1_LINK_NODIR];
  move(nu, 1);
  move(-direction, 1);
  return (A + B);
}

template <class T>
T link1<T>::staples_third_refresh(const vector<T> &vec, int eta, FLOAT alpha2,
                                  FLOAT alpha3) {
  T A;
  T B;
  int dir = direction;
  move_dir(eta);
  A = smearing_second_refresh(vec, dir, alpha2, alpha3);
  move(eta, 1);
  move_dir(dir);
  A = A * smearing_second_refresh(vec, eta, alpha2, alpha3);
  move(dir, 1);
  move_dir(-eta);
  A = A * smearing_second_refresh(vec, dir, alpha2, alpha3);
  move(-eta, 1);
  move(-dir, 1);
  B = smearing_second_refresh(vec, dir, alpha2, alpha3);
  move(-eta, 1);
  move_dir(dir);
  B = B * smearing_second_refresh(vec, eta, alpha2, alpha3);
  move(dir, 1);
  move_dir(eta);
  B = B * smearing_second_refresh(vec, dir, alpha2, alpha3);
  move(eta, 1);
  move(-dir, 1);
  move_dir(dir);
  return (A + B);
}

template <class T>
vector<T> link1<T>::smearing_first(const vector<T> &array, FLOAT alpha3, int mu,
                                   int nu, int rho) {
  vector<T> vec(data_size / 4);
  move_dir(mu);
  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
          go(x, y, z, t);
          vec[PLACE1_LINK_NODIR] = (1 - alpha3) * get_matrix(array);
          for (int i = 1; i < 5; i++) {
            if (i != mu && i != nu && i != rho) {
              vec[PLACE1_LINK_NODIR] = vec[PLACE1_LINK_NODIR] +
                                       alpha3 / 2. * staples_first(array, i);
            }
          }
          vec[PLACE1_LINK_NODIR] = vec[PLACE1_LINK_NODIR].proj();
        }
      }
    }
  }
  return vec;
}

void make_map_first(unordered_map<int, int> &indexes) {
  int key;
  int count = 0;
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      int k = 4;
      if (i != j) {
        key = i * 100 + j * 10 + k;
        indexes[key] = count;
        key = i * 100 + k * 10 + j;
        indexes[key] = count;
        count++;
      }
    }
  }
  int i = 4;
  for (int j = 1; j <= 2; j++) {
    for (int k = j + 1; k <= 3; k++) {
      key = i * 100 + j * 10 + k;
      indexes[key] = count;
      key = i * 100 + k * 10 + j;
      indexes[key] = count;
      count++;
    }
  }
}

template <class T>
vector<vector<T>> link1<T>::smearing_first_full(const vector<T> &array,
                                                FLOAT alpha3) {
  unordered_map<int, int> indexes;
  make_map_first(indexes);
  vector<vector<T>> smearing(9, vector<T>(data_size / 4));
  int key;
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      int k = 4;
      if (i != j) {
        key = i * 100 + j * 10 + k;
        smearing[indexes[key]] = smearing_first(array, alpha3, i, j, k);
      }
    }
  }
  int i = 4;
  for (int j = 1; j <= 2; j++) {
    for (int k = j + 1; k <= 3; k++) {
      key = i * 100 + j * 10 + k;
      smearing[indexes[key]] = smearing_first(array, alpha3, i, j, k);
    }
  }
  return smearing;
}

template <class T>
vector<T> link1<T>::smearing_second(const vector<T> &array,
                                    vector<vector<T>> &smearing_first,
                                    FLOAT alpha2, int mu, int nu) {
  unordered_map<int, int> indexes;
  make_map_first(indexes);
  vector<T> vec(data_size / 4);
  move_dir(mu);
  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
          go(x, y, z, t);
          vec[PLACE1_LINK_NODIR] = (1 - alpha2) * get_matrix(array);
          for (int i = 1; i < 5; i++) {
            if (i != mu && i != nu) {
              vec[PLACE1_LINK_NODIR] =
                  vec[PLACE1_LINK_NODIR] +
                  alpha2 / 4. *
                      staples_second(smearing_first, indexes, i, mu, nu);
            }
          }
          vec[PLACE1_LINK_NODIR] = vec[PLACE1_LINK_NODIR].proj();
        }
      }
    }
  }
  return vec;
}

void make_map_second(unordered_map<int, int> &indexes) {
  int key;
  int count = 0;
  for (int i = 1; i <= 3; i++) {
    key = i * 10 + 4;
    indexes[key] = count;
    count++;
  }
  for (int i = 1; i <= 3; i++) {
    key = 40 + i;
    indexes[key] = count;
    count++;
  }
}

template <class T>
vector<vector<T>> link1<T>::smearing_second_full(
    const vector<T> &array, vector<vector<T>> &smearing_first, FLOAT alpha2) {
  unordered_map<int, int> indexes;
  make_map_second(indexes);
  vector<vector<T>> smearing(6, vector<T>(data_size / 4));
  int key;
  for (int i = 1; i <= 3; i++) {
    key = i * 10 + 4;
    smearing[indexes[key]] =
        smearing_second(array, smearing_first, alpha2, i, 4);
  }
  for (int i = 1; i <= 3; i++) {
    key = 40 + i;
    smearing[indexes[key]] =
        smearing_second(array, smearing_first, alpha2, 4, i);
  }
  return smearing;
}

template <class T>
vector<T> link1<T>::smearing_HYP(const vector<T> &array,
                                 vector<vector<T>> &smearing_second,
                                 FLOAT alpha1) {
  unordered_map<int, int> indexes;
  make_map_second(indexes);
  vector<T> vec(data_size);
  move_dir(4);
  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
          go(x, y, z, t);
          vec[PLACE4_LINK_DIR] = (1 - alpha1) * get_matrix(array);
          for (int i = 1; i < 4; i++) {
            vec[PLACE4_LINK_DIR] =
                vec[PLACE4_LINK_DIR] +
                alpha1 / 6. * staples_third(smearing_second, indexes, i, 4);
          }
          vec[PLACE4_LINK_DIR] = vec[PLACE4_LINK_DIR].proj();
        }
      }
    }
  }
  for (int d = 1; d < 4; d++) {
    move_dir(d);
    for (int t = 0; t < t_size; t++) {
      for (int z = 0; z < z_size; z++) {
        for (int y = 0; y < y_size; y++) {
          for (int x = 0; x < x_size; x++) {
            go(x, y, z, t);
            vec[PLACE4_LINK_DIR] = array[PLACE4_LINK_DIR];
          }
        }
      }
    }
  }
  return vec;
}

template <class T>
vector<T> link1<T>::smearing_APE(const vector<T> &array, FLOAT alpha_APE) {
  vector<T> vec(data_size);
  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
          go(x, y, z, t);
          for (int i = 1; i < 4; i++) {
            move_dir(i);
            vec[PLACE4_LINK_DIR] = (1 - alpha_APE) * array[PLACE4_LINK_DIR];
            for (int d = 1; d < 4; d++) {
              if (d != i) {
                vec[PLACE4_LINK_DIR] =
                    vec[PLACE4_LINK_DIR] +
                    (alpha_APE / 6.) * staples_first(array, d);
              }
            }
            vec[PLACE4_LINK_DIR] = vec[PLACE4_LINK_DIR].proj();
          }
        }
      }
    }
  }
  move_dir(4);
  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
          go(x, y, z, t);
          vec[PLACE4_LINK_DIR] = array[PLACE4_LINK_DIR];
        }
      }
    }
  }
  return vec;
}

template <class T>
T link1<T>::smearing_first_refresh(const vector<T> &vec, int nu, int rho,
                                   FLOAT alpha3) {
  T A;
  A = (1 - alpha3) * get_matrix(vec);
  for (int d = 1; d < 5; d++) {
    if (d != abs(nu) && d != abs(rho) && d != abs(direction)) {
      A = A + alpha3 / 2. * staples_first(vec, d);
    }
  }
  A.proj();
  return A;
}

template <class T>
T link1<T>::smearing_second_refresh(const vector<T> &vec, int nu, FLOAT alpha2,
                                    FLOAT alpha3) {
  T A;
  A = (1 - alpha2) * get_matrix(vec);
  for (int d = 1; d < 5; d++) {
    if (d != abs(nu) && d != abs(direction)) {
      A = A + alpha2 / 4. * staples_second_refresh(vec, d, nu, alpha3);
    }
  }
  A.proj();
  return A;
}

template <class T>
vector<T> link1<T>::smearing_HYP_refresh(data<T> &conf, FLOAT alpha1,
                                         FLOAT alpha2, FLOAT alpha3) {
  vector<T> vec(data_size);
  vec = conf.array;
  T A;
  move_dir(4);
  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
          go(x, y, z, t);
          A = (1 - alpha1) * get_matrix(vec);
          for (int d = 1; d < 4; d++) {
            A = A + alpha1 / 6. * staples_third_refresh(vec, d, alpha2, alpha3);
          }
          vec[PLACE4_LINK_DIR] = A.proj();
        }
      }
    }
  }
  return vec;
}

template <class T>
vector<T> link1<T>::smearing_APE_refresh(data<T> &conf, FLOAT alpha_APE) {
  vector<T> vec(data_size);
  vec = conf.array;
  T A;
  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
          go(x, y, z, t);
          for (int i = 1; i < 4; i++) {
            move_dir(i);
            A = (1 - alpha_APE) * conf.array[PLACE4_LINK_DIR];
            for (int d = 1; d < 4; d++) {
              if (d != i)
                A = A + (alpha_APE / 4.) * staples_first(vec, d);
            }
            vec[PLACE4_LINK_DIR] = A.proj();
          }
        }
      }
    }
  }
  return vec;
}

template <> su2 link1<su2>::stout_omega(data<su2> &conf, FLOAT rho);
template <> su2 link1<su2>::stout_factor(data<su2> &conf, FLOAT rho);

template <> vector<su2> link1<su2>::smearing_stout(data<su2> &conf, FLOAT rho) {
  vector<su2> vec(data_size);
  vec = conf.array;
  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
          go(x, y, z, t);
          for (int i = 1; i < 5; i++) {
            move_dir(i);
            vec[PLACE4_LINK_DIR] =
                stout_factor(conf, rho) * conf.array[PLACE4_LINK_DIR];
          }
        }
      }
    }
  }
  return vec;
}

template <> su2 link1<su2>::stout_factor(data<su2> &conf, FLOAT rho) {
  su2 A;
  su2 B;
  su2 C;
  su2 C1;
  C = stout_omega(conf, rho);
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

template <> su2 link1<su2>::stout_omega(data<su2> &conf, FLOAT rho) {
  int dir = direction;
  su2 A;
  su2 B(0., 0., 0., 0.);
  A = conf.array[PLACE4_LINK_DIR].inverse();
  for (int i = 1; i < 5; i++) {
    if (i != dir) {
      B = B + staples_first(conf.array, i);
    }
  }
  A = B * A * rho;
  return A;
}

template <class T> void link1<T>::gauge_transform(data<T> &conf) {
  T A = conf.array[0];
  T C = A.conj();
  FLOAT a;
  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
          for (int i = 1; i < 5; i++) {
            a = PLACE4_LINK_DIR;
            conf.array[a] = C * conf.array[a] * A;
          }
        }
      }
    }
  }
}

// monopoles

/*template <class T>
FLOAT link1<T>::monopole_plaket(data<T> &conf, int i, int j) {
  int dir = direction;
  move_dir(i);
  float angle = get_angle_abelian(conf.array);
  move(direction, 1);
  move_dir(j);
  angle += get_angle_abelian(conf.array);
  move(direction, 1);
  move_dir(-i);
  angle += get_angle_abelian(conf.array);
  move(direction, 1);
  move_dir(-j);
  angle += get_angle_abelian(conf.array);
  move(direction, 1);
  move_dir(dir);
  while ((angle > Pi) || (angle < -Pi)) {
    if (angle > Pi)
      angle = angle - 2 * Pi;
    if (angle < -Pi)
      angle = angle + 2 * Pi;
  }
  return angle;
}

template <class T> FLOAT link1<T>::get_current(data<T> &conf) {
  FLOAT jj = 0.;
  link1 Lx(coordinate[0], coordinate[1], coordinate[2], coordinate[3],
           lattice_size[0], lattice_size[1], lattice_size[2], lattice_size[3]);
  Lx.move(1, 1);
  link1 Ly(coordinate[0], coordinate[1], coordinate[2], coordinate[3],
           lattice_size[0], lattice_size[1], lattice_size[2], lattice_size[3]);
  Ly.move(2, 1);
  link1 Lz(coordinate[0], coordinate[1], coordinate[2], coordinate[3],
           lattice_size[0], lattice_size[1], lattice_size[2], lattice_size[3]);
  Lz.move(3, 1);
  link1 Lt(coordinate[0], coordinate[1], coordinate[2], coordinate[3],
           lattice_size[0], lattice_size[1], lattice_size[2], lattice_size[3]);
  Lt.move(4, 1);
  if (direction == 4)
    jj = Lx.monopole_plaket(conf, 2, 3) - monopole_plaket(conf, 2, 3) -
         (Ly.monopole_plaket(conf, 1, 3) - monopole_plaket(conf, 1, 3)) +
         Lz.monopole_plaket(conf, 1, 2) - monopole_plaket(conf, 1, 2);
  if (direction == 1)
    jj = -(Lt.monopole_plaket(conf, 2, 3) - monopole_plaket(conf, 2, 3)) +
         (Ly.monopole_plaket(conf, 4, 3) - monopole_plaket(conf, 4, 3)) -
         (Lz.monopole_plaket(conf, 4, 2) - monopole_plaket(conf, 4, 2));
  if (direction == 2)
    jj = Lt.monopole_plaket(conf, 1, 3) - monopole_plaket(conf, 1, 3) -
         (Lx.monopole_plaket(conf, 4, 3) - monopole_plaket(conf, 4, 3)) +
         Lz.monopole_plaket(conf, 4, 1) - monopole_plaket(conf, 4, 1);
  if (direction == 3)
    jj = -(Lt.monopole_plaket(conf, 1, 2) - monopole_plaket(conf, 1, 2)) +
         (Lx.monopole_plaket(conf, 4, 2) - monopole_plaket(conf, 4, 2)) -
         (Ly.monopole_plaket(conf, 4, 1) - monopole_plaket(conf, 4, 1));
  return jj / 2. / Pi;
}

template <class T> int link1<T>::current_test(FLOAT *J) {
  for (int i = 1; i <= 4; i++) {
    move_dir(i);
    if ((J[get_place()] > 0.3) || (J[get_place()] < -0.3)) {
      return i;
    }
  }

  for (int i = 1; i <= 4; i++) {
    move_dir(i);
    move(i, -1);
    if ((J[get_place()] > 0.3) || (J[get_place()] < -0.3)) {
      return -i;
    }
    move(i, 1);
  }
  return 0;
}*/

template class link1<su2>;
template class link1<abelian>;

template ostream &operator<<(ostream &os, const link1<su2> &link);
template ostream &operator<<(ostream &os, const link1<abelian> &link);