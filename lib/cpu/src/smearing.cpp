#include "../include/smearing.h"
#include "../include/link.h"

// #define data_size                                                              \
//   4 * link.lattice_size[0] * link.lattice_size[1] * link.lattice_size[2] *     \
//       link.lattice_size[3]
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
  link.move_dir(eta);
  A = *link.get_matrix(vec);
  link.move(eta, 1);
  link.move_dir(dir);
  A = A * link.get_matrix(vec);
  link.move(dir, 1);
  link.move(eta, -1);
  link.move_dir(eta);
  A = A ^ link.get_matrix(vec);
  link.move(dir, -1);
  link.move(eta, -1);
  B = link.get_matrix(vec)->conj();
  link.move_dir(dir);
  B = B * link.get_matrix(vec);
  link.move(dir, 1);
  link.move_dir(eta);
  B = B * link.get_matrix(vec);
  link.move(eta, 1);
  link.move(dir, -1);
  link.move_dir(dir);
  return (A + B);
}

template <class T>
T staples_second(const std::vector<std::vector<T>> &smearing_first, link1 &link,
                 std::unordered_map<int, int> &indexes, int rho, int mu,
                 int nu) {
  T A;
  T B;
  int a = 0;
  A = smearing_first[indexes[rho * 100 + nu * 10 + mu]][PLACE1_LINK_NODIR];
  link.move(rho, 1);
  A = A * smearing_first[indexes[mu * 100 + rho * 10 + nu]][PLACE1_LINK_NODIR];
  link.move(link.direction, 1);
  link.move(rho, -1);
  A = A * smearing_first[indexes[rho * 100 + nu * 10 + mu]][PLACE1_LINK_NODIR]
              .conj();
  link.move(link.direction, -1);
  link.move(rho, -1);
  B = smearing_first[indexes[rho * 100 + nu * 10 + mu]][PLACE1_LINK_NODIR]
          .conj();
  B = B * smearing_first[indexes[mu * 100 + rho * 10 + nu]][PLACE1_LINK_NODIR];
  link.move(link.direction, 1);
  B = B * smearing_first[indexes[rho * 100 + nu * 10 + mu]][PLACE1_LINK_NODIR];
  link.move(rho, 1);
  link.move(link.direction, -1);
  return (A + B);
}

template <class T>
T staples_second_refresh(const std::vector<T> &vec, link1 &link, int eta,
                         int nu, FLOAT alpha3) {
  T A;
  T B;
  int dir = link.direction;
  link.move_dir(eta);
  A = smearing_first_refresh(vec, link, nu, dir, alpha3);
  link.move(eta, 1);
  link.move_dir(dir);
  A = A * smearing_first_refresh(vec, link, nu, eta, alpha3);
  link.move(dir, 1);
  link.move(eta, -1);
  link.move_dir(eta);
  A = A * smearing_first_refresh(vec, link, nu, dir, alpha3).conj();
  link.move(dir, -1);
  link.move(eta, -1);
  B = smearing_first_refresh(vec, link, nu, dir, alpha3).conj();
  link.move_dir(dir);
  B = B * smearing_first_refresh(vec, link, nu, eta, alpha3);
  link.move(dir, 1);
  link.move_dir(eta);
  B = B * smearing_first_refresh(vec, link, nu, dir, alpha3);
  link.move(eta, 1);
  link.move(dir, -1);
  link.move_dir(dir);
  return (A + B);
}

template <class T>
T staples_third(const std::vector<std::vector<T>> &smearing_second, link1 &link,
                std::unordered_map<int, int> indexes, int nu, int mu) {
  T A;
  T B;
  A = smearing_second[indexes[nu * 10 + mu]][PLACE1_LINK_NODIR];
  link.move(nu, 1);
  A = A * smearing_second[indexes[mu * 10 + nu]][PLACE1_LINK_NODIR];
  link.move(link.direction, 1);
  link.move(nu, -1);
  A = A * smearing_second[indexes[nu * 10 + mu]][PLACE1_LINK_NODIR].conj();
  link.move(link.direction, -1);
  link.move(nu, -1);
  B = smearing_second[indexes[nu * 10 + mu]][PLACE1_LINK_NODIR].conj();
  B = B * smearing_second[indexes[mu * 10 + nu]][PLACE1_LINK_NODIR];
  link.move(link.direction, 1);
  B = B * smearing_second[indexes[nu * 10 + mu]][PLACE1_LINK_NODIR];
  link.move(nu, 1);
  link.move(link.direction, -1);
  return (A + B);
}

template <class T>
T staples_third_refresh(const std::vector<T> &vec, link1 &link, int eta,
                        FLOAT alpha2, FLOAT alpha3) {
  T A;
  T B;
  int dir = link.direction;
  link.move_dir(eta);
  A = smearing_second_refresh(vec, link, dir, alpha2, alpha3);
  link.move(eta, 1);
  link.move_dir(dir);
  A = A * smearing_second_refresh(vec, link, eta, alpha2, alpha3);
  link.move(dir, 1);
  link.move(eta, -1);
  link.move_dir(eta);
  A = A * smearing_second_refresh(vec, link, dir, alpha2, alpha3).conj();
  link.move(dir, -1);
  link.move(eta, -1);
  B = smearing_second_refresh(vec, link, dir, alpha2, alpha3).conj();
  link.move_dir(dir);
  B = B * smearing_second_refresh(vec, link, eta, alpha2, alpha3);
  link.move(dir, 1);
  link.move_dir(eta);
  B = B * smearing_second_refresh(vec, link, dir, alpha2, alpha3);
  link.move(eta, 1);
  link.move(dir, -1);
  link.move_dir(dir);
  return (A + B);
}

template <class T>
std::vector<T> smearing_first(const std::vector<T> &array, FLOAT alpha3, int mu,
                              int nu, int rho) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<T> vec(data_size / 4);
  link.move_dir(mu);
  SPACE_ITER_START
  vec[PLACE1_LINK_NODIR] = (1 - alpha3) * *link.get_matrix(array);
  for (int i = 0; i < 4; i++) {
    if (i != mu && i != nu && i != rho) {
      vec[PLACE1_LINK_NODIR] =
          vec[PLACE1_LINK_NODIR] + alpha3 / 2. * staples_first(array, link, i);
    }
  }
  vec[PLACE1_LINK_NODIR] = vec[PLACE1_LINK_NODIR].proj();
  SPACE_ITER_END
  return vec;
}

void make_map_first(std::unordered_map<int, int> &indexes) {
  int key;
  int count = 0;
  for (int i = 0; i <= 2; i++) {
    for (int j = 0; j <= 2; j++) {
      int k = 3;
      if (i != j) {
        key = i * 100 + j * 10 + k;
        indexes[key] = count;
        key = i * 100 + k * 10 + j;
        indexes[key] = count;
        count++;
      }
    }
  }
  int i = 3;
  for (int j = 0; j <= 1; j++) {
    for (int k = j + 1; k <= 2; k++) {
      key = i * 100 + j * 10 + k;
      indexes[key] = count;
      key = i * 100 + k * 10 + j;
      indexes[key] = count;
      count++;
    }
  }
}

template <class T>
std::vector<std::vector<T>> smearing_first_full(const std::vector<T> &array,
                                                FLOAT alpha3) {
  std::unordered_map<int, int> indexes;
  make_map_first(indexes);
  std::vector<std::vector<T>> smearing(9, std::vector<T>(data_size / 4));
  int key;
  for (int i = 0; i <= 2; i++) {
    for (int j = 0; j <= 2; j++) {
      int k = 3;
      if (i != j) {
        key = i * 100 + j * 10 + k;
        smearing[indexes[key]] = smearing_first(array, alpha3, i, j, k);
      }
    }
  }
  int i = 3;
  for (int j = 0; j <= 1; j++) {
    for (int k = j + 1; k <= 2; k++) {
      key = i * 100 + j * 10 + k;
      smearing[indexes[key]] = smearing_first(array, alpha3, i, j, k);
    }
  }
  return smearing;
}

template <class T>
std::vector<T> smearing_second(const std::vector<T> &array,
                               std::vector<std::vector<T>> &smearing_first,
                               FLOAT alpha2, int mu, int nu) {
  link1 link(x_size, y_size, z_size, t_size);
  std::unordered_map<int, int> indexes;
  make_map_first(indexes);
  std::vector<T> vec(data_size / 4);
  link.move_dir(mu);
  SPACE_ITER_START
  vec[PLACE1_LINK_NODIR] = (1 - alpha2) * *link.get_matrix(array);
  for (int i = 0; i < 4; i++) {
    if (i != mu && i != nu) {
      vec[PLACE1_LINK_NODIR] =
          vec[PLACE1_LINK_NODIR] +
          alpha2 / 4. *
              staples_second(smearing_first, link, indexes, i, mu, nu);
    }
  }
  vec[PLACE1_LINK_NODIR] = vec[PLACE1_LINK_NODIR].proj();
  SPACE_ITER_END
  return vec;
}

void make_map_second(std::unordered_map<int, int> &indexes) {
  int key;
  int count = 0;
  for (int i = 0; i <= 2; i++) {
    key = i * 10 + 3;
    indexes[key] = count;
    count++;
  }
  for (int i = 0; i <= 2; i++) {
    key = 30 + i;
    indexes[key] = count;
    count++;
  }
}

template <class T>
std::vector<std::vector<T>>
smearing_second_full(const std::vector<T> &array,
                     std::vector<std::vector<T>> &smearing_first,
                     FLOAT alpha2) {
  std::unordered_map<int, int> indexes;
  make_map_second(indexes);
  std::vector<std::vector<T>> smearing(6, std::vector<T>(data_size / 4));
  int key;
  for (int i = 0; i <= 2; i++) {
    key = i * 10 + 3;
    smearing[indexes[key]] =
        smearing_second(array, smearing_first, alpha2, i, 3);
  }
  for (int i = 0; i <= 2; i++) {
    key = 30 + i;
    smearing[indexes[key]] =
        smearing_second(array, smearing_first, alpha2, 3, i);
  }
  return smearing;
}

template <class T>
std::vector<T> smearing_HYP(const std::vector<T> &array,
                            std::vector<std::vector<T>> &smearing_second,
                            FLOAT alpha1) {
  link1 link(x_size, y_size, z_size, t_size);
  std::unordered_map<int, int> indexes;
  make_map_second(indexes);
  std::vector<T> vec(data_size);
  link.move_dir(3);
  SPACE_ITER_START
  vec[PLACE4_LINK_DIR] = (1 - alpha1) * *link.get_matrix(array);
  for (int i = 0; i < 3; i++) {
    vec[PLACE4_LINK_DIR] =
        vec[PLACE4_LINK_DIR] +
        alpha1 / 6. * staples_third(smearing_second, link, indexes, i, 3);
  }
  vec[PLACE4_LINK_DIR] = vec[PLACE4_LINK_DIR].proj();
  SPACE_ITER_END
  for (int d = 0; d < 3; d++) {
    link.move_dir(d);
    SPACE_ITER_START
    vec[PLACE4_LINK_DIR] = array[PLACE4_LINK_DIR];
  }
  SPACE_ITER_END
  return vec;
}

template <class T>
std::vector<T> smearing_APE(const std::vector<T> &array, FLOAT alpha_APE) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<T> vec(data_size);
  SPACE_ITER_START
  link.go(x, y, z, t);
  for (int i = 0; i < 3; i++) {
    link.move_dir(i);
    vec[PLACE4_LINK_DIR] = (1 - alpha_APE) * array[PLACE4_LINK_DIR];
    for (int d = 0; d < 3; d++) {
      if (d != i) {
        vec[PLACE4_LINK_DIR] = vec[PLACE4_LINK_DIR] +
                               (alpha_APE / 6.) * staples_first(array, link, d);
      }
    }
    vec[PLACE4_LINK_DIR] = vec[PLACE4_LINK_DIR].proj();
  }
  SPACE_ITER_END
  link.move_dir(3);
  SPACE_ITER_START
  link.go(x, y, z, t);
  vec[PLACE4_LINK_DIR] = array[PLACE4_LINK_DIR];
  SPACE_ITER_END
  return vec;
}

template <class T>
std::vector<T> smearing1_APE(const std::vector<T> &array, FLOAT alpha_APE) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<T> vec(data_size);
  SPACE_ITER_START
  link.go(x, y, z, t);
  for (int i = 0; i < 3; i++) {
    link.move_dir(i);
    vec[PLACE4_LINK_DIR] = array[PLACE4_LINK_DIR];
    for (int d = 0; d < 3; d++) {
      if (d != i) {
        vec[PLACE4_LINK_DIR] =
            vec[PLACE4_LINK_DIR] + alpha_APE * staples_first(array, link, d);
      }
    }
    vec[PLACE4_LINK_DIR] = vec[PLACE4_LINK_DIR].proj();
  }
  SPACE_ITER_END
  link.move_dir(3);
  SPACE_ITER_START
  link.go(x, y, z, t);
  vec[PLACE4_LINK_DIR] = array[PLACE4_LINK_DIR];
  SPACE_ITER_END
  return vec;
}

template <class T>
std::map<std::tuple<int, int>, std::vector<T>>
smearing_APE_2d(const std::vector<T> &array, FLOAT alpha_APE) {
  std::map<std::tuple<int, int>, std::vector<T>> smeared;

  link1 link(x_size, y_size, z_size, t_size);

  // not smeared direction
  for (int i = 0; i < 3; i++) {
    // 2 other smeared directions
    for (int j = 0; j < 3; j++) {

      if (i != j) {
        std::vector<T> vec(x_size * y_size * z_size * t_size);

        link.move_dir(j);

        // for (int t = 0; t < t_size; t++) {
        // link.move_dir(j);

        // SPACE_ITER_START_3D

        SPACE_ITER_START

        vec[link.place / 4] = (1 - alpha_APE) * array[link.place + j];

        for (int d = 0; d < 3; d++) {
          if (d != i && d != j) {
            vec[link.place / 4] =
                vec[link.place / 4] +
                (alpha_APE / 2.) * staples_first(array, link, d);
          }
        }
        vec[link.place / 4] = vec[link.place / 4].proj();

        SPACE_ITER_END

        // SPACE_ITER_END_3D
        // }
        smeared[std::tuple<int, int>{i, j}] = std::move(vec);
      }
    }
  }
  return smeared;
}

template <class T>
void smearing_APE_2d_continue(
    std::map<std::tuple<int, int>, std::vector<T>> &smeared, FLOAT alpha_APE) {
  for (int i = 0; i < 3; i++) {
    smearing_APE_2d_continue_plane(smeared, i, alpha_APE);
  }
}

template <class T>
void smearing_APE_2d_continue_plane(
    std::map<std::tuple<int, int>, std::vector<T>> &smeared, int mu,
    FLOAT alpha_APE) {
  std::vector<int> directions;

  link1 link(x_size, y_size, z_size, t_size);

  std::vector<std::vector<T>> vec(
      2, std::vector<T>(x_size * y_size * z_size * t_size));

  int nu1, nu2;

  for (int i = 0; i < 2; i++) {
    if (i != mu) {

      for (int k = i + 1; k < 3; k++) {
        if (k != mu) {
          nu1 = i;
          nu2 = k;
        }
      }
    }
  }

  SPACE_ITER_START

  vec[0][link.place / 4] =
      (1 - alpha_APE) * smeared[std::tuple<int, int>{mu, nu1}][link.place / 4];

  vec[0][link.place / 4] =
      vec[0][link.place / 4] +
      (alpha_APE / 2.) *
          staples_2d_continue(smeared[std::tuple<int, int>{mu, nu1}],
                              smeared[std::tuple<int, int>{mu, nu2}], link, nu1,
                              nu2);

  vec[0][link.place / 4] = vec[0][link.place / 4].proj();

  SPACE_ITER_END

  SPACE_ITER_START

  vec[1][link.place / 4] =
      (1 - alpha_APE) * smeared[std::tuple<int, int>{mu, nu2}][link.place / 4];

  vec[1][link.place / 4] =
      vec[1][link.place / 4] +
      (alpha_APE / 2.) *
          staples_2d_continue(smeared[std::tuple<int, int>{mu, nu2}],
                              smeared[std::tuple<int, int>{mu, nu1}], link, nu2,
                              nu1);

  vec[1][link.place / 4] = vec[1][link.place / 4].proj();

  SPACE_ITER_END

  smeared[std::tuple<int, int>{mu, nu1}] = std::move(vec[0]);
  smeared[std::tuple<int, int>{mu, nu2}] = std::move(vec[1]);
}

template <class T>
T staples_2d_continue(std::vector<T> &array1, std::vector<T> &array2,
                      link1 &link, int i, int j) {
  T A;
  T B;
  A = array2[link.place / 4];
  link.move(j, 1);
  A = A * array1[link.place / 4];
  link.move(i, 1);
  link.move(j, -1);
  A = A ^ &array2[link.place / 4];
  link.move(i, -1);
  link.move(j, -1);
  B = array2[link.place / 4].conj();
  B = B * array1[link.place / 4];
  link.move(i, 1);
  B = B * array2[link.place / 4];
  link.move(j, 1);
  link.move(i, -1);
  return (A + B);
}

template <class T>
T smearing_first_refresh(const std::vector<T> &vec, link1 &link, int nu,
                         int rho, FLOAT alpha3) {
  T A;
  A = (1 - alpha3) * *link.get_matrix(vec);
  for (int d = 0; d < 4; d++) {
    if (d != nu && d != rho && d != link.direction) {
      A = A + alpha3 / 2. * staples_first(vec, link, d);
    }
  }
  A.proj();
  return A;
}

template <class T>
T smearing_second_refresh(const std::vector<T> &vec, link1 &link, int nu,
                          FLOAT alpha2, FLOAT alpha3) {
  T A;
  A = (1 - alpha2) * *link.get_matrix(vec);
  for (int d = 0; d < 4; d++) {
    if (d != nu && d != link.direction) {
      A = A + alpha2 / 4. * staples_second_refresh(vec, link, d, nu, alpha3);
    }
  }
  A.proj();
  return A;
}

template <class T>
std::vector<T> smearing_HYP_refresh(data<T> &conf, FLOAT alpha1, FLOAT alpha2,
                                    FLOAT alpha3) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<T> vec(data_size);
  vec = conf.array;
  T A;
  link.move_dir(3);
  SPACE_ITER_START
  A = (1 - alpha1) * *link.get_matrix(vec);
  for (int d = 0; d < 3; d++) {
    A = A + alpha1 / 6. * staples_third_refresh(vec, link, d, alpha2, alpha3);
  }
  vec[PLACE4_LINK_DIR] = A.proj();
  SPACE_ITER_END
  return vec;
}

template <class T>
std::vector<T> smearing_APE_refresh(data<T> &conf, FLOAT alpha_APE) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<T> vec(data_size);
  vec = conf.array;
  T A;
  SPACE_ITER_START
  for (int i = 1; i < 4; i++) {
    link.move_dir(i);
    A = (1 - alpha_APE) * conf.array[PLACE4_LINK_DIR];
    for (int d = 0; d < 3; d++) {
      if (d != i)
        A = A + (alpha_APE / 4.) * staples_first(vec, link, d);
    }
    vec[PLACE4_LINK_DIR] = A.proj();
  }
  SPACE_ITER_END
  return vec;
}

std::vector<su2> smearing_stout(data<su2> &conf, FLOAT rho) {
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

su2 stout_factor(data<su2> &conf, link1 &link, FLOAT rho) {
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

su2 stout_omega(data<su2> &conf, link1 &link, FLOAT rho) {
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

// specialications

// su2
template su2 staples_first(const std::vector<su2> &vec, link1 &link, int eta);
template su2 staples_second(const std::vector<std::vector<su2>> &smearing_first,
                            link1 &link, std::unordered_map<int, int> &indexes,
                            int rho, int mu, int nu);
template su2
staples_second_refresh(const std::vector<su2> &vec, link1 &link, int eta,
                       int nu,
                       FLOAT alpha3); // staples for refreshing
                                      // algorythm(refresh link every step)
template su2 staples_third(const std::vector<std::vector<su2>> &smearing_second,
                           link1 &link, std::unordered_map<int, int> indexes,
                           int nu, int mu);
template su2 staples_third_refresh(const std::vector<su2> &vec, link1 &link,
                                   int eta, FLOAT alpha2, FLOAT alpha3);
template std::vector<su2> smearing_first(const std::vector<su2> &array,
                                         FLOAT alpha3, int mu, int nu, int rho);
template std::vector<std::vector<su2>>
smearing_first_full(const std::vector<su2> &array, FLOAT alpha3);
template std::vector<su2>
smearing_second(const std::vector<su2> &array,
                std::vector<std::vector<su2>> &smearing_first, FLOAT alpha2,
                int mu, int nu);
template std::vector<std::vector<su2>>
smearing_second_full(const std::vector<su2> &array,
                     std::vector<std::vector<su2>> &smearing_first,
                     FLOAT alpha2);
template std::vector<su2>
smearing_HYP(const std::vector<su2> &array,
             std::vector<std::vector<su2>> &smearing_second, FLOAT alpha1);
template std::vector<su2> smearing_APE(const std::vector<su2> &array,
                                       FLOAT alpha_APE);
template std::vector<su2> smearing1_APE(const std::vector<su2> &array,
                                        FLOAT alpha_APE);
template std::map<std::tuple<int, int>, std::vector<su2>>
smearing_APE_2d(const std::vector<su2> &array, FLOAT alpha_APE);
template void smearing_APE_2d_continue(
    std::map<std::tuple<int, int>, std::vector<su2>> &smeared, FLOAT alpha_APE);
template void smearing_APE_2d_continue_plane(
    std::map<std::tuple<int, int>, std::vector<su2>> &smeared, int mu,
    FLOAT alpha_APE);
template su2 staples_2d_continue(std::vector<su2> &array1,
                                 std::vector<su2> &array2, link1 &link, int i,
                                 int j);
template su2 smearing_first_refresh(const std::vector<su2> &vec, link1 &link,
                                    int nu, int rho,
                                    FLOAT alpha3); // refresh link every step
template su2 smearing_second_refresh(const std::vector<su2> &vec, link1 &link,
                                     int nu, FLOAT alpha2,
                                     FLOAT alpha3); // refresh link every step
template std::vector<su2>
smearing_HYP_refresh(data<su2> &conf, FLOAT alpha1, FLOAT alpha2,
                     FLOAT alpha3); // refresh link every step
template std::vector<su2>
smearing_APE_refresh(data<su2> &conf,
                     FLOAT alpha_APE); // refresh link every step

// abelian
template abelian staples_first(const std::vector<abelian> &vec, link1 &link,
                               int eta);
template abelian
staples_second(const std::vector<std::vector<abelian>> &smearing_first,
               link1 &link, std::unordered_map<int, int> &indexes, int rho,
               int mu, int nu);
template abelian
staples_second_refresh(const std::vector<abelian> &vec, link1 &link, int eta,
                       int nu,
                       FLOAT alpha3); // staples for refreshing
                                      // algorythm(refresh link every step)
template abelian
staples_third(const std::vector<std::vector<abelian>> &smearing_second,
              link1 &link, std::unordered_map<int, int> indexes, int nu,
              int mu);
template abelian staples_third_refresh(const std::vector<abelian> &vec,
                                       link1 &link, int eta, FLOAT alpha2,
                                       FLOAT alpha3);
template std::vector<abelian> smearing_first(const std::vector<abelian> &array,
                                             FLOAT alpha3, int mu, int nu,
                                             int rho);
template std::vector<std::vector<abelian>>
smearing_first_full(const std::vector<abelian> &array, FLOAT alpha3);
template std::vector<abelian>
smearing_second(const std::vector<abelian> &array,
                std::vector<std::vector<abelian>> &smearing_first, FLOAT alpha2,
                int mu, int nu);
template std::vector<std::vector<abelian>>
smearing_second_full(const std::vector<abelian> &array,
                     std::vector<std::vector<abelian>> &smearing_first,
                     FLOAT alpha2);
template std::vector<abelian>
smearing_HYP(const std::vector<abelian> &array,
             std::vector<std::vector<abelian>> &smearing_second, FLOAT alpha1);
template std::vector<abelian> smearing_APE(const std::vector<abelian> &array,
                                           FLOAT alpha_APE);
template std::vector<abelian> smearing1_APE(const std::vector<abelian> &array,
                                            FLOAT alpha_APE);
template std::map<std::tuple<int, int>, std::vector<abelian>>
smearing_APE_2d(const std::vector<abelian> &array, FLOAT alpha_APE);
template void smearing_APE_2d_continue(
    std::map<std::tuple<int, int>, std::vector<abelian>> &smeared,
    FLOAT alpha_APE);
template void smearing_APE_2d_continue_plane(
    std::map<std::tuple<int, int>, std::vector<abelian>> &smeared, int mu,
    FLOAT alpha_APE);
template abelian staples_2d_continue(std::vector<abelian> &array1,
                                     std::vector<abelian> &array2, link1 &link,
                                     int i, int j);
template abelian
smearing_first_refresh(const std::vector<abelian> &vec, link1 &link, int nu,
                       int rho,
                       FLOAT alpha3); // refresh link every step
template abelian
smearing_second_refresh(const std::vector<abelian> &vec, link1 &link, int nu,
                        FLOAT alpha2,
                        FLOAT alpha3); // refresh link every step
template std::vector<abelian>
smearing_HYP_refresh(data<abelian> &conf, FLOAT alpha1, FLOAT alpha2,
                     FLOAT alpha3); // refresh link every step
template std::vector<abelian>
smearing_APE_refresh(data<abelian> &conf,
                     FLOAT alpha_APE); // refresh link every step