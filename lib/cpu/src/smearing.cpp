#include "../include/smearing.h"
#include "../include/link.h"
#include "../include/matrix.h"

#include <algorithm>
#include <execution>
#include <omp.h>
#include <utility>
#include <vector>

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "sys/sysinfo.h"
#include "sys/types.h"

// int getValue() { // Note: this value is in KB!
//   FILE *file = fopen("/proc/self/status", "r");
//   int result = -1;
//   char line[128];

//   while (fgets(line, 128, file) != NULL) {
//     if (strncmp(line, "VmRSS:", 6) == 0) {
//       result = ParseLine(line);
//       break;
//     }
//   }
//   fclose(file);
//   return result;
// }

struct sysinfo memInfo;

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
                         int nu, double alpha3) {
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
                        double alpha2, double alpha3) {
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
std::vector<T> smearing_first(const std::vector<T> &array, double alpha3,
                              int mu, int nu, int rho) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<T> vec(data_size / 4);
  link.move_dir(mu);
  SPACE_ITER_START
  vec[PLACE1_LINK_NODIR] = (1 - alpha3) * array[link.place + mu];
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
                                                double alpha3) {
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
                               double alpha2, int mu, int nu) {
  link1 link(x_size, y_size, z_size, t_size);
  std::unordered_map<int, int> indexes;
  make_map_first(indexes);
  std::vector<T> vec(data_size / 4);
  link.move_dir(mu);
  SPACE_ITER_START
  vec[link.place / 4] = (1 - alpha2) * array[link.place + mu];
  for (int i = 0; i < 4; i++) {
    if (i != mu && i != nu) {
      vec[link.place / 4] =
          vec[link.place / 4] +
          alpha2 / 4. *
              staples_second(smearing_first, link, indexes, i, mu, nu);
    }
  }
  vec[link.place / 4] = vec[link.place / 4].proj();
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
                     double alpha2) {
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
                            double alpha1) {
  link1 link(x_size, y_size, z_size, t_size);
  std::unordered_map<int, int> indexes;
  make_map_second(indexes);
  std::vector<T> vec(data_size);
  link.move_dir(3);
  SPACE_ITER_START
  vec[link.place + 3] = (1 - alpha1) * array[link.place + 3];
  for (int i = 0; i < 3; i++) {
    vec[link.place + 3] =
        vec[link.place + 3] +
        alpha1 / 6. * staples_third(smearing_second, link, indexes, i, 3);
  }
  vec[link.place + 3] = vec[link.place + 3].proj();
  SPACE_ITER_END
  for (int d = 0; d < 3; d++) {
    SPACE_ITER_START
    vec[link.place + d] = array[link.place + d];
    SPACE_ITER_END
  }
  return vec;
}

template <class T>
std::vector<T> smearing_APE(const std::vector<T> &array, double alpha_APE) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<T> vec(data_size);
  SPACE_ITER_START
  link.go(x, y, z, t);
  for (int i = 0; i < 3; i++) {
    link.move_dir(i);
    vec[link.place + i] = (1 - alpha_APE) * array[link.place + i];
    for (int d = 0; d < 3; d++) {
      if (d != i) {
        vec[link.place + i] = vec[link.place + i] +
                              (alpha_APE / 6.) * staples_first(array, link, d);
      }
    }
    vec[link.place + i] = vec[link.place + i].proj();
  }
  SPACE_ITER_END
  SPACE_ITER_START
  link.go(x, y, z, t);
  vec[link.place + 3] = array[link.place + 3];
  SPACE_ITER_END
  return vec;
}

template <class T>
std::vector<T> smearing1_APE(const std::vector<T> &array, double alpha_APE) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<T> vec(data_size);
  SPACE_ITER_START
  link.go(x, y, z, t);
  for (int i = 0; i < 3; i++) {
    link.move_dir(i);
    vec[link.place + i] = array[link.place + i];
    for (int d = 0; d < 3; d++) {
      if (d != i) {
        vec[link.place + i] =
            vec[link.place + i] + alpha_APE * staples_first(array, link, d);
      }
    }
    vec[link.place + i] = vec[link.place + i].proj();
  }
  SPACE_ITER_END
  SPACE_ITER_START
  link.go(x, y, z, t);
  vec[link.place + 3] = array[link.place + 3];
  SPACE_ITER_END
  return vec;
}

template <class T>
T smearing_first_refresh(const std::vector<T> &vec, link1 &link, int nu,
                         int rho, double alpha3) {
  T A;
  A = (1 - alpha3) * vec[link.place + link.direction];
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
                          double alpha2, double alpha3) {
  T A;
  A = (1 - alpha2) * vec[link.place + link.direction];
  for (int d = 0; d < 4; d++) {
    if (d != nu && d != link.direction) {
      A = A + alpha2 / 4. * staples_second_refresh(vec, link, d, nu, alpha3);
    }
  }
  A.proj();
  return A;
}

template <class T>
std::vector<T> smearing_HYP_refresh(data<T> &conf, double alpha1, double alpha2,
                                    double alpha3) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<T> vec(data_size);
  vec = conf.array;
  T A;
  link.move_dir(3);
  SPACE_ITER_START
  A = (1 - alpha1) * vec[link.place + 3];
  for (int d = 0; d < 3; d++) {
    A = A + alpha1 / 6. * staples_third_refresh(vec, link, d, alpha2, alpha3);
  }
  vec[PLACE4_LINK_DIR] = A.proj();
  SPACE_ITER_END
  return vec;
}

template <class T>
std::vector<T> smearing_APE_refresh(data<T> &conf, double alpha_APE) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<T> vec(data_size);
  vec = conf.array;
  T A;
  SPACE_ITER_START
  for (int i = 1; i < 4; i++) {
    link.move_dir(i);
    A = conf.array[PLACE4_LINK_DIR];
    for (int d = 0; d < 3; d++) {
      if (d != i)
        A = A + alpha_APE * staples_first(vec, link, d);
    }
    vec[PLACE4_LINK_DIR] = A.proj();
  }
  SPACE_ITER_END
  return vec;
}

std::vector<su2> smearing_stout(data<su2> &conf, double rho) {
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

su2 stout_factor(data<su2> &conf, link1 &link, double rho) {
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

su2 stout_omega(data<su2> &conf, link1 &link, double rho) {
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

template <class T>
void smearing_plane_minor(std::vector<T> &smeared,
                          const std::vector<T> &conf_mu,
                          const std::vector<T> &conf_nu, int size_mu1,
                          int size_mu2, int size_nu1, int size_nu2,
                          double alpha) {
  int data_size1 = x_size * y_size * z_size * t_size;

  T bracket;

#pragma omp parallel for collapse(3) private(bracket) firstprivate(alpha)
  for (int k = 0; k < data_size1; k += size_nu2) {
    for (int i = 0; i < size_nu2; i += size_mu2) {
      for (int j = 0; j < size_mu2; j++) {
        if (i < size_nu2 - size_nu1) {
          bracket = conf_nu[i + k + j] * conf_mu[i + k + j + size_nu1];
        } else
          bracket =
              conf_nu[i + k + j] * conf_mu[i + k + j - size_nu2 + size_nu1];
        if (j < size_mu2 - size_mu1) {
          smeared[i + k + j] =
              smeared[i + k + j] +
              alpha * (bracket ^ conf_nu[i + k + j + size_mu1]);
        } else
          smeared[i + k + j] =
              smeared[i + k + j] +
              alpha * (bracket ^ conf_nu[i + k + j - size_mu2 + size_mu1]);

        if (i >= size_nu1) {
          bracket =
              conf_nu[i + k + j - size_nu1] % conf_mu[i + k + j - size_nu1];
          if (j < size_mu2 - size_mu1)
            smeared[i + k + j] =
                smeared[i + k + j] +
                alpha * (bracket * conf_nu[i + k + j + size_mu1 - size_nu1]);
          else
            smeared[i + k + j] =
                smeared[i + k + j] +
                alpha * (bracket *
                         conf_nu[i + k + j - size_mu2 + size_mu1 - size_nu1]);
        } else {
          bracket = conf_nu[i + k + j + size_nu2 - size_nu1] %
                    conf_mu[i + k + j + size_nu2 - size_nu1];
          if (j < size_mu2 - size_mu1)
            smeared[i + k + j] =
                smeared[i + k + j] +
                alpha * (bracket *
                         conf_nu[i + k + j + size_nu2 - size_nu1 + size_mu1]);
          else
            smeared[i + k + j] =
                smeared[i + k + j] +
                alpha * (bracket * conf_nu[i + k + j + size_nu2 - size_nu1 -
                                           size_mu2 + size_mu1]);
        }
      }
    }
  }
}

template <class T>
void smearing_plane_minor_start(std::vector<T> &smeared,
                                const std::vector<T> &conf_mu,
                                const std::vector<T> &conf_nu, int size_mu1,
                                int size_mu2, int size_nu1, int size_nu2,
                                double alpha) {
  int data_size1 = x_size * y_size * z_size * t_size;

  T bracket;

#pragma omp parallel for collapse(3) private(bracket) firstprivate(alpha)
  for (int k = 0; k < data_size1; k += size_nu2) {
    for (int i = 0; i < size_nu2; i += size_mu2) {
      for (int j = 0; j < size_mu2; j++) {
        if (i < size_nu2 - size_nu1) {
          bracket = conf_nu[i + k + j] * conf_mu[i + k + j + size_nu1];
        } else
          bracket =
              conf_nu[i + k + j] * conf_mu[i + k + j - size_nu2 + size_nu1];
        if (j < size_mu2 - size_mu1) {
          smeared[i + k + j] =
              conf_mu[i + k + j] +
              alpha * (bracket ^ conf_nu[i + k + j + size_mu1]);
        } else
          smeared[i + k + j] =
              conf_mu[i + k + j] +
              alpha * (bracket ^ conf_nu[i + k + j - size_mu2 + size_mu1]);
        if (i >= size_nu1) {
          bracket =
              conf_nu[i + k + j - size_nu1] % conf_mu[i + k + j - size_nu1];
          if (j < size_mu2 - size_mu1)
            smeared[i + k + j] =
                smeared[i + k + j] +
                alpha * (bracket * conf_nu[i + k + j + size_mu1 - size_nu1]);
          else
            smeared[i + k + j] =
                smeared[i + k + j] +
                alpha * (bracket *
                         conf_nu[i + k + j - size_mu2 + size_mu1 - size_nu1]);
        } else {
          bracket = conf_nu[i + k + j + size_nu2 - size_nu1] %
                    conf_mu[i + k + j + size_nu2 - size_nu1];
          if (j < size_mu2 - size_mu1)
            smeared[i + k + j] =
                smeared[i + k + j] +
                alpha * (bracket *
                         conf_nu[i + k + j + size_nu2 - size_nu1 + size_mu1]);
          else
            smeared[i + k + j] =
                smeared[i + k + j] +
                alpha * (bracket * conf_nu[i + k + j + size_nu2 - size_nu1 -
                                           size_mu2 + size_mu1]);
        }
      }
    }
  }
}

template <class T>
void smearing_plane_minor_start_proj(std::vector<T> &smeared,
                                     const std::vector<T> &conf_mu,
                                     const std::vector<T> &conf_nu,
                                     int size_mu1, int size_mu2, int size_nu1,
                                     int size_nu2, double alpha) {
  int data_size1 = x_size * y_size * z_size * t_size;

  T bracket;

#pragma omp parallel for collapse(3) private(bracket) firstprivate(alpha)
  for (int k = 0; k < data_size1; k += size_nu2) {
    for (int i = 0; i < size_nu2; i += size_mu2) {
      for (int j = 0; j < size_mu2; j++) {
        if (i < size_nu2 - size_nu1) {
          bracket = conf_nu[i + k + j] * conf_mu[i + k + j + size_nu1];
        } else
          bracket =
              conf_nu[i + k + j] * conf_mu[i + k + j - size_nu2 + size_nu1];
        if (j < size_mu2 - size_mu1) {
          smeared[i + k + j] =
              conf_mu[i + k + j] +
              alpha * (bracket ^ conf_nu[i + k + j + size_mu1]);
        } else
          smeared[i + k + j] =
              conf_mu[i + k + j] +
              alpha * (bracket ^ conf_nu[i + k + j - size_mu2 + size_mu1]);
        if (i >= size_nu1) {
          bracket =
              conf_nu[i + k + j - size_nu1] % conf_mu[i + k + j - size_nu1];
          if (j < size_mu2 - size_mu1)
            smeared[i + k + j] =
                (smeared[i + k + j] +
                 alpha * (bracket * conf_nu[i + k + j + size_mu1 - size_nu1]))
                    .proj();
          else
            smeared[i + k + j] =
                (smeared[i + k + j] +
                 alpha * (bracket *
                          conf_nu[i + k + j - size_mu2 + size_mu1 - size_nu1]))
                    .proj();
        } else {
          bracket = conf_nu[i + k + j + size_nu2 - size_nu1] %
                    conf_mu[i + k + j + size_nu2 - size_nu1];
          if (j < size_mu2 - size_mu1)
            smeared[i + k + j] =
                (smeared[i + k + j] +
                 alpha * (bracket *
                          conf_nu[i + k + j + size_nu2 - size_nu1 + size_mu1]))
                    .proj();
          else
            smeared[i + k + j] =
                (smeared[i + k + j] +
                 alpha * (bracket * conf_nu[i + k + j + size_nu2 - size_nu1 -
                                            size_mu2 + size_mu1]))
                    .proj();
        }
      }
    }
  }
}

template <class T>
void smearing_plane_minor_end(std::vector<T> &smeared,
                              const std::vector<T> &conf_mu,
                              const std::vector<T> &conf_nu, int size_mu1,
                              int size_mu2, int size_nu1, int size_nu2,
                              double alpha) {
  int data_size1 = x_size * y_size * z_size * t_size;

  T bracket;

#pragma omp parallel for collapse(3) private(bracket) firstprivate(alpha)
  for (int k = 0; k < data_size1; k += size_nu2) {
    for (int i = 0; i < size_nu2; i += size_mu2) {
      for (int j = 0; j < size_mu2; j++) {
        if (i < size_nu2 - size_nu1) {
          bracket = conf_nu[i + k + j] * conf_mu[i + k + j + size_nu1];
        } else
          bracket =
              conf_nu[i + k + j] * conf_mu[i + k + j - size_nu2 + size_nu1];
        if (j < size_mu2 - size_mu1) {
          smeared[i + k + j] =
              smeared[i + k + j] +
              alpha * (bracket ^ conf_nu[i + k + j + size_mu1]);
        } else
          smeared[i + k + j] =
              smeared[i + k + j] +
              alpha * (bracket ^ conf_nu[i + k + j - size_mu2 + size_mu1]);

        if (i >= size_nu1) {
          bracket =
              conf_nu[i + k + j - size_nu1] % conf_mu[i + k + j - size_nu1];
          if (j < size_mu2 - size_mu1)
            smeared[i + k + j] =
                (smeared[i + k + j] +
                 alpha * (bracket * conf_nu[i + k + j + size_mu1 - size_nu1]))
                    .proj();
          else
            smeared[i + k + j] =
                (smeared[i + k + j] +
                 alpha * (bracket *
                          conf_nu[i + k + j - size_mu2 + size_mu1 - size_nu1]))
                    .proj();
        } else {
          bracket = conf_nu[i + k + j + size_nu2 - size_nu1] %
                    conf_mu[i + k + j + size_nu2 - size_nu1];
          if (j < size_mu2 - size_mu1)
            smeared[i + k + j] =
                (smeared[i + k + j] +
                 alpha * (bracket *
                          conf_nu[i + k + j + size_nu2 - size_nu1 + size_mu1]))
                    .proj();
          else
            smeared[i + k + j] =
                (smeared[i + k + j] +
                 alpha * (bracket * conf_nu[i + k + j + size_nu2 - size_nu1 -
                                            size_mu2 + size_mu1]))
                    .proj();
        }
      }
    }
  }
}

template <class T>
void smearing_plane_major(std::vector<T> &smeared,
                          const std::vector<T> &conf_mu,
                          const std::vector<T> &conf_nu, int size_mu1,
                          int size_mu2, int size_nu1, int size_nu2,
                          double alpha) {
  int data_size1 = x_size * y_size * z_size * t_size;

  T bracket;

#pragma omp parallel for collapse(3) private(bracket) firstprivate(alpha)
  for (int k = 0; k < data_size1; k += size_mu2) {
    for (int i = 0; i < size_mu2; i += size_nu2) {
      for (int j = 0; j < size_nu2; j++) {
        if (j < size_nu2 - size_nu1)
          bracket = conf_nu[i + k + j] * conf_mu[i + k + j + size_nu1];
        else
          bracket =
              conf_nu[i + k + j] * conf_mu[i + k + j - size_nu2 + size_nu1];
        if (i < size_mu2 - size_mu1)
          smeared[i + k + j] =
              smeared[i + k + j] +
              alpha * (bracket ^ conf_nu[i + k + j + size_mu1]);
        else
          smeared[i + k + j] =
              smeared[i + k + j] +
              alpha * (bracket ^ conf_nu[i + k + j - size_mu2 + size_mu1]);

        if (j >= size_nu1) {
          bracket =
              conf_nu[i + k + j - size_nu1] % conf_mu[i + k + j - size_nu1];
          if (i < size_mu2 - size_mu1) {
            smeared[i + k + j] =
                smeared[i + k + j] +
                alpha * (bracket * conf_nu[i + k + j + size_mu1 - size_nu1]);

          } else
            smeared[i + k + j] =
                smeared[i + k + j] +
                alpha * (bracket *
                         conf_nu[i + k + j - size_mu2 + size_mu1 - size_nu1]);
        } else {
          bracket = conf_nu[i + k + j + size_nu2 - size_nu1] %
                    conf_mu[i + k + j + size_nu2 - size_nu1];
          if (i < size_mu2 - size_mu1)
            smeared[i + k + j] =
                smeared[i + k + j] +
                alpha * (bracket *
                         conf_nu[i + k + j + size_nu2 - size_nu1 + size_mu1]);
          else
            smeared[i + k + j] =
                smeared[i + k + j] +
                alpha * (bracket * conf_nu[i + k + j + size_nu2 - size_nu1 -
                                           size_mu2 + size_mu1]);
        }
      }
    }
  }
}

template <class T>
void smearing_plane_major_start(std::vector<T> &smeared,
                                const std::vector<T> &conf_mu,
                                const std::vector<T> &conf_nu, int size_mu1,
                                int size_mu2, int size_nu1, int size_nu2,
                                double alpha) {
  int data_size1 = x_size * y_size * z_size * t_size;

  T bracket;

#pragma omp parallel for collapse(3) private(bracket) firstprivate(alpha)
  for (int k = 0; k < data_size1; k += size_mu2) {
    for (int i = 0; i < size_mu2; i += size_nu2) {
      for (int j = 0; j < size_nu2; j++) {
        if (j < size_nu2 - size_nu1)
          bracket = conf_nu[i + k + j] * conf_mu[i + k + j + size_nu1];
        else
          bracket =
              conf_nu[i + k + j] * conf_mu[i + k + j - size_nu2 + size_nu1];
        if (i < size_mu2 - size_mu1)
          smeared[i + k + j] =
              conf_mu[i + k + j] +
              alpha * (bracket ^ conf_nu[i + k + j + size_mu1]);
        else
          smeared[i + k + j] =
              conf_mu[i + k + j] +
              alpha * (bracket ^ conf_nu[i + k + j - size_mu2 + size_mu1]);

        if (j >= size_nu1) {
          bracket =
              conf_nu[i + k + j - size_nu1] % conf_mu[i + k + j - size_nu1];
          if (i < size_mu2 - size_mu1) {
            smeared[i + k + j] =
                smeared[i + k + j] +
                alpha * (bracket * conf_nu[i + k + j + size_mu1 - size_nu1]);

          } else
            smeared[i + k + j] =
                smeared[i + k + j] +
                alpha * (bracket *
                         conf_nu[i + k + j - size_mu2 + size_mu1 - size_nu1]);
        } else {
          bracket = conf_nu[i + k + j + size_nu2 - size_nu1] %
                    conf_mu[i + k + j + size_nu2 - size_nu1];
          if (i < size_mu2 - size_mu1)
            smeared[i + k + j] =
                smeared[i + k + j] +
                alpha * (bracket *
                         conf_nu[i + k + j + size_nu2 - size_nu1 + size_mu1]);
          else
            smeared[i + k + j] =
                smeared[i + k + j] +
                alpha * (bracket * conf_nu[i + k + j + size_nu2 - size_nu1 -
                                           size_mu2 + size_mu1]);
        }
      }
    }
  }
}

template <class T>
void smearing_plane_major_start_proj(std::vector<T> &smeared,
                                     const std::vector<T> &conf_mu,
                                     const std::vector<T> &conf_nu,
                                     int size_mu1, int size_mu2, int size_nu1,
                                     int size_nu2, double alpha) {
  int data_size1 = x_size * y_size * z_size * t_size;

  T bracket;

#pragma omp parallel for collapse(3) private(bracket) firstprivate(alpha)
  for (int k = 0; k < data_size1; k += size_mu2) {
    for (int i = 0; i < size_mu2; i += size_nu2) {
      for (int j = 0; j < size_nu2; j++) {
        if (j < size_nu2 - size_nu1)
          bracket = conf_nu[i + k + j] * conf_mu[i + k + j + size_nu1];
        else
          bracket =
              conf_nu[i + k + j] * conf_mu[i + k + j - size_nu2 + size_nu1];
        if (i < size_mu2 - size_mu1)
          smeared[i + k + j] =
              conf_mu[i + k + j] +
              alpha * (bracket ^ conf_nu[i + k + j + size_mu1]);
        else
          smeared[i + k + j] =
              conf_mu[i + k + j] +
              alpha * (bracket ^ conf_nu[i + k + j - size_mu2 + size_mu1]);

        if (j >= size_nu1) {
          bracket =
              conf_nu[i + k + j - size_nu1] % conf_mu[i + k + j - size_nu1];
          if (i < size_mu2 - size_mu1) {
            smeared[i + k + j] =
                (smeared[i + k + j] +
                 alpha * (bracket * conf_nu[i + k + j + size_mu1 - size_nu1]))
                    .proj();

          } else
            smeared[i + k + j] =
                (smeared[i + k + j] +
                 alpha * (bracket *
                          conf_nu[i + k + j - size_mu2 + size_mu1 - size_nu1]))
                    .proj();
        } else {
          bracket = conf_nu[i + k + j + size_nu2 - size_nu1] %
                    conf_mu[i + k + j + size_nu2 - size_nu1];
          if (i < size_mu2 - size_mu1)
            smeared[i + k + j] =
                (smeared[i + k + j] +
                 alpha * (bracket *
                          conf_nu[i + k + j + size_nu2 - size_nu1 + size_mu1]))
                    .proj();
          else
            smeared[i + k + j] =
                (smeared[i + k + j] +
                 alpha * (bracket * conf_nu[i + k + j + size_nu2 - size_nu1 -
                                            size_mu2 + size_mu1]))
                    .proj();
        }
      }
    }
  }
}

template <class T>
void smearing_plane_major_end(std::vector<T> &smeared,
                              const std::vector<T> &conf_mu,
                              const std::vector<T> &conf_nu, int size_mu1,
                              int size_mu2, int size_nu1, int size_nu2,
                              double alpha) {
  int data_size1 = x_size * y_size * z_size * t_size;

  T bracket;

#pragma omp parallel for collapse(3) private(bracket) firstprivate(alpha)
  for (int k = 0; k < data_size1; k += size_mu2) {
    for (int i = 0; i < size_mu2; i += size_nu2) {
      for (int j = 0; j < size_nu2; j++) {
        if (j < size_nu2 - size_nu1)
          bracket = conf_nu[i + k + j] * conf_mu[i + k + j + size_nu1];
        else
          bracket =
              conf_nu[i + k + j] * conf_mu[i + k + j - size_nu2 + size_nu1];
        if (i < size_mu2 - size_mu1)
          smeared[i + k + j] =
              smeared[i + k + j] +
              alpha * (bracket ^ conf_nu[i + k + j + size_mu1]);
        else
          smeared[i + k + j] =
              smeared[i + k + j] +
              alpha * (bracket ^ conf_nu[i + k + j - size_mu2 + size_mu1]);

        if (j >= size_nu1) {
          bracket =
              conf_nu[i + k + j - size_nu1] % conf_mu[i + k + j - size_nu1];
          if (i < size_mu2 - size_mu1) {
            smeared[i + k + j] =
                (smeared[i + k + j] +
                 alpha * (bracket * conf_nu[i + k + j + size_mu1 - size_nu1]))
                    .proj();

          } else
            smeared[i + k + j] =
                (smeared[i + k + j] +
                 alpha * (bracket *
                          conf_nu[i + k + j - size_mu2 + size_mu1 - size_nu1]))
                    .proj();
        } else {
          bracket = conf_nu[i + k + j + size_nu2 - size_nu1] %
                    conf_mu[i + k + j + size_nu2 - size_nu1];
          if (i < size_mu2 - size_mu1)
            smeared[i + k + j] =
                (smeared[i + k + j] +
                 alpha * (bracket *
                          conf_nu[i + k + j + size_nu2 - size_nu1 + size_mu1]))
                    .proj();
          else
            smeared[i + k + j] =
                (smeared[i + k + j] +
                 alpha * (bracket * conf_nu[i + k + j + size_nu2 - size_nu1 -
                                            size_mu2 + size_mu1]))
                    .proj();
        }
      }
    }
  }
}

template <class T>
void smearing_APE_parallel(std::vector<std::vector<T>> &conf, double alpha) {
  std::vector<std::vector<T>> smeared(3, std::vector<T>(conf[0].size()));

  smearing_plane_minor_start(smeared[0], conf[0], conf[1], 1, x_size, x_size,
                             x_size * y_size, alpha);
  smearing_plane_minor_end(smeared[0], conf[0], conf[2], 1, x_size,
                           x_size * y_size, x_size * y_size * z_size, alpha);
  smearing_plane_major_start(smeared[1], conf[1], conf[0], x_size,
                             x_size * y_size, 1, x_size, alpha);
  smearing_plane_minor_end(smeared[1], conf[1], conf[2], x_size,
                           x_size * y_size, x_size * y_size,
                           x_size * y_size * z_size, alpha);
  smearing_plane_major_start(smeared[2], conf[2], conf[0], x_size * y_size,
                             x_size * y_size * z_size, 1, x_size, alpha);
  smearing_plane_major_end(smeared[2], conf[2], conf[1], x_size * y_size,
                           x_size * y_size * z_size, x_size, x_size * y_size,
                           alpha);

  for (int i = 0; i < 3; i++) {
    conf[i] = std::move(smeared[i]);
  }
  for (int i = 0; i < smeared.size(); i++) {
    smeared[i].clear();
    smeared[i].shrink_to_fit();
  }
}

std::map<std::tuple<int, int>, int> indices_map_APE_2d() {
  std::map<std::tuple<int, int>, int> indices_map;

  int count = 0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i != j) {
        indices_map[std::tuple<int, int>(i, j)] = count;
        count++;
      }
    }
  }
  return indices_map;
}

template <class T>
std::vector<std::vector<T>>
smearing_APE_2d_initial(std::vector<std::vector<T>> &conf, double alpha) {
  int vector_size = x_size * y_size * z_size * t_size;
  std::vector<std::vector<T>> smeared(6, std::vector<T>(vector_size));
  std::map<std::tuple<int, int>, int> indices_map = indices_map_APE_2d();

  std::vector<int> steps = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  // i is quark line direction, j is smeared direction, k is parallel to them
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i != j) {
        for (int k = 0; k < 3; k++) {
          if (i != k && j != k) {
            if (j < k) {
              smearing_plane_minor_start_proj(
                  smeared[indices_map[std::tuple<int, int>(i, j)]], conf[j],
                  conf[k], steps[j], steps[j + 1], steps[k], steps[k + 1],
                  alpha);
            } else {
              smearing_plane_major_start_proj(
                  smeared[indices_map[std::tuple<int, int>(i, j)]], conf[j],
                  conf[k], steps[j], steps[j + 1], steps[k], steps[k + 1],
                  alpha);
            }
          }
        }
      }
    }
  }
  return smeared;
}

template <class T>
void smearing_APE_2d(std::vector<std::vector<T>> &conf, double alpha) {
  int vector_size = x_size * y_size * z_size * t_size;
  std::vector<std::vector<T>> smeared(6, std::vector<T>(vector_size));
  std::map<std::tuple<int, int>, int> indices_map = indices_map_APE_2d();

  std::vector<int> steps = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  // i is quark line direction, j is smeared direction, k is parallel to them
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i != j) {
        for (int k = 0; k < 3; k++) {
          if (i != k && j != k) {
            if (j < k) {
              smearing_plane_minor_start_proj(
                  smeared[indices_map[std::tuple<int, int>(i, j)]],
                  conf[indices_map[std::tuple<int, int>(i, j)]],
                  conf[indices_map[std::tuple<int, int>(i, k)]], steps[j],
                  steps[j + 1], steps[k], steps[k + 1], alpha);
            } else {
              smearing_plane_major_start_proj(
                  smeared[indices_map[std::tuple<int, int>(i, j)]],
                  conf[indices_map[std::tuple<int, int>(i, j)]],
                  conf[indices_map[std::tuple<int, int>(i, k)]], steps[j],
                  steps[j + 1], steps[k], steps[k + 1], alpha);
            }
          }
        }
      }
    }
  }
  for (int i = 0; i < 6; i++) {
    conf[i] = std::move(smeared[i]);
  }
}

template <class T>
void smearing_plane_HYP_minor(std::vector<T> &smeared, std::vector<T> &conf_mu,
                              std::vector<T> &conf_nu, int size_mu1,
                              int size_mu2, int size_nu1, int size_nu2,
                              double alpha, double divisor) {
  int data_size1 = x_size * y_size * z_size * t_size;

  T bracket;

#pragma omp parallel for collapse(3) private(bracket)
  for (int k = 0; k < data_size1; k += size_nu2) {
    for (int i = 0; i < size_nu2; i += size_mu2) {
      for (int j = 0; j < size_mu2; j++) {
        if (i < size_nu2 - size_nu1) {
          bracket = conf_nu[i + k + j] * conf_mu[i + k + j + size_nu1];
        } else
          bracket =
              conf_nu[i + k + j] * conf_mu[i + k + j - size_nu2 + size_nu1];
        if (j < size_mu2 - size_mu1) {

          smeared[i + k + j] =
              smeared[i + k + j] +
              (alpha / divisor) * (bracket ^ conf_nu[i + k + j + size_mu1]);
        } else
          smeared[i + k + j] =
              smeared[i + k + j] +
              (alpha / divisor) *
                  (bracket ^ conf_nu[i + k + j - size_mu2 + size_mu1]);

        if (i >= size_nu1) {
          bracket =
              conf_nu[i + k + j - size_nu1] % conf_mu[i + k + j - size_nu1];
          if (j < size_mu2 - size_mu1)
            smeared[i + k + j] =
                smeared[i + k + j] +
                (alpha / divisor) *
                    (bracket * conf_nu[i + k + j + size_mu1 - size_nu1]);
          else
            smeared[i + k + j] =
                smeared[i + k + j] +
                (alpha / divisor) *
                    (bracket *
                     conf_nu[i + k + j - size_mu2 + size_mu1 - size_nu1]);
        } else {
          bracket = conf_nu[i + k + j + size_nu2 - size_nu1] %
                    conf_mu[i + k + j + size_nu2 - size_nu1];
          if (j < size_mu2 - size_mu1)
            smeared[i + k + j] =
                smeared[i + k + j] +
                (alpha / divisor) *
                    (bracket *
                     conf_nu[i + k + j + size_nu2 - size_nu1 + size_mu1]);
          else
            smeared[i + k + j] =
                smeared[i + k + j] +
                (alpha / divisor) *
                    (bracket * conf_nu[i + k + j + size_nu2 - size_nu1 -
                                       size_mu2 + size_mu1]);
        }
      }
    }
  }
}

template <class T>
void smearing_plane_HYP_major(std::vector<T> &smeared, std::vector<T> &conf_mu,
                              std::vector<T> &conf_nu, int size_mu1,
                              int size_mu2, int size_nu1, int size_nu2,
                              double alpha, double divisor) {
  int data_size1 = x_size * y_size * z_size * t_size;

  T bracket;

#pragma omp parallel for collapse(3) private(bracket)
  for (int k = 0; k < data_size1; k += size_mu2) {
    for (int i = 0; i < size_mu2; i += size_nu2) {
      for (int j = 0; j < size_nu2; j++) {
        if (j < size_nu2 - size_nu1)
          bracket = conf_nu[i + k + j] * conf_mu[i + k + j + size_nu1];
        else
          bracket =
              conf_nu[i + k + j] * conf_mu[i + k + j - size_nu2 + size_nu1];
        if (i < size_mu2 - size_mu1)
          smeared[i + k + j] =
              smeared[i + k + j] +
              (alpha / divisor) * (bracket ^ conf_nu[i + k + j + size_mu1]);
        else
          smeared[i + k + j] =
              smeared[i + k + j] +
              (alpha / divisor) *
                  (bracket ^ conf_nu[i + k + j - size_mu2 + size_mu1]);

        if (j >= size_nu1) {
          bracket =
              conf_nu[i + k + j - size_nu1] % conf_mu[i + k + j - size_nu1];
          if (i < size_mu2 - size_mu1) {
            smeared[i + k + j] =
                smeared[i + k + j] +
                (alpha / divisor) *
                    (bracket * conf_nu[i + k + j + size_mu1 - size_nu1]);

          } else
            smeared[i + k + j] =
                smeared[i + k + j] +
                (alpha / divisor) *
                    (bracket *
                     conf_nu[i + k + j - size_mu2 + size_mu1 - size_nu1]);
        } else {
          bracket = conf_nu[i + k + j + size_nu2 - size_nu1] %
                    conf_mu[i + k + j + size_nu2 - size_nu1];
          if (i < size_mu2 - size_mu1)
            smeared[i + k + j] =
                smeared[i + k + j] +
                (alpha / divisor) *
                    (bracket *
                     conf_nu[i + k + j + size_nu2 - size_nu1 + size_mu1]);
          else
            smeared[i + k + j] =
                smeared[i + k + j] +
                (alpha / divisor) *
                    (bracket * conf_nu[i + k + j + size_nu2 - size_nu1 -
                                       size_mu2 + size_mu1]);
        }
      }
    }
  }
}

std::map<std::tuple<int, int, int>, int> indices_map_HYP1() {
  std::map<std::tuple<int, int, int>, int> indices_map;
  int index_counter = 0;
  for (int nu = 0; nu < 2; nu++) {
    for (int sigma = nu + 1; sigma < 3; sigma++) {
      indices_map[std::tuple<int, int, int>(3, nu, sigma)] = index_counter;
      indices_map[std::tuple<int, int, int>(3, sigma, nu)] = index_counter;
      index_counter++;
    }
  }

  for (int nu = 0; nu < 3; nu++) {
    for (int sigma = 0; sigma < 3; sigma++) {
      if (nu != sigma) {
        indices_map[std::tuple<int, int, int>(nu, 3, sigma)] = index_counter;
        indices_map[std::tuple<int, int, int>(nu, sigma, 3)] = index_counter;
        index_counter++;
      }
    }
  }

  return indices_map;
}

std::map<std::tuple<int, int>, int> indices_map_HYP2() {
  std::map<std::tuple<int, int>, int> indices_map;
  int index_counter = 0;
  for (int nu = 0; nu < 3; nu++) {
    indices_map[std::tuple<int, int>(3, nu)] = index_counter;
    index_counter++;
  }

  for (int nu = 0; nu < 3; nu++) {
    indices_map[std::tuple<int, int>(nu, 3)] = index_counter;
    index_counter++;
  }

  return indices_map;
}

void show_mem_use() {
  long long totalPhysMem;
  long long physMemUsed;
  long long totalVirtualMem;
  long long virtualMemUsed;

  sysinfo(&memInfo);
  totalPhysMem = memInfo.totalram;
  totalPhysMem *= memInfo.mem_unit;
  std::cout << "total RAM1: " << totalPhysMem << std::endl;
  physMemUsed = memInfo.totalram - memInfo.freeram;
  physMemUsed *= memInfo.mem_unit;
  std::cout << "used RAM1: " << physMemUsed << std::endl;
  totalVirtualMem = memInfo.totalram;
  totalVirtualMem += memInfo.totalswap;
  totalVirtualMem *= memInfo.mem_unit;
  std::cout << "virtual RAM1: " << totalVirtualMem << std::endl;
  virtualMemUsed = memInfo.totalram - memInfo.freeram;
  virtualMemUsed += memInfo.totalswap - memInfo.freeswap;
  virtualMemUsed *= memInfo.mem_unit;
  std::cout << "virtual RAM1 used: " << virtualMemUsed << std::endl;
}

template <class T>
void smearing_HYP_new(std::vector<std::vector<T>> &conf, double alpha1,
                      double alpha2, double alpha3) {

  double start_time;
  double end_time;
  double calculation_time;

  std::vector<int> steps = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  std::map<std::tuple<int, int, int>, int> indices_map1 = indices_map_HYP1();
  std::vector<std::vector<T>> links1(9);

  std::cout << "before first step" << std::endl;
  show_mem_use();

  // first step

  start_time = omp_get_wtime();

  links1[0] = conf[3];
  std::for_each(links1[0].begin(), links1[0].end(),
                [alpha3](T &A) { A = (1 - alpha3) * A; });
  links1[1] = links1[0];
  links1[2] = links1[0];

  bool if_first_time = true;
  int index_tmp;

  for (int nu = 0; nu < 3; nu++) {
    if_first_time = true;
    for (int sigma = 0; sigma < 3; sigma++) {
      if (nu != sigma) {
        if (if_first_time) {
          index_tmp = indices_map1[std::tuple<int, int, int>(nu, 3, sigma)];
          links1[index_tmp] = conf[nu];
          std::for_each(links1[index_tmp].begin(), links1[index_tmp].end(),
                        [alpha3](T &A) { A = (1 - alpha3) * A; });
          if_first_time = false;
        } else
          links1[indices_map1[std::tuple<int, int, int>(nu, 3, sigma)]] =
              links1[index_tmp];
      }
    }
  }

  end_time = omp_get_wtime();
  calculation_time = end_time - start_time;
  std::cout << "creating vectors 1: " << calculation_time << std::endl;

  std::cout << "matrix created" << std::endl;
  std::cout << links1[indices_map1[std::tuple<int, int, int>(1, 0, 3)]][0]
            << std::endl;

  start_time = omp_get_wtime();

  for (int nu = 0; nu < 2; nu++) {
    for (int sigma = nu + 1; sigma < 3; sigma++) {
      for (int rho = 0; rho < 3; rho++) {
        if (rho != nu && rho != sigma) {
          smearing_plane_HYP_major(
              links1[indices_map1[std::tuple<int, int, int>(3, nu, sigma)]],
              conf[3], conf[rho], steps[3], steps[4], steps[rho],
              steps[rho + 1], alpha3, 2);
        }
      }
    }
  }
  std::cout << "before first step calculations" << std::endl;
  show_mem_use();

  for (int nu = 0; nu < 3; nu++) {
    for (int sigma = 0; sigma < 3; sigma++) {
      if (nu != sigma) {
        for (int rho = 0; rho < 3; rho++) {
          if (rho != nu && rho != sigma) {
            if (nu < rho)
              smearing_plane_HYP_minor(
                  links1[indices_map1[std::tuple<int, int, int>(nu, 3, sigma)]],
                  conf[nu], conf[rho], steps[nu], steps[nu + 1], steps[rho],
                  steps[rho + 1], alpha3, 2);
            if (nu > rho)
              smearing_plane_HYP_major(
                  links1[indices_map1[std::tuple<int, int, int>(nu, 3, sigma)]],
                  conf[nu], conf[rho], steps[nu], steps[nu + 1], steps[rho],
                  steps[rho + 1], alpha3, 2);
          }
        }
      }
    }
  }
  end_time = omp_get_wtime();
  calculation_time = end_time - start_time;
  std::cout << "calculations 1: " << calculation_time << std::endl;

  start_time = omp_get_wtime();

  for (int i = 0; i < 9; i++) {
#pragma omp parallel for firstprivate(i)
    for (int j = 0; j < links1[i].size(); j++) {
      links1[i][j] = links1[i][j].proj();
    }
  }
  end_time = omp_get_wtime();
  calculation_time = end_time - start_time;
  std::cout << "projection 1: " << calculation_time << std::endl;

  for (auto it = indices_map1.begin(); it != indices_map1.end(); it++) {
    std::cout << "matrix at index " << std::get<0>(it->first)
              << std::get<1>(it->first) << std::get<2>(it->first) << std::endl;
    std::cout << links1[it->second][0] << std::endl;
  }

  std::cout << "before second step" << std::endl;
  show_mem_use();

  // second_step

  start_time = omp_get_wtime();

  std::vector<std::vector<T>> links2(6);
  std::map<std::tuple<int, int>, int> indices_map2 = indices_map_HYP2();

  links2[0] = conf[3];
  std::for_each(links2[0].begin(), links2[0].end(),
                [alpha2](T &A) { A = (1 - alpha2) * A; });
  links2[1] = links2[0];
  links2[2] = links2[0];

  for (int nu = 0; nu < 3; nu++) {
    index_tmp = indices_map2[std::tuple<int, int>(nu, 3)];
    links2[index_tmp] = conf[nu];
    std::for_each(links2[index_tmp].begin(), links2[index_tmp].end(),
                  [alpha2](T &A) { A = (1 - alpha2) * A; });
  }
  end_time = omp_get_wtime();
  calculation_time = end_time - start_time;
  std::cout << "vectors creation 2: " << calculation_time << std::endl;

  start_time = omp_get_wtime();

  for (int nu = 0; nu < 3; nu++) {
    for (int sigma = 0; sigma < 3; sigma++) {
      if (sigma != nu)
        smearing_plane_HYP_major(
            links2[indices_map2[std::tuple<int, int>(3, nu)]],
            links1[indices_map1[std::tuple<int, int, int>(3, nu, sigma)]],
            links1[indices_map1[std::tuple<int, int, int>(sigma, 3, nu)]],
            steps[3], steps[4], steps[sigma], steps[sigma + 1], alpha2, 4);
    }
  }

  std::cout << "before second step calculations" << std::endl;
  show_mem_use();

  for (int nu = 0; nu < 3; nu++) {
    for (int sigma = 0; sigma < 3; sigma++) {
      if (sigma != nu) {
        if (nu < sigma)
          smearing_plane_HYP_minor(
              links2[indices_map2[std::tuple<int, int>(nu, 3)]],
              links1[indices_map1[std::tuple<int, int, int>(nu, 3, sigma)]],
              links1[indices_map1[std::tuple<int, int, int>(sigma, nu, 3)]],
              steps[nu], steps[nu + 1], steps[sigma], steps[sigma + 1], alpha2,
              4);
        if (nu > sigma)
          smearing_plane_HYP_major(
              links2[indices_map2[std::tuple<int, int>(nu, 3)]],
              links1[indices_map1[std::tuple<int, int, int>(nu, 3, sigma)]],
              links1[indices_map1[std::tuple<int, int, int>(sigma, nu, 3)]],
              steps[nu], steps[nu + 1], steps[sigma], steps[sigma + 1], alpha2,
              4);
      }
    }
  }

  end_time = omp_get_wtime();
  calculation_time = end_time - start_time;
  std::cout << "calculations 2: " << calculation_time << std::endl;

  start_time = omp_get_wtime();

  for (int i = 0; i < 6; i++) {
#pragma omp parallel for firstprivate(i)
    for (int j = 0; j < links2[i].size(); j++) {
      links2[i][j] = links2[i][j].proj();
    }
  }
  end_time = omp_get_wtime();
  calculation_time = end_time - start_time;
  std::cout << "projection 2: " << calculation_time << std::endl;
  // for (int nu = 0; nu < 3; nu++) {
  //   std::cout << links2[indices_map2[std::tuple<int, int>(nu, 3)]][0]
  //             << std::endl;
  //   std::cout << links2[indices_map2[std::tuple<int, int>(3, nu)]][0]
  //             << std::endl;
  // }
  std::cout << "before third step" << std::endl;
  show_mem_use();

  start_time = omp_get_wtime();

  std::vector<T> smeared = conf[3];

  std::for_each(smeared.begin(), smeared.end(),
                [alpha1](T &A) { A = (1 - alpha1) * A; });

  end_time = omp_get_wtime();
  calculation_time = end_time - start_time;
  std::cout << "creaction 3: " << calculation_time << std::endl;

  start_time = omp_get_wtime();

  for (int nu = 0; nu < 3; nu++) {
    smearing_plane_HYP_major(
        smeared, links2[indices_map2[std::tuple<int, int>(3, nu)]],
        links2[indices_map2[std::tuple<int, int>(nu, 3)]], steps[3], steps[4],
        steps[nu], steps[nu + 1], alpha1, 6);
  }

  end_time = omp_get_wtime();
  calculation_time = end_time - start_time;
  std::cout << "calculations 3: " << calculation_time << std::endl;

  start_time = omp_get_wtime();

#pragma omp parallel for
  for (int i = 0; i < smeared.size(); i++) {
    smeared[i] = smeared[i].proj();
  }

  end_time = omp_get_wtime();
  calculation_time = end_time - start_time;
  std::cout << "projection 3: " << calculation_time << std::endl;

  conf[3] = smeared;
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
void make_step1(std::vector<std::vector<T>> &links1,
                std::vector<std::vector<T>> &conf,
                std::vector<std::tuple<int, int, int>> &indices,
                std::map<std::tuple<int, int, int>, int> &indices_map,
                double alpha) {
  std::vector<int> steps = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};
  // for (int i = 0; i < indices.size(); i++) {
  //   std::cout << "indices full " << std::get<0>(indices[i])
  //             << std::get<1>(indices[i]) << std::get<2>(indices[i])
  //             << std::endl;
  // }
  for (int i = 0; i < indices.size(); i++) {
    // std::cout << "indices " << i << std::endl;
    // std::cout << "index " << std::get<0>(indices[i]) <<
    // std::get<1>(indices[i])
    //           << std::get<2>(indices[i]);
    // if (indices_map.count(indices[i]) == 0) {
    //   std::cout << " not present" << std::endl;
    // } else {
    //   std::cout << " present" << std::endl;
    // }
    if (indices_map.count(indices[i]) == 0) {
      // std::cout << "calculate indices " << std::get<0>(indices[i])
      //           << std::get<1>(indices[i]) << std::get<2>(indices[i])
      //           << std::endl;
      indices_map[indices[i]] = links1.size();
      links1.push_back(
          HYP_initialize_vector(conf[std::get<0>(indices[i])], alpha));
      // std::cout << "created matrix" << std::endl;
      // std::cout << links1.back()[0] << std::endl;
      for (int eta = 0; eta < 3; eta++) {
        // std::cout << "eta " << eta << std::endl;
        if (eta != std::get<0>(indices[i]) && eta != std::get<1>(indices[i]) &&
            eta != std::get<2>(indices[i])) {
          // std::cout << "calculate with eta " << eta << " indices "
          //           << std::get<0>(indices[i]) << std::get<1>(indices[i])
          //           << std::get<2>(indices[i]) << std::endl;
          if (eta > std::get<0>(indices[i])) {
            // std::cout << "plane minor" << std::endl;
            smearing_plane_HYP_minor(links1.back(),
                                     conf[std::get<0>(indices[i])], conf[eta],
                                     steps[std::get<0>(indices[i])],
                                     steps[std::get<0>(indices[i]) + 1],
                                     steps[eta], steps[eta + 1], alpha, 2);
          }
          if (eta < std::get<0>(indices[i])) {
            // std::cout << "plane major" << std::endl;
            smearing_plane_HYP_major(links1.back(),
                                     conf[std::get<0>(indices[i])], conf[eta],
                                     steps[std::get<0>(indices[i])],
                                     steps[std::get<0>(indices[i]) + 1],
                                     steps[eta], steps[eta + 1], alpha, 2);
          }
        }
      }
      // std::cout << "step1 projection" << std::endl;
#pragma omp parallel for
      for (int i = 0; i < links1.back().size(); i++) {
        links1.back()[i] = links1.back()[i].proj();
      }
    }
  }
  // for (auto it = indices_map.begin(); it != indices_map.end(); it++) {
  //   std::cout << "matrix of index " << std::get<0>(it->first)
  //             << std::get<1>(it->first) << std::get<2>(it->first) <<
  //             std::endl;
  //   std::cout << links1[it->second][0] << std::endl;
  // }
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
        // if (indices_map1.count(std::tuple<int, int, int>(3, nu, rho)) == 0) {
        //   std::cout << "problem: index not present " << 3 << nu << rho
        //             << std::endl;
        // }
      } else {
        index = indices_map1[std::tuple<int, int, int>(3, rho, nu)];
        // if (indices_map1.count(std::tuple<int, int, int>(3, rho, nu)) == 0) {
        //   std::cout << "problem: index not present " << 3 << rho << nu
        //             << std::endl;
        // }
      }
      // if (indices_map1.count(std::tuple<int, int, int>(rho, nu, 3)) == 0) {
      //   std::cout << "problem: index not present " << rho << nu << 3
      //             << std::endl;
      // }
      smearing_plane_HYP_major(
          links2[0], links1[index],
          links1[indices_map1[std::tuple<int, int, int>(rho, nu, 3)]], steps[3],
          steps[4], steps[rho], steps[rho + 1], alpha2, 4);
    }
  }
#pragma omp parallel for
  for (int i = 0; i < links2[0].size(); i++) {
    links2[0][i] = links2[0][i].proj();
  }
  for (int rho = 0; rho < 3; rho++) {
    if (rho != nu) {
      if (nu < rho) {
        // if (indices_map1.count(std::tuple<int, int, int>(nu, rho, 3)) == 0) {
        //   std::cout << "problem: index not present " << nu << rho << 3
        //             << std::endl;
        // }
        // if (indices_map1.count(std::tuple<int, int, int>(rho, nu, 3)) == 0) {
        //   std::cout << "problem: index not present " << rho << nu << 3
        //             << std::endl;
        // }
        smearing_plane_HYP_minor(
            links2[1],
            links1[indices_map1[std::tuple<int, int, int>(nu, rho, 3)]],
            links1[indices_map1[std::tuple<int, int, int>(rho, nu, 3)]],
            steps[nu], steps[nu + 1], steps[rho], steps[rho + 1], alpha2, 4);
      }
      if (nu > rho) {
        // if (indices_map1.count(std::tuple<int, int, int>(nu, rho, 3)) == 0) {
        //   std::cout << "problem: index not present " << nu << rho << 3
        //             << std::endl;
        // }
        // if (indices_map1.count(std::tuple<int, int, int>(rho, nu, 3)) == 0) {
        //   std::cout << "problem: index not present " << rho << nu << 3
        //             << std::endl;
        // }
        smearing_plane_HYP_major(
            links2[1],
            links1[indices_map1[std::tuple<int, int, int>(nu, rho, 3)]],
            links1[indices_map1[std::tuple<int, int, int>(rho, nu, 3)]],
            steps[nu], steps[nu + 1], steps[rho], steps[rho + 1], alpha2, 4);
      }
    }
  }
#pragma omp parallel for
  for (int i = 0; i < links2[1].size(); i++) {
    links2[1][i] = links2[1][i].proj();
  }
  // std::cout << "step2 matrices" << std::endl;
  // std::cout << links2[0][0] << std::endl;
  // std::cout << links2[1][0] << std::endl;
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
void smearing_HYP_parallel(std::vector<std::vector<T>> &conf, double alpha1,
                           double alpha2, double alpha3) {
  std::vector<int> steps = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};
  std::vector<std::vector<std::tuple<int, int, int>>> indices3 =
      make_indices3();
  // for (int i = 0; i < 3; i++) {
  // std::cout << "indices " << i << ":" << std::endl;
  // for (int j = 0; j < indices3[i].size(); j++) {
  //   std::cout << "indices full " << std::get<0>(indices3[i][j])
  //             << std::get<1>(indices3[i][j]) << std::get<2>(indices3[i][j])
  //             << std::endl;
  // }
  // }
  // std::cout << "initializing vectors" << std::endl;
  std::vector<T> smeared = HYP_initialize_vector(conf[3], alpha1);
  std::vector<std::vector<T>> links2(2);
  std::vector<std::vector<T>> links1;
  std::map<std::tuple<int, int, int>, int> indices_map3;
  std::vector<std::tuple<int, int, int>> deleted_indices;
  std::vector<std::tuple<int, int, int>> delete_indices;
  for (int nu3 = 0; nu3 < 3; nu3++) {
    std::cout << "before first step" << std::endl;
    show_mem_use();
    make_step1(links1, conf, indices3[nu3], indices_map3, alpha3);
    std::cout << "before second step" << std::endl;
    show_mem_use();
    make_step2(links2, links1, conf, indices_map3, nu3, alpha2);
    std::cout << "before deleting vectors" << std::endl;
    show_mem_use();
    delete_indices = indices_to_delete(indices3, deleted_indices, nu3);
    delete_vectors(delete_indices, indices_map3, links1);
    std::cout << "before third step" << std::endl;
    show_mem_use();
    smearing_plane_HYP_major(smeared, links2[0], links2[1], steps[3], steps[4],
                             steps[nu3], steps[nu3 + 1], alpha1, 6);
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
                       double alpha3); // staples for refreshing
                                       // algorythm(refresh link every step)
template su2 staples_third(const std::vector<std::vector<su2>> &smearing_second,
                           link1 &link, std::unordered_map<int, int> indexes,
                           int nu, int mu);
template su2 staples_third_refresh(const std::vector<su2> &vec, link1 &link,
                                   int eta, double alpha2, double alpha3);
template std::vector<su2> smearing_first(const std::vector<su2> &array,
                                         double alpha3, int mu, int nu,
                                         int rho);
template std::vector<std::vector<su2>>
smearing_first_full(const std::vector<su2> &array, double alpha3);
template std::vector<su2>
smearing_second(const std::vector<su2> &array,
                std::vector<std::vector<su2>> &smearing_first, double alpha2,
                int mu, int nu);
template std::vector<std::vector<su2>>
smearing_second_full(const std::vector<su2> &array,
                     std::vector<std::vector<su2>> &smearing_first,
                     double alpha2);
template std::vector<su2>
smearing_HYP(const std::vector<su2> &array,
             std::vector<std::vector<su2>> &smearing_second, double alpha1);
template std::vector<su2> smearing_APE(const std::vector<su2> &array,
                                       double alpha_APE);
template std::vector<su2> smearing1_APE(const std::vector<su2> &array,
                                        double alpha_APE);
template su2 smearing_first_refresh(const std::vector<su2> &vec, link1 &link,
                                    int nu, int rho,
                                    double alpha3); // refresh link every step
template su2 smearing_second_refresh(const std::vector<su2> &vec, link1 &link,
                                     int nu, double alpha2,
                                     double alpha3); // refresh link every step
template std::vector<su2>
smearing_HYP_refresh(data<su2> &conf, double alpha1, double alpha2,
                     double alpha3); // refresh link every step
template std::vector<su2>
smearing_APE_refresh(data<su2> &conf,
                     double alpha_APE); // refresh link every step

template std::vector<std::vector<su2>>
separate_smearing(std::vector<su2> &conf);

template void smearing_plane_minor(std::vector<su2> &smeared,
                                   const std::vector<su2> &conf_mu,
                                   const std::vector<su2> &conf_nu,
                                   int size_mu1, int size_mu2, int size_nu1,
                                   int size_nu2, double alpha);
template void smearing_plane_minor_start(std::vector<su2> &smeared,
                                         const std::vector<su2> &conf_mu,
                                         const std::vector<su2> &conf_nu,
                                         int size_mu1, int size_mu2,
                                         int size_nu1, int size_nu2,
                                         double alpha);

template void smearing_plane_major(std::vector<su2> &smeared,
                                   const std::vector<su2> &conf_mu,
                                   const std::vector<su2> &conf_nu,
                                   int size_mu1, int size_mu2, int size_nu1,
                                   int size_nu2, double alpha);
template void smearing_plane_major_end(std::vector<su2> &smeared,
                                       const std::vector<su2> &conf_mu,
                                       const std::vector<su2> &conf_nu,
                                       int size_mu1, int size_mu2, int size_nu1,
                                       int size_nu2, double alpha);

template void smearing_APE_parallel(std::vector<std::vector<su2>> &conf,
                                    double alpha);
template std::vector<std::vector<su2>>
smearing_APE_2d_initial(std::vector<std::vector<su2>> &conf, double alpha);
template void smearing_APE_2d(std::vector<std::vector<su2>> &conf,
                              double alpha);

template void smearing_plane_HYP_minor(std::vector<su2> &smeared,
                                       std::vector<su2> &conf_mu,
                                       std::vector<su2> &conf_nu, int size_mu1,
                                       int size_mu2, int size_nu1, int size_nu2,
                                       double alpha, double divisor);

template void smearing_plane_HYP_major(std::vector<su2> &smeared,
                                       std::vector<su2> &conf_mu,
                                       std::vector<su2> &conf_nu, int size_mu1,
                                       int size_mu2, int size_nu1, int size_nu2,
                                       double alpha, double divisor);

template void smearing_HYP_new(std::vector<std::vector<su2>> &conf,
                               double alpha1, double alpha2, double alpha3);
template void smearing_HYP_parallel(std::vector<std::vector<su2>> &conf,
                                    double alpha1, double alpha2,
                                    double alpha3);

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
                       double alpha3); // staples for refreshing
                                       // algorythm(refresh link every step)
template abelian
staples_third(const std::vector<std::vector<abelian>> &smearing_second,
              link1 &link, std::unordered_map<int, int> indexes, int nu,
              int mu);
template abelian staples_third_refresh(const std::vector<abelian> &vec,
                                       link1 &link, int eta, double alpha2,
                                       double alpha3);
template std::vector<abelian> smearing_first(const std::vector<abelian> &array,
                                             double alpha3, int mu, int nu,
                                             int rho);
template std::vector<std::vector<abelian>>
smearing_first_full(const std::vector<abelian> &array, double alpha3);
template std::vector<abelian>
smearing_second(const std::vector<abelian> &array,
                std::vector<std::vector<abelian>> &smearing_first,
                double alpha2, int mu, int nu);
template std::vector<std::vector<abelian>>
smearing_second_full(const std::vector<abelian> &array,
                     std::vector<std::vector<abelian>> &smearing_first,
                     double alpha2);
template std::vector<abelian>
smearing_HYP(const std::vector<abelian> &array,
             std::vector<std::vector<abelian>> &smearing_second, double alpha1);
template std::vector<abelian> smearing_APE(const std::vector<abelian> &array,
                                           double alpha_APE);
template std::vector<abelian> smearing1_APE(const std::vector<abelian> &array,
                                            double alpha_APE);
template abelian
smearing_first_refresh(const std::vector<abelian> &vec, link1 &link, int nu,
                       int rho,
                       double alpha3); // refresh link every step
template abelian
smearing_second_refresh(const std::vector<abelian> &vec, link1 &link, int nu,
                        double alpha2,
                        double alpha3); // refresh link every step
template std::vector<abelian>
smearing_HYP_refresh(data<abelian> &conf, double alpha1, double alpha2,
                     double alpha3); // refresh link every step
template std::vector<abelian>
smearing_APE_refresh(data<abelian> &conf,
                     double alpha_APE); // refresh link every step

template std::vector<std::vector<abelian>>
separate_smearing(std::vector<abelian> &conf);

template void smearing_plane_minor(std::vector<abelian> &smeared,
                                   const std::vector<abelian> &conf_mu,
                                   const std::vector<abelian> &conf_nu,
                                   int size_mu1, int size_mu2, int size_nu1,
                                   int size_nu2, double alpha);
template void smearing_plane_minor_start(std::vector<abelian> &smeared,
                                         const std::vector<abelian> &conf_mu,
                                         const std::vector<abelian> &conf_nu,
                                         int size_mu1, int size_mu2,
                                         int size_nu1, int size_nu2,
                                         double alpha);

template void smearing_plane_major(std::vector<abelian> &smeared,
                                   const std::vector<abelian> &conf_mu,
                                   const std::vector<abelian> &conf_nu,
                                   int size_mu1, int size_mu2, int size_nu1,
                                   int size_nu2, double alpha);
template void smearing_plane_major_end(std::vector<abelian> &smeared,
                                       const std::vector<abelian> &conf_mu,
                                       const std::vector<abelian> &conf_nu,
                                       int size_mu1, int size_mu2, int size_nu1,
                                       int size_nu2, double alpha);

template void smearing_APE_parallel(std::vector<std::vector<abelian>> &conf,
                                    double alpha);
template std::vector<std::vector<abelian>>
smearing_APE_2d_initial(std::vector<std::vector<abelian>> &conf, double alpha);
template void smearing_APE_2d(std::vector<std::vector<abelian>> &conf,
                              double alpha);

template void smearing_plane_HYP_minor(std::vector<abelian> &smeared,
                                       std::vector<abelian> &conf_mu,
                                       std::vector<abelian> &conf_nu,
                                       int size_mu1, int size_mu2, int size_nu1,
                                       int size_nu2, double alpha,
                                       double divisor);

template void smearing_plane_HYP_major(std::vector<abelian> &smeared,
                                       std::vector<abelian> &conf_mu,
                                       std::vector<abelian> &conf_nu,
                                       int size_mu1, int size_mu2, int size_nu1,
                                       int size_nu2, double alpha,
                                       double divisor);

template void smearing_HYP_new(std::vector<std::vector<abelian>> &conf,
                               double alpha1, double alpha2, double alpha3);
template void smearing_HYP_parallel(std::vector<std::vector<abelian>> &conf,
                                    double alpha1, double alpha2,
                                    double alpha3);

// su3
template su3 staples_first(const std::vector<su3> &vec, link1 &link, int eta);
template su3 staples_second(const std::vector<std::vector<su3>> &smearing_first,
                            link1 &link, std::unordered_map<int, int> &indexes,
                            int rho, int mu, int nu);
template su3
staples_second_refresh(const std::vector<su3> &vec, link1 &link, int eta,
                       int nu,
                       double alpha3); // staples for refreshing
                                       // algorythm(refresh link every step)
template su3 staples_third(const std::vector<std::vector<su3>> &smearing_second,
                           link1 &link, std::unordered_map<int, int> indexes,
                           int nu, int mu);
template su3 staples_third_refresh(const std::vector<su3> &vec, link1 &link,
                                   int eta, double alpha2, double alpha3);
template std::vector<su3> smearing_first(const std::vector<su3> &array,
                                         double alpha3, int mu, int nu,
                                         int rho);
template std::vector<std::vector<su3>>
smearing_first_full(const std::vector<su3> &array, double alpha3);
template std::vector<su3>
smearing_second(const std::vector<su3> &array,
                std::vector<std::vector<su3>> &smearing_first, double alpha2,
                int mu, int nu);
template std::vector<std::vector<su3>>
smearing_second_full(const std::vector<su3> &array,
                     std::vector<std::vector<su3>> &smearing_first,
                     double alpha2);
template std::vector<su3>
smearing_HYP(const std::vector<su3> &array,
             std::vector<std::vector<su3>> &smearing_second, double alpha1);
template std::vector<su3> smearing_APE(const std::vector<su3> &array,
                                       double alpha_APE);
template std::vector<su3> smearing1_APE(const std::vector<su3> &array,
                                        double alpha_APE);
template su3 smearing_first_refresh(const std::vector<su3> &vec, link1 &link,
                                    int nu, int rho,
                                    double alpha3); // refresh link every step
template su3 smearing_second_refresh(const std::vector<su3> &vec, link1 &link,
                                     int nu, double alpha2,
                                     double alpha3); // refresh link every step
template std::vector<su3>
smearing_HYP_refresh(data<su3> &conf, double alpha1, double alpha2,
                     double alpha3); // refresh link every step
template std::vector<su3>
smearing_APE_refresh(data<su3> &conf,
                     double alpha_APE); // refresh link every step

template std::vector<std::vector<su3>>
separate_smearing(std::vector<su3> &conf);

template void smearing_plane_minor(std::vector<su3> &smeared,
                                   const std::vector<su3> &conf_mu,
                                   const std::vector<su3> &conf_nu,
                                   int size_mu1, int size_mu2, int size_nu1,
                                   int size_nu2, double alpha);
template void smearing_plane_minor_start(std::vector<su3> &smeared,
                                         const std::vector<su3> &conf_mu,
                                         const std::vector<su3> &conf_nu,
                                         int size_mu1, int size_mu2,
                                         int size_nu1, int size_nu2,
                                         double alpha);

template void smearing_plane_major(std::vector<su3> &smeared,
                                   const std::vector<su3> &conf_mu,
                                   const std::vector<su3> &conf_nu,
                                   int size_mu1, int size_mu2, int size_nu1,
                                   int size_nu2, double alpha);
template void smearing_plane_major_end(std::vector<su3> &smeared,
                                       const std::vector<su3> &conf_mu,
                                       const std::vector<su3> &conf_nu,
                                       int size_mu1, int size_mu2, int size_nu1,
                                       int size_nu2, double alpha);

template void smearing_APE_parallel(std::vector<std::vector<su3>> &conf,
                                    double alpha);

std::vector<std::vector<su3>>
smearing_APE_2d_initial(std::vector<std::vector<su3>> &conf, double alpha);
template void smearing_APE_2d(std::vector<std::vector<su3>> &conf,
                              double alpha);

template void smearing_plane_HYP_minor(std::vector<su3> &smeared,
                                       std::vector<su3> &conf_mu,
                                       std::vector<su3> &conf_nu, int size_mu1,
                                       int size_mu2, int size_nu1, int size_nu2,
                                       double alpha, double divisor);

template void smearing_plane_HYP_major(std::vector<su3> &smeared,
                                       std::vector<su3> &conf_mu,
                                       std::vector<su3> &conf_nu, int size_mu1,
                                       int size_mu2, int size_nu1, int size_nu2,
                                       double alpha, double divisor);

template void smearing_HYP_new(std::vector<std::vector<su3>> &conf,
                               double alpha1, double alpha2, double alpha3);
template void smearing_HYP_parallel(std::vector<std::vector<su3>> &conf,
                                    double alpha1, double alpha2,
                                    double alpha3);

// su3_abelian
template su3_abelian staples_first(const std::vector<su3_abelian> &vec,
                                   link1 &link, int eta);
template su3_abelian
staples_second(const std::vector<std::vector<su3_abelian>> &smearing_first,
               link1 &link, std::unordered_map<int, int> &indexes, int rho,
               int mu, int nu);
template su3_abelian
staples_second_refresh(const std::vector<su3_abelian> &vec, link1 &link,
                       int eta, int nu,
                       double alpha3); // staples for refreshing
                                       // algorythm(refresh link every step)
template su3_abelian
staples_third(const std::vector<std::vector<su3_abelian>> &smearing_second,
              link1 &link, std::unordered_map<int, int> indexes, int nu,
              int mu);
template su3_abelian staples_third_refresh(const std::vector<su3_abelian> &vec,
                                           link1 &link, int eta, double alpha2,
                                           double alpha3);
template std::vector<su3_abelian>
smearing_first(const std::vector<su3_abelian> &array, double alpha3, int mu,
               int nu, int rho);
template std::vector<std::vector<su3_abelian>>
smearing_first_full(const std::vector<su3_abelian> &array, double alpha3);
template std::vector<su3_abelian>
smearing_second(const std::vector<su3_abelian> &array,
                std::vector<std::vector<su3_abelian>> &smearing_first,
                double alpha2, int mu, int nu);
template std::vector<std::vector<su3_abelian>>
smearing_second_full(const std::vector<su3_abelian> &array,
                     std::vector<std::vector<su3_abelian>> &smearing_first,
                     double alpha2);
template std::vector<su3_abelian>
smearing_HYP(const std::vector<su3_abelian> &array,
             std::vector<std::vector<su3_abelian>> &smearing_second,
             double alpha1);
template std::vector<su3_abelian>
smearing_APE(const std::vector<su3_abelian> &array, double alpha_APE);
template std::vector<su3_abelian>
smearing1_APE(const std::vector<su3_abelian> &array, double alpha_APE);
template su3_abelian
smearing_first_refresh(const std::vector<su3_abelian> &vec, link1 &link, int nu,
                       int rho,
                       double alpha3); // refresh link every step
template su3_abelian
smearing_second_refresh(const std::vector<su3_abelian> &vec, link1 &link,
                        int nu, double alpha2,
                        double alpha3); // refresh link every step
template std::vector<su3_abelian>
smearing_HYP_refresh(data<su3_abelian> &conf, double alpha1, double alpha2,
                     double alpha3); // refresh link every step
template std::vector<su3_abelian>
smearing_APE_refresh(data<su3_abelian> &conf,
                     double alpha_APE); // refresh link every step

template std::vector<std::vector<su3_abelian>>
separate_smearing(std::vector<su3_abelian> &conf);

template void smearing_plane_minor(std::vector<su3_abelian> &smeared,
                                   const std::vector<su3_abelian> &conf_mu,
                                   const std::vector<su3_abelian> &conf_nu,
                                   int size_mu1, int size_mu2, int size_nu1,
                                   int size_nu2, double alpha);
template void smearing_plane_minor_start(
    std::vector<su3_abelian> &smeared, const std::vector<su3_abelian> &conf_mu,
    const std::vector<su3_abelian> &conf_nu, int size_mu1, int size_mu2,
    int size_nu1, int size_nu2, double alpha);

template void smearing_plane_major(std::vector<su3_abelian> &smeared,
                                   const std::vector<su3_abelian> &conf_mu,
                                   const std::vector<su3_abelian> &conf_nu,
                                   int size_mu1, int size_mu2, int size_nu1,
                                   int size_nu2, double alpha);
template void smearing_plane_major_end(std::vector<su3_abelian> &smeared,
                                       const std::vector<su3_abelian> &conf_mu,
                                       const std::vector<su3_abelian> &conf_nu,
                                       int size_mu1, int size_mu2, int size_nu1,
                                       int size_nu2, double alpha);

template void smearing_APE_parallel(std::vector<std::vector<su3_abelian>> &conf,
                                    double alpha);

std::vector<std::vector<su3_abelian>>
smearing_APE_2d_initial(std::vector<std::vector<su3_abelian>> &conf,
                        double alpha);
template void smearing_APE_2d(std::vector<std::vector<su3_abelian>> &conf,
                              double alpha);

template void smearing_plane_HYP_minor(std::vector<su3_abelian> &smeared,
                                       std::vector<su3_abelian> &conf_mu,
                                       std::vector<su3_abelian> &conf_nu,
                                       int size_mu1, int size_mu2, int size_nu1,
                                       int size_nu2, double alpha,
                                       double divisor);

template void smearing_plane_HYP_major(std::vector<su3_abelian> &smeared,
                                       std::vector<su3_abelian> &conf_mu,
                                       std::vector<su3_abelian> &conf_nu,
                                       int size_mu1, int size_mu2, int size_nu1,
                                       int size_nu2, double alpha,
                                       double divisor);

template void smearing_HYP_new(std::vector<std::vector<su3_abelian>> &conf,
                               double alpha1, double alpha2, double alpha3);
template void smearing_HYP_parallel(std::vector<std::vector<su3_abelian>> &conf,
                                    double alpha1, double alpha2,
                                    double alpha3);