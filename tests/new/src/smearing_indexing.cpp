#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../include/indexing.h"

#include <omp.h>
#include <vector>

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

using namespace std;

template <class T>
std::vector<std::vector<T>> separate_smearing_unchanged(std::vector<T> &conf) {
  int data_size = x_size * y_size * z_size * t_size;
  std::vector<std::vector<T>> result(4, std::vector<T>(data_size));

  link1 link(x_size, y_size, z_size, t_size);

  for (int mu = 0; mu < 4; ++mu) {

    SPACE_ITER_START

    result[mu][link.place / 4] = conf[link.place + mu];

    SPACE_ITER_END
  }

  return result;
}

template <class T>
void smearing_plane_minor_test1(std::vector<T> &smeared,
                                std::vector<T> &conf_mu,
                                std::vector<T> &conf_nu, int size_mu1,
                                int size_mu2, int size_nu1, int size_nu2,
                                double alpha) {
  int data_size = x_size * y_size * z_size * t_size;

  T bracket;

#pragma omp parallel for collapse(3) private(bracket) num_threads(4)
  for (int k = 0; k < data_size; k += size_nu2) {
    for (int i = 0; i < size_nu2; i += size_mu2) {
      for (int j = 0; j < size_mu2; j++) {
        if (i < size_nu2 - size_nu1) {
          // if (size_mu1 == 1 && size_nu1 == x_size && k == 0 && i == 0 &&
          //     j == 0) {
          //   std::cout << "smearing_plane_minor A1 " << conf_nu[i + k + j]
          //             << std::endl;
          // }
          // if (size_mu1 == 1 && size_nu1 == x_size && k == 0 && i == 0 &&
          //     j == 0) {
          //   std::cout << "smearing_plane_minor A2 "
          //             << conf_mu[i + k + j + size_nu1] << std::endl;
          // }
          bracket = conf_nu[i + k + j] * conf_mu[i + k + j + size_nu1];
        } else
          bracket =
              conf_nu[i + k + j] * conf_mu[i + k + j - size_nu2 + size_nu1];
        if (j < size_mu2 - size_mu1) {

          // if (size_mu1 == 1 && size_nu1 == x_size && k == 0 && i == 0 &&
          //     j == 0) {
          //   std::cout << "smearing_plane_minor A3 "
          //             << conf_nu[i + k + j + size_mu1] << std::endl;
          // }

          smeared[i + k + j] =
              smeared[i + k + j] +
              alpha * (bracket ^ conf_nu[i + k + j + size_mu1]);
        } else
          smeared[i + k + j] =
              smeared[i + k + j] +
              alpha * (bracket ^ conf_nu[i + k + j - size_mu2 + size_mu1]);

        // if (size_mu1 == 1 && size_nu1 == x_size && k == 0 && i == 0 && j ==
        // 0) {
        //   std::cout << "smearing_plane_minor bracket test1 "
        //             << (bracket ^ conf_nu[i + k + j + size_mu1]) <<
        //             std::endl;
        // }

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
        // if (size_mu1 == 1 && size_nu1 == x_size && k == 0 && i == 0 && j ==
        // 0) {
        //   std::cout << "smearing_plane_minor bracket test2 "
        //             << (bracket *
        //                 conf_nu[i + k + j + size_nu2 - size_nu1 + size_mu1])
        //             << std::endl;
        // }
      }
    }
  }
}

template <class T>
void smearing_plane_major_test1(std::vector<T> &smeared,
                                std::vector<T> &conf_mu,
                                std::vector<T> &conf_nu, int size_mu1,
                                int size_mu2, int size_nu1, int size_nu2,
                                double alpha) {
  int data_size = x_size * y_size * z_size * t_size;

  T bracket;

#pragma omp parallel for collapse(3) private(bracket) num_threads(4)
  for (int k = 0; k < data_size; k += size_mu2) {
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
void smearing_APE_test1(std::vector<std::vector<T>> &conf, double alpha) {
  int data_size = x_size * y_size * z_size * t_size;
  std::vector<std::vector<T>> smeared(3);

  double start_time;
  double end_time;
  double search_time;

  start_time = omp_get_wtime();

  for (int i = 0; i < 3; i++) {
    smeared[i] = conf[i];
  }

  smearing_plane_minor_test1(smeared[0], conf[0], conf[1], 1, x_size, x_size,
                             x_size * y_size, alpha);
  smearing_plane_minor_test1(smeared[0], conf[0], conf[2], 1, x_size,
                             x_size * y_size, x_size * y_size * z_size, alpha);
  smearing_plane_major_test1(smeared[1], conf[1], conf[0], x_size,
                             x_size * y_size, 1, x_size, alpha);
  smearing_plane_minor_test1(smeared[1], conf[1], conf[2], x_size,
                             x_size * y_size, x_size * y_size,
                             x_size * y_size * z_size, alpha);
  smearing_plane_major_test1(smeared[2], conf[2], conf[0], x_size * y_size,
                             x_size * y_size * z_size, 1, x_size, alpha);
  smearing_plane_major_test1(smeared[2], conf[2], conf[1], x_size * y_size,
                             x_size * y_size * z_size, x_size, x_size * y_size,
                             alpha);

  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < smeared[j].size(); i++) {
      smeared[j][i] = smeared[j][i].proj();
    }
  }

  for (int i = 0; i < 3; i++) {
    conf[i] = smeared[i];
  }

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "smearing_APE_test1 time: " << search_time << std::endl;
}

// abelian
template std::vector<std::vector<abelian>>
separate_smearing_unchanged(std::vector<abelian> &conf);

template void smearing_plane_minor_test1(std::vector<abelian> &smeared,
                                         std::vector<abelian> &conf_mu,
                                         std::vector<abelian> &conf_nu,
                                         int size_mu1, int size_mu2,
                                         int size_nu1, int size_nu2,
                                         double alpha);

template void smearing_plane_major_test1(std::vector<abelian> &smeared,
                                         std::vector<abelian> &conf_mu,
                                         std::vector<abelian> &conf_nu,
                                         int size_mu1, int size_mu2,
                                         int size_nu1, int size_nu2,
                                         double alpha);

template void smearing_APE_test1(std::vector<std::vector<abelian>> &conf,
                                 double alpha);

// su2
template std::vector<std::vector<su2>>
separate_smearing_unchanged(std::vector<su2> &conf);

template void smearing_plane_minor_test1(std::vector<su2> &smeared,
                                         std::vector<su2> &conf_mu,
                                         std::vector<su2> &conf_nu,
                                         int size_mu1, int size_mu2,
                                         int size_nu1, int size_nu2,
                                         double alpha);

template void smearing_plane_major_test1(std::vector<su2> &smeared,
                                         std::vector<su2> &conf_mu,
                                         std::vector<su2> &conf_nu,
                                         int size_mu1, int size_mu2,
                                         int size_nu1, int size_nu2,
                                         double alpha);

template void smearing_APE_test1(std::vector<std::vector<su2>> &conf,
                                 double alpha);

// su3
template std::vector<std::vector<su3>>
separate_smearing_unchanged(std::vector<su3> &conf);

template void smearing_plane_minor_test1(std::vector<su3> &smeared,
                                         std::vector<su3> &conf_mu,
                                         std::vector<su3> &conf_nu,
                                         int size_mu1, int size_mu2,
                                         int size_nu1, int size_nu2,
                                         double alpha);

template void smearing_plane_major_test1(std::vector<su3> &smeared,
                                         std::vector<su3> &conf_mu,
                                         std::vector<su3> &conf_nu,
                                         int size_mu1, int size_mu2,
                                         int size_nu1, int size_nu2,
                                         double alpha);

template void smearing_APE_test1(std::vector<std::vector<su3>> &conf,
                                 double alpha);