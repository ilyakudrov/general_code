#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../include/indexing.h"

#include <algorithm>
#include <execution>
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

template <class T>
void smearing_plane_HYP_minor_test1(std::vector<T> &smeared,
                                    std::vector<T> &conf_mu,
                                    std::vector<T> &conf_nu, int size_mu1,
                                    int size_mu2, int size_nu1, int size_nu2,
                                    double alpha, double divisor) {
  int data_size = x_size * y_size * z_size * t_size;

  T bracket;

#pragma omp parallel for collapse(3) private(bracket) num_threads(4)
  for (int k = 0; k < data_size; k += size_nu2) {
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
void smearing_plane_HYP_major_test1(std::vector<T> &smeared,
                                    std::vector<T> &conf_mu,
                                    std::vector<T> &conf_nu, int size_mu1,
                                    int size_mu2, int size_nu1, int size_nu2,
                                    double alpha, double divisor) {
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

std::map<std::tuple<int, int, int>, int> make_map_HYP1() {
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

std::map<std::tuple<int, int>, int> make_map_HYP2() {
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

template <class T>
void smearing_HYP_test1(std::vector<std::vector<T>> &conf, double alpha1,
                        double alpha2, double alpha3) {
  int data_size = x_size * y_size * z_size * t_size;

  double start_time;
  double end_time;
  double search_time;

  std::vector<int> steps = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  std::map<std::tuple<int, int, int>, int> indices_map1 = make_map_HYP1();
  std::vector<std::vector<T>> links1(9);

  // first step

  links1[0] = conf[3];
  std::for_each(links1[0].begin(), links1[0].end(),
                [alpha3](T &A) { A = (1 - alpha3) * A; });
  links1[1] = links1[0];
  links1[2] = links1[0];

  bool if_first_time = true;
  int index_tmp;

  start_time = omp_get_wtime();

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

  for (int nu = 0; nu < 2; nu++) {
    for (int sigma = nu + 1; sigma < 3; sigma++) {
      for (int rho = 0; rho < 3; rho++) {
        if (rho != nu && rho != sigma) {
          smearing_plane_HYP_major_test1(
              links1[indices_map1[std::tuple<int, int, int>(3, nu, sigma)]],
              conf[3], conf[rho], steps[3], steps[4], steps[rho],
              steps[rho + 1], alpha3, 2);
        }
      }
    }
  }

  for (int nu = 0; nu < 3; nu++) {
    for (int sigma = 0; sigma < 3; sigma++) {
      if (nu != sigma) {
        for (int rho = 0; rho < 3; rho++) {
          if (rho != nu && rho != sigma) {
            if (nu < rho)
              smearing_plane_HYP_minor_test1(
                  links1[indices_map1[std::tuple<int, int, int>(nu, 3, sigma)]],
                  conf[nu], conf[rho], steps[nu], steps[nu + 1], steps[rho],
                  steps[rho + 1], alpha3, 2);
            if (nu > rho)
              smearing_plane_HYP_major_test1(
                  links1[indices_map1[std::tuple<int, int, int>(nu, 3, sigma)]],
                  conf[nu], conf[rho], steps[nu], steps[nu + 1], steps[rho],
                  steps[rho + 1], alpha3, 2);
          }
        }
      }
    }
  }

  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < links1[i].size(); j++) {
      links1[i][j] = links1[i][j].proj();
    }
  }

  // second_step

  std::vector<std::vector<T>> links2(6);
  std::map<std::tuple<int, int>, int> indices_map2 = make_map_HYP2();

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

  for (int nu = 0; nu < 3; nu++) {
    for (int sigma = 0; sigma < 3; sigma++) {
      if (sigma != nu)
        smearing_plane_HYP_major_test1(
            links2[indices_map2[std::tuple<int, int>(3, nu)]],
            links1[indices_map1[std::tuple<int, int, int>(3, nu, sigma)]],
            links1[indices_map1[std::tuple<int, int, int>(sigma, 3, nu)]],
            steps[3], steps[4], steps[sigma], steps[sigma + 1], alpha2, 4);
    }
  }

  for (int nu = 0; nu < 3; nu++) {
    for (int sigma = 0; sigma < 3; sigma++) {
      if (sigma != nu) {
        if (nu < sigma)
          smearing_plane_HYP_minor_test1(
              links2[indices_map2[std::tuple<int, int>(nu, 3)]],
              links1[indices_map1[std::tuple<int, int, int>(nu, 3, sigma)]],
              links1[indices_map1[std::tuple<int, int, int>(sigma, nu, 3)]],
              steps[nu], steps[nu + 1], steps[sigma], steps[sigma + 1], alpha2,
              4);
        if (nu > sigma)
          smearing_plane_HYP_major_test1(
              links2[indices_map2[std::tuple<int, int>(nu, 3)]],
              links1[indices_map1[std::tuple<int, int, int>(nu, 3, sigma)]],
              links1[indices_map1[std::tuple<int, int, int>(sigma, nu, 3)]],
              steps[nu], steps[nu + 1], steps[sigma], steps[sigma + 1], alpha2,
              4);
      }
    }
  }

  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < links2[i].size(); j++) {
      links2[i][j] = links2[i][j].proj();
    }
  }

  // third step

  std::vector<T> smeared = conf[3];

  std::for_each(smeared.begin(), smeared.end(),
                [alpha1](T &A) { A = (1 - alpha1) * A; });

  for (int nu = 0; nu < 3; nu++) {
    smearing_plane_HYP_major_test1(
        smeared, links2[indices_map2[std::tuple<int, int>(3, nu)]],
        links2[indices_map2[std::tuple<int, int>(nu, 3)]], steps[3], steps[4],
        steps[nu], steps[nu + 1], alpha1, 6);
  }

  for (int i = 0; i < smeared.size(); i++) {
    smeared[i] = smeared[i].proj();
  }

  conf[3] = smeared;

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "smearing_HYP_test1 time: " << search_time << std::endl;
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

template void smearing_plane_HYP_minor_test1(std::vector<abelian> &smeared,
                                             std::vector<abelian> &conf_mu,
                                             std::vector<abelian> &conf_nu,
                                             int size_mu1, int size_mu2,
                                             int size_nu1, int size_nu2,
                                             double alpha, double divisor);

template void smearing_plane_HYP_major_test1(std::vector<abelian> &smeared,
                                             std::vector<abelian> &conf_mu,
                                             std::vector<abelian> &conf_nu,
                                             int size_mu1, int size_mu2,
                                             int size_nu1, int size_nu2,
                                             double alpha, double divisor);

template void smearing_HYP_test1(std::vector<std::vector<abelian>> &conf,
                                 double alpha1, double alpha2, double alpha3);

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

template void smearing_plane_HYP_minor_test1(std::vector<su2> &smeared,
                                             std::vector<su2> &conf_mu,
                                             std::vector<su2> &conf_nu,
                                             int size_mu1, int size_mu2,
                                             int size_nu1, int size_nu2,
                                             double alpha, double divisor);

template void smearing_plane_HYP_major_test1(std::vector<su2> &smeared,
                                             std::vector<su2> &conf_mu,
                                             std::vector<su2> &conf_nu,
                                             int size_mu1, int size_mu2,
                                             int size_nu1, int size_nu2,
                                             double alpha, double divisor);

template void smearing_HYP_test1(std::vector<std::vector<su2>> &conf,
                                 double alpha1, double alpha2, double alpha3);

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

template void smearing_plane_HYP_minor_test1(std::vector<su3> &smeared,
                                             std::vector<su3> &conf_mu,
                                             std::vector<su3> &conf_nu,
                                             int size_mu1, int size_mu2,
                                             int size_nu1, int size_nu2,
                                             double alpha, double divisor);

template void smearing_plane_HYP_major_test1(std::vector<su3> &smeared,
                                             std::vector<su3> &conf_mu,
                                             std::vector<su3> &conf_nu,
                                             int size_mu1, int size_mu2,
                                             int size_nu1, int size_nu2,
                                             double alpha, double divisor);

template void smearing_HYP_test1(std::vector<std::vector<su3>> &conf,
                                 double alpha1, double alpha2, double alpha3);