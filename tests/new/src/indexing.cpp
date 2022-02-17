#include "../include/indexing.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"

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

template <class T>
double plaket_plane(std::vector<T> &conf_mu, std::vector<T> &conf_nu, int size1,
                    int size2) {
  int data_size = x_size * y_size * z_size * t_size;
  std::vector<T> plakets(data_size);
  for (int i = 0; i < data_size; i += t_size) {
    for (int j = i; j < i + t_size - 1; j++) {
      plakets[j] = conf_mu[j] * conf_nu[j + 1];
    }
    plakets[i + t_size - 1] = conf_mu[i + t_size - 1] * conf_nu[i];
  }
  for (int i = 0; i < data_size; i += size2) {
    for (int j = i; j < i + size2 - size1; j++) {
      plakets[j] = plakets[j] * conf_mu[j + size1].conj();
    }
    for (int j = i; j < i + size1; j++) {
      plakets[j + size2 - size1] =
          plakets[j + size2 - size1] * conf_mu[j].conj();
    }
  }
  double result = 0;
  for (int i = 0; i < data_size; i++) {
    result += (plakets[i] * conf_nu[i].conj()).tr();
  }

  return result / data_size;
}

template <class T>
double plaket_time_test(std::vector<std::vector<T>> &separated) {

  double result = 0;
  result += plaket_plane(separated[3], separated[0], t_size, x_size * t_size);
  result += plaket_plane(separated[3], separated[1], x_size * t_size,
                         y_size * x_size * t_size);
  result += plaket_plane(separated[3], separated[2], y_size * x_size * t_size,
                         y_size * x_size * t_size * z_size);
  return result / 3;
}

template <class T>
double plaket_time_test1(std::vector<std::vector<T>> &separated) {
  int data_size = x_size * y_size * z_size * t_size;
  std::vector<std::vector<T>> plakets(3, std::vector<T>(data_size));
  for (int i = 0; i < data_size; i += t_size) {
    for (int j = i; j < i + t_size - 1; j++) {
      for (int mu = 0; mu < 3; mu++) {
        plakets[mu][j] = separated[3][j] * separated[mu][j + 1];
      }
    }
    for (int mu = 0; mu < 3; mu++) {
      plakets[mu][i + t_size - 1] =
          separated[mu][i + t_size - 1] * separated[mu][i];
    }
  }
  std::vector<int> sizes1 = {t_size, x_size * t_size, y_size * x_size * t_size};
  std::vector<int> sizes2 = {x_size * t_size, y_size * x_size * t_size,
                             y_size * x_size * t_size * z_size};

  for (int mu = 0; mu < 3; mu++) {
    for (int i = 0; i < data_size; i += sizes2[mu]) {
      for (int j = i; j < i + sizes2[mu] - sizes1[mu]; j++) {
        plakets[mu][j] = plakets[mu][j] * separated[mu][j + sizes1[mu]].conj();
      }
      for (int j = i; j < i + sizes1[mu]; j++) {
        plakets[mu][j + sizes2[mu] - sizes1[mu]] =
            plakets[mu][j + sizes2[mu] - sizes1[mu]] * separated[mu][j].conj();
      }
    }
  }
  double result = 0;
  for (int i = 0; i < data_size; i++) {
    for (int mu = 0; mu < 3; mu++) {
      result += (plakets[mu][i] * separated[mu][i].conj()).tr();
    }
  }

  return result / (3 * data_size);
}

template <class T>
std::vector<std::vector<T>> separate_3(std::vector<T> &conf, int mu) {
  int data_size = x_size * y_size * z_size * t_size;
  std::vector<std::vector<T>> result(4, std::vector<T>(data_size));

  link1 link(x_size, y_size, z_size, t_size);

  for (int mu = 0; mu < 4; ++mu) {
    SPACE_ITER_START
    result[mu]
          [link.coordinate[2] * link.lattice_size[3] * link.lattice_size[0] *
               link.lattice_size[1] +
           link.coordinate[1] * link.lattice_size[3] * link.lattice_size[0] +
           link.coordinate[0] * link.lattice_size[3] + link.coordinate[3]] =
              conf[link.place + mu];

    SPACE_ITER_END
  }

  return result;
}

// su2
template double plaket_plane(std::vector<su2> &conf_mu,
                             std::vector<su2> &conf_nu, int size1, int size2);

template double plaket_time_test(std::vector<std::vector<su2>> &separated);
template double plaket_time_test1(std::vector<std::vector<su2>> &separated);

template std::vector<std::vector<su2>> separate_3(std::vector<su2> &conf,
                                                  int mu);

// abelian

template double plaket_plane(std::vector<abelian> &conf_mu,
                             std::vector<abelian> &conf_nu, int size1,
                             int size2);

template double plaket_time_test(std::vector<std::vector<abelian>> &separated);
template double plaket_time_test1(std::vector<std::vector<abelian>> &separated);

template std::vector<std::vector<abelian>>
separate_3(std::vector<abelian> &conf, int mu);

// su3
template double plaket_plane(std::vector<su3_full> &conf_mu,
                             std::vector<su3_full> &conf_nu, int size1,
                             int size2);

template double plaket_time_test(std::vector<std::vector<su3_full>> &separated);
template double
plaket_time_test1(std::vector<std::vector<su3_full>> &separated);

template std::vector<std::vector<su3_full>>
separate_3(std::vector<su3_full> &conf, int mu);
