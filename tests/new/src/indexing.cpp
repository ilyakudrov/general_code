#include "../include/indexing.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"

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
      plakets[j] = plakets[j] ^ &conf_mu[j + size1];
    }
    for (int j = i; j < i + size1; j++) {
      plakets[j + size2 - size1] = plakets[j + size2 - size1] ^ &conf_mu[j];
    }
  }
  double result = 0;
  for (int i = 0; i < data_size; i++) {
    result += (plakets[i] ^ &conf_nu[i]).tr();
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
          separated[3][i + t_size - 1] * separated[mu][i];
    }
  }
  std::vector<int> sizes1 = {t_size, x_size * t_size, y_size * x_size * t_size};
  std::vector<int> sizes2 = {x_size * t_size, y_size * x_size * t_size,
                             y_size * x_size * t_size * z_size};

  for (int mu = 0; mu < 3; mu++) {
    for (int i = 0; i < data_size; i += sizes2[mu]) {
      for (int j = i; j < i + sizes2[mu] - sizes1[mu]; j++) {
        plakets[mu][j] = plakets[mu][j] ^ &separated[3][j + sizes1[mu]];
      }
      for (int j = i; j < i + sizes1[mu]; j++) {
        plakets[mu][j + sizes2[mu] - sizes1[mu]] =
            plakets[mu][j + sizes2[mu] - sizes1[mu]] ^ &separated[3][j];
      }
    }
  }
  double result = 0;
  for (int i = 0; i < data_size; i++) {
    for (int mu = 0; mu < 3; mu++) {
      result += (plakets[mu][i] ^ &separated[mu][i]).tr();
    }
  }

  return result / (3 * data_size);
}

template <class T>
double plaket_time_test2(std::vector<std::vector<T>> &separated) {
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
          separated[3][i + t_size - 1] * separated[mu][i];
    }
  }
  std::vector<int> sizes1 = {t_size, x_size * t_size, y_size * x_size * t_size};
  std::vector<int> sizes2 = {x_size * t_size, y_size * x_size * t_size,
                             y_size * x_size * t_size * z_size};

  for (int i = 0; i < data_size; i += sizes2[0]) {
    for (int t = i; t < i + sizes2[0] - sizes1[0]; t++) {
      plakets[0][t] = plakets[0][t] ^ &separated[3][t + sizes1[0]];
      plakets[1][t] = plakets[1][t] ^ &separated[3][t + sizes1[1]];
      plakets[2][t] = plakets[2][t] ^ &separated[3][t + sizes1[2]];
    }
    for (int t = i + sizes2[0] - sizes1[0]; t < i + sizes2[0]; t++) {
      plakets[1][t] = plakets[1][t] ^ &separated[3][t + sizes1[1]];
      plakets[2][t] = plakets[2][t] ^ &separated[3][t + sizes1[2]];
    }
    for (int t = i; t < i + sizes1[0]; t++) {
      plakets[0][t + sizes2[0] - sizes1[0]] =
          plakets[0][t + sizes2[0] - sizes1[0]] ^ &separated[3][t];
    }
  }

  double result = 0;
  for (int i = 0; i < data_size; i++) {
    for (int mu = 0; mu < 3; mu++) {
      result += (plakets[mu][i] ^ &separated[mu][i]).tr();
    }
  }

  return result / (3 * data_size);
}

template <class T>
double plaket_plane1(std::vector<T> &conf_mu, std::vector<T> &conf_nu,
                     int size1, int size2) {
  int data_size = x_size * y_size * z_size * t_size;
  std::vector<T> plakets(data_size);
  double result = 0;
  for (int k = 0; k < data_size; k += size2) {
    for (int i = k; i < k + size2; i += t_size) {
      for (int j = i; j < i + t_size - 1; j++) {
        plakets[j] = conf_mu[j] * conf_nu[j + 1];
      }
      plakets[i + t_size - 1] = conf_mu[i + t_size - 1] * conf_nu[i];
    }
    for (int j = k; j < k + size2 - size1; j++) {
      plakets[j] = plakets[j] ^ &conf_mu[j + size1];
    }
    for (int j = k; j < k + size1; j++) {
      plakets[j + size2 - size1] = plakets[j + size2 - size1] ^ &conf_mu[j];
    }
    for (int i = k; i < k + size2; i++) {
      result += (plakets[i] ^ &conf_nu[i]).tr();
    }
  }

  return result / data_size;
}

template <class T>
double plaket_time_test3(std::vector<std::vector<T>> &separated) {

  double result = 0;
  result += plaket_plane1(separated[3], separated[0], t_size, x_size * t_size);
  result += plaket_plane1(separated[3], separated[1], x_size * t_size,
                          y_size * x_size * t_size);
  result += plaket_plane1(separated[3], separated[2], y_size * x_size * t_size,
                          y_size * x_size * t_size * z_size);
  return result / 3;
}

template <class T>
double plaket_plane2(std::vector<T> &conf_mu, std::vector<T> &conf_nu,
                     int size1, int size2) {
  int data_size = x_size * y_size * z_size * t_size;

  std::vector<T> plakets(size2);
  double result = 0;

  for (int k = 0; k < data_size; k += size2) {
    for (int i = 0; i < size2; i += t_size) {
      for (int j = 0; j < t_size - 1; j++) {
        plakets[j + i] = conf_mu[j + i + k] * conf_nu[j + i + k + 1];
      }
      plakets[i + t_size - 1] = conf_mu[i + k + t_size - 1] * conf_nu[i + k];
    }
    for (int j = 0; j < size2 - size1; j++) {
      plakets[j] = plakets[j] ^ &conf_mu[j + k + size1];
    }
    for (int j = 0; j < size1; j++) {
      plakets[j + size2 - size1] = plakets[j + size2 - size1] ^ &conf_mu[j + k];
    }
    for (int i = 0; i < size2; i++) {
      result += plakets[i].multiply_tr(&conf_nu[i + k]);
    }
  }

  return result / data_size;
}

template <class T>
double plaket_time_test4(std::vector<std::vector<T>> &separated) {

  double result = 0;
  result += plaket_plane2(separated[3], separated[0], t_size, x_size * t_size);
  result += plaket_plane2(separated[3], separated[1], x_size * t_size,
                          y_size * x_size * t_size);
  result += plaket_plane2(separated[3], separated[2], y_size * x_size * t_size,
                          y_size * x_size * t_size * z_size);
  return result / 3;
}

template <class T> double plaket_time_test5(std::vector<T> &conf) {
  link1 link(x_size, y_size, z_size, t_size);
  int data_size = x_size * y_size * z_size * t_size * 4;

  std::vector<T> plakets(x_size * 3);
  int place;

  double result = 0;

  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        link.coordinate[1] = y;
        link.coordinate[2] = z;
        link.coordinate[3] = t;
        link.update(1);
        link.update(2);
        link.update(3);
        for (int x = 0; x < x_size; x++) {
          link.coordinate[0] = x;
          link.update(0);

          place = link.place;
          for (int mu = 0; mu < 3; mu++) {
            link.move(mu, 1);
            plakets[x * 3 + mu] = conf[place + mu] * conf[link.place + 3];
            link.move(mu, -1);
          }
        }
        for (int x = 0; x < x_size; x++) {
          link.coordinate[0] = x;
          link.update(0);

          place = link.place;
          link.move(3, 1);
          for (int mu = 0; mu < 3; mu++) {
            plakets[x * 3 + mu] = plakets[x * 3 + mu] ^ &conf[link.place + mu];
          }
          link.move(3, -1);
        }
        for (int x = 0; x < x_size; x++) {
          link.coordinate[0] = x;
          link.update(0);

          for (int mu = 0; mu < 3; mu++) {
            result += plakets[x * 3 + mu].multiply_tr(&conf[link.place + 3]);
          }
        }
      }
    }
  }
  return result / (x_size * y_size * z_size * t_size * 3);
}

template <class T> double plaket_time_test6(std::vector<T> &conf) {
  link1 link(x_size, y_size, z_size, t_size);
  int data_size = x_size * y_size * z_size * t_size * 4;

  std::vector<T> plakets(x_size * 3);
  int place;

  double result = 0;

  for (int mu = 0; mu < 3; mu++) {
    for (int t = 0; t < t_size; t++) {
      for (int z = 0; z < z_size; z++) {
        for (int y = 0; y < y_size; y++) {
          link.coordinate[1] = y;
          link.coordinate[2] = z;
          link.coordinate[3] = t;
          link.update(1);
          link.update(2);
          link.update(3);
          for (int x = 0; x < x_size; x++) {
            link.coordinate[0] = x;
            link.update(0);

            place = link.place;
            link.move(mu, 1);
            plakets[x * 3 + mu] = conf[place + mu] * conf[link.place + 3];
            link.move(mu, -1);
          }
          for (int x = 0; x < x_size; x++) {
            link.coordinate[0] = x;
            link.update(0);

            place = link.place;
            link.move(3, 1);
            plakets[x * 3 + mu] = plakets[x * 3 + mu] ^ &conf[link.place + mu];
            link.move(3, -1);
          }
          for (int x = 0; x < x_size; x++) {
            link.coordinate[0] = x;
            link.update(0);

            result += plakets[x * 3 + mu].multiply_tr(&conf[link.place + 3]);
          }
        }
      }
    }
  }
  return result / (x_size * y_size * z_size * t_size * 3);
}

template <class T>
double plaket_plane3(std::vector<T> &conf, int mu, int nu, int size_mu1,
                     int size_mu2, int size_nu1, int size_nu2) {
  int data_size = x_size * y_size * z_size * t_size;

  std::vector<T> plakets(size_mu2);
  double result = 0;

  for (int i = 0; i < data_size; i += size_mu2) {
    for (int n = 0; n < size_mu2; n += size_nu2) {
      for (int j = 0; j < size_nu2 - size_nu1; ++j) {
        plakets[j + n] =
            conf[4 * (j + n + i) + nu] * conf[4 * (j + n + i + size_nu1) + mu];
      }
      for (int j = 0; j < size_nu1; ++j) {
        plakets[j + n + size_nu2 - size_nu1] =
            conf[4 * (j + n + i + size_nu2 - size_nu1) + nu] *
            conf[4 * (j + n + i) + mu];
      }
    }

    for (int j = 0; j < size_mu2 - size_mu1; ++j) {
      plakets[j] = plakets[j] ^ &conf[4 * (j + i + size_mu1) + nu];
    }
    for (int j = 0; j < size_mu1; ++j) {
      plakets[j + size_mu2 - size_mu1] =
          plakets[j + size_mu2 - size_mu1] ^ &conf[4 * (j + i) + nu];
    }

    for (int j = 0; j < size_mu2; ++j) {
      result += plakets[j].multiply_tr(&conf[4 * (j + i) + mu]);
    }
  }

  return result / data_size;
}

template <class T> double plaket_time_test7(std::vector<T> &conf) {

  double result = 0;
  result += plaket_plane3(conf, 3, 0, x_size * y_size * z_size,
                          x_size * y_size * z_size * t_size, 1, x_size);
  result +=
      plaket_plane3(conf, 3, 1, x_size * y_size * z_size,
                    x_size * y_size * z_size * t_size, x_size, x_size * y_size);
  result += plaket_plane3(conf, 3, 2, x_size * y_size * z_size,
                          x_size * y_size * z_size * t_size, x_size * y_size,
                          x_size * y_size * z_size);
  return result / 3;
}

template <class T>
double plaket_plane4(std::vector<T> &conf_mu, std::vector<T> &conf_nu,
                     int size_mu1, int size_mu2, int size_nu1, int size_nu2) {
  int data_size = x_size * y_size * z_size * t_size;

  std::vector<T> plakets(size_nu2);
  double result = 0;

  for (int k = 0; k < data_size; k += size_nu2) {
    for (int i = 0; i < size_nu2; i += size_mu2) {
      for (int j = 0; j < size_mu2 - size_mu1; j++) {
        plakets[j + i] = conf_mu[j + i + k] * conf_nu[j + i + k + size_mu1];
      }
      for (int j = 0; j < size_mu1; j++) {
        plakets[j + i + size_mu2 - size_mu1] =
            conf_mu[j + i + k + size_mu2 - size_mu1] * conf_nu[j + i + k];
      }
    }
    for (int j = 0; j < size_nu2 - size_nu1; j++) {
      plakets[j] = plakets[j] ^ &conf_mu[j + k + size_nu1];
    }
    for (int j = 0; j < size_nu1; j++) {
      plakets[j + size_nu2 - size_nu1] =
          plakets[j + size_nu2 - size_nu1] ^ &conf_mu[j + k];
    }
    for (int i = 0; i < size_nu2; i++) {
      result += plakets[i].multiply_tr(&conf_nu[i + k]);
    }
  }

  return result / data_size;
}

template <class T>
double plaket_time_test8(std::vector<std::vector<T>> &separated) {

  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  double result = 0;
  result += plaket_plane4(separated[3], separated[0], 1, t_size, t_size,
                          x_size * t_size);
  result += plaket_plane4(separated[3], separated[1], 1, t_size,
                          x_size * t_size, y_size * x_size * t_size);
  result += plaket_plane4(separated[3], separated[2], 1, t_size,
                          y_size * x_size * t_size,
                          y_size * x_size * t_size * z_size);
  result += plaket_plane4(separated[0], separated[1], t_size, x_size * t_size,
                          x_size * t_size, y_size * x_size * t_size);
  result += plaket_plane4(separated[0], separated[2], t_size, x_size * t_size,
                          y_size * x_size * t_size,
                          y_size * x_size * t_size * z_size);
  result += plaket_plane4(separated[1], separated[2], x_size * t_size,
                          y_size * x_size * t_size, y_size * x_size * t_size,
                          y_size * x_size * t_size * z_size);
  return result / 6;
}

template <class T>
double plaket_plane5(std::vector<T> &conf_mu, std::vector<T> &conf_nu,
                     int size_mu1, int size_mu2, int size_nu1, int size_nu2,
                     int step) {
  int data_size = x_size * y_size * z_size * t_size;

  std::vector<T> plakets(step);
  // T plakets;
  double result = 0;

  for (int k = 0; k < data_size; k += size_nu2) {
    for (int i = k; i < k + size_nu2; i += size_mu2) {
      for (int j = i; j < i + size_mu2; j += step) {
        for (int m = 0; m < step; m++) {
          if (j + m < i + size_mu2 - size_mu1)
            plakets[m] = conf_mu[j + m] * conf_nu[j + m + size_mu1];
          else
            plakets[m] = conf_mu[j + m] * conf_nu[j + m - size_mu2 + size_mu1];
          // }
          // for (int m = 0; m < step; m++) {
          if (j + m < k + size_nu2 - size_nu1)
            plakets[m] = plakets[m] ^ &conf_mu[j + m + size_nu1];
          else
            plakets[m] = plakets[m] ^ &conf_mu[j + m - size_nu2 + size_nu1];
          // }
          // for (int m = 0; m < step; m++) {
          result += plakets[m].multiply_tr(&conf_nu[j + m]);
        }
      }
    }
  }

  return result / data_size;
}

template <class T>
double plaket_time_test9(std::vector<std::vector<T>> &separated, int step) {

  double result = 0;
  result += plaket_plane5(separated[3], separated[0], 1, t_size, t_size,
                          x_size * t_size, 1);
  result += plaket_plane5(separated[3], separated[1], 1, t_size,
                          x_size * t_size, y_size * x_size * t_size, 1);
  result += plaket_plane5(separated[3], separated[2], 1, t_size,
                          y_size * x_size * t_size,
                          y_size * x_size * t_size * z_size, 1);
  result += plaket_plane5(separated[0], separated[1], t_size, x_size * t_size,
                          x_size * t_size, y_size * x_size * t_size, step);
  result += plaket_plane5(separated[0], separated[2], t_size, x_size * t_size,
                          y_size * x_size * t_size,
                          y_size * x_size * t_size * z_size, step);
  result += plaket_plane5(separated[1], separated[2], x_size * t_size,
                          y_size * x_size * t_size, y_size * x_size * t_size,
                          y_size * x_size * t_size * z_size, step);
  return result / 6;
}

template <class T>
double plaket_plane6(std::vector<T> &conf_mu, std::vector<T> &conf_nu,
                     int size_mu1, int size_mu2, int size_nu1, int size_nu2,
                     int thread_num) {
  int data_size = x_size * y_size * z_size * t_size;

  T plakets;
  double result = 0;

#pragma omp parallel for collapse(3) private(plakets) reduction(+ : result) num_threads(thread_num)
  for (int k = 0; k < data_size; k += size_nu2) {
    for (int i = 0; i < size_nu2; i += size_mu2) {
      for (int j = 0; j < size_mu2; j++) {
        if (j < size_mu2 - size_mu1)
          plakets = conf_mu[i + k + j] * conf_nu[i + k + j + size_mu1];
        else
          plakets =
              conf_mu[i + k + j] * conf_nu[i + k + j - size_mu2 + size_mu1];
        if (i + j < size_nu2 - size_nu1)
          plakets = plakets ^ &conf_mu[i + k + j + size_nu1];
        else
          plakets = plakets ^ &conf_mu[i + k + j - size_nu2 + size_nu1];
        result += plakets.multiply_tr(&conf_nu[i + k + j]);
      }
    }
  }

  return result / data_size;
}

template <class T>
double plaket_time_test10(std::vector<std::vector<T>> &separated,
                          int thread_num) {

  double start_time;
  double end_time;
  double search_time;

  start_time = omp_get_wtime();

  double result = 0;
  result += plaket_plane6(separated[3], separated[0], 1, t_size, t_size,
                          x_size * t_size, thread_num);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "plaket_time_test10 time 1: " << search_time << std::endl;

  start_time = omp_get_wtime();

  result +=
      plaket_plane6(separated[3], separated[1], 1, t_size, x_size * t_size,
                    y_size * x_size * t_size, thread_num);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "plaket_time_test10 time 2: " << search_time << std::endl;

  start_time = omp_get_wtime();

  result += plaket_plane6(separated[3], separated[2], 1, t_size,
                          y_size * x_size * t_size,
                          y_size * x_size * t_size * z_size, thread_num);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "plaket_time_test10 time 3: " << search_time << std::endl;

  start_time = omp_get_wtime();

  result +=
      plaket_plane6(separated[0], separated[1], t_size, x_size * t_size,
                    x_size * t_size, y_size * x_size * t_size, thread_num);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "plaket_time_test10 time 4: " << search_time << std::endl;

  start_time = omp_get_wtime();

  result += plaket_plane6(separated[0], separated[2], t_size, x_size * t_size,
                          y_size * x_size * t_size,
                          y_size * x_size * t_size * z_size, thread_num);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "plaket_time_test10 time 5: " << search_time << std::endl;

  start_time = omp_get_wtime();

  result += plaket_plane6(separated[1], separated[2], x_size * t_size,
                          y_size * x_size * t_size, y_size * x_size * t_size,
                          y_size * x_size * t_size * z_size, thread_num);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "plaket_time_test10 time 6: " << search_time << std::endl;

  return result / 6;
}

template <class T>
std::vector<std::vector<T>> separate_3(std::vector<T> &conf) {
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
template double plaket_time_test2(std::vector<std::vector<su2>> &separated);

template std::vector<std::vector<su2>> separate_3(std::vector<su2> &conf);
template double plaket_plane1(std::vector<su2> &conf_mu,
                              std::vector<su2> &conf_nu, int size1, int size2);
template double plaket_time_test3(std::vector<std::vector<su2>> &separated);

template double plaket_plane2(std::vector<su2> &conf_mu,
                              std::vector<su2> &conf_nu, int size1, int size2);
template double plaket_time_test4(std::vector<std::vector<su2>> &separated);
template double plaket_time_test5(std::vector<su2> &conf);
template double plaket_time_test6(std::vector<su2> &conf);
template double plaket_plane3(std::vector<su2> &conf, int mu, int nu,
                              int size_mu1, int size_mu2, int size_nu1,
                              int size_nu2);

template double plaket_time_test7(std::vector<su2> &conf);

template double plaket_plane4(std::vector<su2> &conf_mu,
                              std::vector<su2> &conf_nu, int size_mu1,
                              int size_mu2, int size_nu1, int size_nu2);
template double plaket_time_test8(std::vector<std::vector<su2>> &separated);

template double plaket_plane5(std::vector<su2> &conf_mu,
                              std::vector<su2> &conf_nu, int size_mu1,
                              int size_mu2, int size_nu1, int size_nu2,
                              int step);
template double plaket_time_test9(std::vector<std::vector<su2>> &separated,
                                  int step);

template double plaket_plane6(std::vector<su2> &conf_mu,
                              std::vector<su2> &conf_nu, int size_mu1,
                              int size_mu2, int size_nu1, int size_nu2,
                              int thread_num);
template double plaket_time_test10(std::vector<std::vector<su2>> &separated,
                                   int step);

// abelian

template double plaket_plane(std::vector<abelian> &conf_mu,
                             std::vector<abelian> &conf_nu, int size1,
                             int size2);

template double plaket_time_test(std::vector<std::vector<abelian>> &separated);
template double plaket_time_test1(std::vector<std::vector<abelian>> &separated);
template double plaket_time_test2(std::vector<std::vector<abelian>> &separated);

template double plaket_plane1(std::vector<abelian> &conf_mu,
                              std::vector<abelian> &conf_nu, int size1,
                              int size2);
template double plaket_time_test3(std::vector<std::vector<abelian>> &separated);

template double plaket_plane2(std::vector<abelian> &conf_mu,
                              std::vector<abelian> &conf_nu, int size1,
                              int size2);
template double plaket_time_test4(std::vector<std::vector<abelian>> &separated);

template double plaket_time_test5(std::vector<abelian> &conf);
template double plaket_time_test6(std::vector<abelian> &conf);

template double plaket_plane3(std::vector<abelian> &conf, int mu, int nu,
                              int size_mu1, int size_mu2, int size_nu1,
                              int size_nu2);

template double plaket_time_test7(std::vector<abelian> &conf);

template double plaket_plane4(std::vector<abelian> &conf_mu,
                              std::vector<abelian> &conf_nu, int size_mu1,
                              int size_mu2, int size_nu1, int size_nu2);
template double plaket_time_test8(std::vector<std::vector<abelian>> &separated);

template double plaket_plane5(std::vector<abelian> &conf_mu,
                              std::vector<abelian> &conf_nu, int size_mu1,
                              int size_mu2, int size_nu1, int size_nu2,
                              int step);
template double plaket_time_test9(std::vector<std::vector<abelian>> &separated,
                                  int step);

template double plaket_plane6(std::vector<abelian> &conf_mu,
                              std::vector<abelian> &conf_nu, int size_mu1,
                              int size_mu2, int size_nu1, int size_nu2,
                              int thread_num);
template double plaket_time_test10(std::vector<std::vector<abelian>> &separated,
                                   int step);

template std::vector<std::vector<abelian>>
separate_3(std::vector<abelian> &conf);

// su3_full
template double plaket_plane(std::vector<su3_full> &conf_mu,
                             std::vector<su3_full> &conf_nu, int size1,
                             int size2);

template double plaket_time_test(std::vector<std::vector<su3_full>> &separated);
template double
plaket_time_test1(std::vector<std::vector<su3_full>> &separated);
template double
plaket_time_test2(std::vector<std::vector<su3_full>> &separated);

template double plaket_plane1(std::vector<su3_full> &conf_mu,
                              std::vector<su3_full> &conf_nu, int size1,
                              int size2);
template double
plaket_time_test3(std::vector<std::vector<su3_full>> &separated);

template std::vector<std::vector<su3_full>>
separate_3(std::vector<su3_full> &conf);

template double plaket_plane2(std::vector<su3_full> &conf_mu,
                              std::vector<su3_full> &conf_nu, int size1,
                              int size2);
template double
plaket_time_test4(std::vector<std::vector<su3_full>> &separated);

template double plaket_time_test5(std::vector<su3_full> &conf);
template double plaket_time_test6(std::vector<su3_full> &conf);

template double plaket_plane3(std::vector<su3_full> &conf, int mu, int nu,
                              int size_mu1, int size_mu2, int size_nu1,
                              int size_nu2);

template double plaket_time_test7(std::vector<su3_full> &conf);

template double plaket_plane4(std::vector<su3_full> &conf_mu,
                              std::vector<su3_full> &conf_nu, int size_mu1,
                              int size_mu2, int size_nu1, int size_nu2);
template double
plaket_time_test8(std::vector<std::vector<su3_full>> &separated);

template double plaket_plane5(std::vector<su3_full> &conf_mu,
                              std::vector<su3_full> &conf_nu, int size_mu1,
                              int size_mu2, int size_nu1, int size_nu2,
                              int step);
template double plaket_time_test9(std::vector<std::vector<su3_full>> &separated,
                                  int step);

template double plaket_plane6(std::vector<su3_full> &conf_mu,
                              std::vector<su3_full> &conf_nu, int size_mu1,
                              int size_mu2, int size_nu1, int size_nu2,
                              int thread_num);
template double
plaket_time_test10(std::vector<std::vector<su3_full>> &separated, int step);

// su3
template double plaket_plane(std::vector<su3> &conf_mu,
                             std::vector<su3> &conf_nu, int size1, int size2);

template double plaket_time_test(std::vector<std::vector<su3>> &separated);
template double plaket_time_test1(std::vector<std::vector<su3>> &separated);
template double plaket_time_test2(std::vector<std::vector<su3>> &separated);

template double plaket_plane1(std::vector<su3> &conf_mu,
                              std::vector<su3> &conf_nu, int size1, int size2);
template double plaket_time_test3(std::vector<std::vector<su3>> &separated);

template std::vector<std::vector<su3>> separate_3(std::vector<su3> &conf);

template double plaket_plane2(std::vector<su3> &conf_mu,
                              std::vector<su3> &conf_nu, int size1, int size2);
template double plaket_time_test4(std::vector<std::vector<su3>> &separated);

template double plaket_time_test5(std::vector<su3> &conf);
template double plaket_time_test6(std::vector<su3> &conf);

template double plaket_plane3(std::vector<su3> &conf, int mu, int nu,
                              int size_mu1, int size_mu2, int size_nu1,
                              int size_nu2);

template double plaket_time_test7(std::vector<su3> &conf);

template double plaket_plane4(std::vector<su3> &conf_mu,
                              std::vector<su3> &conf_nu, int size_mu1,
                              int size_mu2, int size_nu1, int size_nu2);
template double plaket_time_test8(std::vector<std::vector<su3>> &separated);

template double plaket_plane5(std::vector<su3> &conf_mu,
                              std::vector<su3> &conf_nu, int size_mu1,
                              int size_mu2, int size_nu1, int size_nu2,
                              int step);
template double plaket_time_test9(std::vector<std::vector<su3>> &separated,
                                  int step);

template double plaket_plane6(std::vector<su3> &conf_mu,
                              std::vector<su3> &conf_nu, int size_mu1,
                              int size_mu2, int size_nu1, int size_nu2,
                              int thread_num);
template double plaket_time_test10(std::vector<std::vector<su3>> &separated,
                                   int step);
