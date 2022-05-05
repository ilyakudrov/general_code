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
std::vector<std::vector<T>> separate_wilson(std::vector<T> &conf) {
  int data_size = x_size * y_size * z_size * t_size;
  std::vector<std::vector<T>> result(4, std::vector<T>(data_size));

  link1 link(x_size, y_size, z_size, t_size);

  vector<int> multiplier_new(4);
  int place_new;

  for (int mu = 0; mu < 4; ++mu) {

    multiplier_new[mu] = 1;
    for (int nu = 1; nu < 4; ++nu) {
      multiplier_new[(nu + mu) % 4] = multiplier_new[(nu + mu - 1) % 4] *
                                      link.lattice_size[(nu + mu - 1) % 4];
    }

    SPACE_ITER_START

    place_new = 0;
    for (int i = 0; i < 4; i++) {
      place_new += multiplier_new[i] * link.coordinate[i];
    }

    result[mu][place_new] = conf[link.place + mu];

    SPACE_ITER_END
  }

  return result;
}

template <class T>
std::vector<T> wilson_lines_test1(std::vector<T> separated, int length,
                                  int step) {
  int data_size = x_size * y_size * z_size * t_size;

  std::vector<T> wilson_lines(data_size);
  T A;

  for (int i = 0; i < data_size; i += step) {
    for (int j = i; j < i + step - length + 1; j++) {
      A = T();
      for (int k = j; k < j + length; k++) {
        A = A * separated[k];
      }
      wilson_lines[j] = A;
    }

    for (int j = i + step - length + 1; j < i + step; j++) {
      A = T();
      for (int k = j; k < i + step; k++) {
        A = A * separated[k];
      }
      for (int k = i; k < i + (length - (i + step - j)); k++) {
        A = A * separated[k];
      }
      wilson_lines[j] = A;
    }
  }

  return wilson_lines;
}

template <class T>
std::vector<T> wilson_lines_test2(std::vector<T> separated, int length,
                                  int step) {
  int data_size = x_size * y_size * z_size * t_size;

  std::vector<T> wilson_lines(data_size);
  T A;

  for (int i = 0; i < data_size; i += step) {
    A = T();
    for (int j = i; j < i + length; j++) {
      A = A * separated[j];
    }
    wilson_lines[i] = A;
    for (int j = i + 1; j < i + step - length + 1; j++) {
      A = separated[j - 1] % A;
      A = A * separated[j + length - 1];
      wilson_lines[j] = A;
    }

    for (int j = i + step - length + 1; j < i + step; j++) {
      A = separated[j - 1] % A;
      A = A * separated[j - step + length - 1];
      wilson_lines[j] = A;
    }
  }

  return wilson_lines;
}

template <class T>
std::vector<T> wilson_lines_increase_test(std::vector<T> separated,
                                          std::vector<T> wilson_lines_old,
                                          int length, int step) {
  int data_size = x_size * y_size * z_size * t_size;

  std::vector<T> wilson_lines(data_size);

  for (int i = 0; i < data_size; i += step) {
    for (int j = i; j < i + step - length + 1; j++) {
      wilson_lines[j] = wilson_lines_old[j] * separated[j + length - 1];
    }
    for (int j = i + step - length + 1; j < i + step; j++) {
      wilson_lines[j] = wilson_lines_old[j] * separated[j - step + length - 1];
    }
  }

  return wilson_lines;
}

template <class T>
std::vector<std::vector<T>> separate_wilson_unchanged(std::vector<T> &conf) {
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

// with unchanged separated conf
template <class T>
std::vector<T> wilson_lines_test3(std::vector<T> separated, int length,
                                  int size1, int size2) {
  int data_size = x_size * y_size * z_size * t_size;

  std::vector<T> wilson_lines(data_size);
  T A;

  for (int i = 0; i < data_size; i += size2) {
    for (int k = i; k < i + size1; k++) {
      A = T();
      for (int j = k; j < k + length * size1; j += size1) {
        A = A * separated[j];
      }
      wilson_lines[k] = A;
      for (int j = k + size1; j < k + size2 - (length - 1) * size1;
           j += size1) {
        A = separated[j - size1] % A;
        A = A * separated[j + (length - 1) * size1];
        wilson_lines[j] = A;
      }

      for (int j = k + size2 - (length - 1) * size1; j < k + size2;
           j += size1) {
        A = separated[j - size1] % A;
        A = A * separated[j - size2 + (length - 1) * size1];
        wilson_lines[j] = A;
      }
    }
  }

  return wilson_lines;
}

template <class T>
double wilson_plane_test1(std::vector<T> &wilson_lines_mu,
                          std::vector<T> &wilson_lines_nu, int size_mu1,
                          int size_mu2, int size_nu1, int size_nu2,
                          int length_mu, int length_nu, int thread_num) {
  int data_size = x_size * y_size * z_size * t_size;

  T plakets;
  double result = 0;

#pragma omp parallel for collapse(3) private(plakets) reduction(+ : result) num_threads(thread_num)
  for (int k = 0; k < data_size; k += size_nu2) {
    for (int i = 0; i < size_nu2; i += size_mu2) {
      for (int j = 0; j < size_mu2; j++) {
        if (j < size_mu2 - length_mu * size_mu1)
          plakets = wilson_lines_mu[i + k + j] *
                    wilson_lines_nu[i + k + j + length_mu * size_mu1];
        else
          plakets =
              wilson_lines_mu[i + k + j] *
              wilson_lines_nu[i + k + j - size_mu2 + length_mu * size_mu1];
        if (i + j < size_nu2 - length_nu * size_nu1)
          plakets = plakets ^ wilson_lines_mu[i + k + j + length_nu * size_nu1];
        else
          plakets =
              plakets ^
              wilson_lines_mu[i + k + j - size_nu2 + length_nu * size_nu1];
        result += plakets.multiply_tr(wilson_lines_nu[i + k + j]);
      }
    }
  }

  return result / data_size;
}

template <class T>
double wilson_loop_test_time(std::vector<std::vector<T>> &wilson_lines,
                             int length_R, int length_T, int thread_num) {

  double start_time;
  double end_time;
  double search_time;

  double total_time = 0;

  start_time = omp_get_wtime();

  double result = 0;
  result += wilson_plane_test1(
      wilson_lines[0], wilson_lines[3], 1, x_size, x_size * y_size * z_size,
      x_size * y_size * z_size * t_size, length_R, length_T, thread_num);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  total_time += search_time;
  std::cout << "wilson_loop_time time 0: " << search_time << std::endl;

  start_time = omp_get_wtime();

  result += wilson_plane_test1(wilson_lines[1], wilson_lines[3], x_size,
                               x_size * y_size, x_size * y_size * z_size,
                               x_size * y_size * z_size * t_size, length_R,
                               length_T, thread_num);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  total_time += search_time;
  std::cout << "wilson_loop_time time 1: " << search_time << std::endl;

  start_time = omp_get_wtime();

  result += wilson_plane_test1(
      wilson_lines[2], wilson_lines[3], x_size * y_size,
      x_size * y_size * z_size, x_size * y_size * z_size,
      x_size * y_size * z_size * t_size, length_R, length_T, thread_num);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  total_time += search_time;
  std::cout << "wilson_loop_time time 2: " << search_time << std::endl;
  std::cout << "wilson_loop_time total time: " << total_time << std::endl;

  return result / 3;
}

// abelian
template std::vector<std::vector<abelian>>
separate_wilson(std::vector<abelian> &conf);

template std::vector<abelian> wilson_lines_test1(std::vector<abelian> separated,
                                                 int length, int step);

template std::vector<abelian> wilson_lines_test2(std::vector<abelian> separated,
                                                 int length, int step);

template std::vector<abelian>
wilson_lines_increase_test(std::vector<abelian> separaed,
                           std::vector<abelian> wilson_lines_old, int length,
                           int step);

template std::vector<std::vector<abelian>>
separate_wilson_unchanged(std::vector<abelian> &conf);

template std::vector<abelian> wilson_lines_test3(std::vector<abelian> separated,
                                                 int length, int size1,
                                                 int size2);

template double wilson_plane_test1(std::vector<abelian> &wilson_lines_mu,
                                   std::vector<abelian> &wilson_lines_nu,
                                   int size_mu1, int size_mu2, int size_nu1,
                                   int size_nu2, int length_mu, int length_nu,
                                   int thread_num);

template double
wilson_loop_test_time(std::vector<std::vector<abelian>> &wilson_lines,
                      int length_R, int length_T, int thread_num);

// su2
template std::vector<std::vector<su2>> separate_wilson(std::vector<su2> &conf);

template std::vector<su2> wilson_lines_test1(std::vector<su2> separated,
                                             int length, int step);

template std::vector<su2> wilson_lines_test2(std::vector<su2> separated,
                                             int length, int step);

template std::vector<su2>
wilson_lines_increase_test(std::vector<su2> separated,
                           std::vector<su2> wilson_lines_old, int length,
                           int step);

template std::vector<std::vector<su2>>
separate_wilson_unchanged(std::vector<su2> &conf);

template std::vector<su2> wilson_lines_test3(std::vector<su2> separated,
                                             int length, int size1, int size2);

template double wilson_plane_test1(std::vector<su2> &wilson_lines_mu,
                                   std::vector<su2> &wilson_lines_nu,
                                   int size_mu1, int size_mu2, int size_nu1,
                                   int size_nu2, int length_mu, int length_nu,
                                   int thread_num);

template double
wilson_loop_test_time(std::vector<std::vector<su2>> &wilson_lines, int length_R,
                      int length_T, int thread_num);

// su3
template std::vector<std::vector<su3>> separate_wilson(std::vector<su3> &conf);

template std::vector<su3> wilson_lines_test1(std::vector<su3> separated,
                                             int length, int step);

template std::vector<su3> wilson_lines_test2(std::vector<su3> separated,
                                             int length, int step);

template std::vector<su3>
wilson_lines_increase_test(std::vector<su3> separated,
                           std::vector<su3> wilson_lines_old, int length,
                           int step);

template std::vector<std::vector<su3>>
separate_wilson_unchanged(std::vector<su3> &conf);

template std::vector<su3> wilson_lines_test3(std::vector<su3> separated,
                                             int length, int size1, int size2);

template double wilson_plane_test1(std::vector<su3> &wilson_lines_mu,
                                   std::vector<su3> &wilson_lines_nu,
                                   int size_mu1, int size_mu2, int size_nu1,
                                   int size_nu2, int length_mu, int length_nu,
                                   int thread_num);

template double
wilson_loop_test_time(std::vector<std::vector<su3>> &wilson_lines, int length_R,
                      int length_T, int thread_num);
