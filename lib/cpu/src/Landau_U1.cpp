#include "../include/Landau_U1.h"
#include "../include/data.h"
#include "../include/indexing.h"
#include "../include/link.h"
#include "../include/matrix.h"

#include <cmath>
#include <complex>
#include <map>
#include <math.h>
#include <omp.h>
#include <random>
#include <tuple>
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

std::vector<double> convert_to_angles(const std::vector<su2> &conf_su2) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  std::vector<double> conf_abelian;
  conf_abelian.reserve(data_size);

  for (int i = 0; i < data_size; i++) {
    conf_abelian.push_back(atan2(conf_su2[i].a3, conf_su2[i].a0));
  }

  return conf_abelian;
}

std::vector<double> convert_to_angles(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf_su2) {
  DataPatternLexicographical data_pattern(conf_su2.lat_dim);
  int data_size = data_pattern.get_data_size();
  std::vector<double> conf_angles(data_size);
  for (int i = 0; i < data_size; i++) {
    conf_angles[i] = atan2(conf_su2[i].a3, conf_su2[i].a0);
  }
  return conf_angles;
}

std::vector<abelian> convert_to_abelian(const std::vector<su2> &conf_su2) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  std::vector<abelian> conf_abelian;
  conf_abelian.reserve(data_size);

  for (int i = 0; i < data_size; i++) {
    conf_abelian.push_back(abelian(1, atan2(conf_su2[i].a3, conf_su2[i].a0)));
  }

  return conf_abelian;
}

std::vector<std::complex<double>>
convert_to_complex(const std::vector<su2> &conf_su2) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  std::vector<std::complex<double>> conf_complex;
  conf_complex.reserve(data_size);

  double module;

  for (int i = 0; i < data_size; i++) {
    module =
        sqrt(conf_su2[i].a0 * conf_su2[i].a0 + conf_su2[i].a3 * conf_su2[i].a3);
    conf_complex.push_back(
        std::complex<double>(conf_su2[i].a0 / module, conf_su2[i].a3 / module));
  }

  return conf_complex;
}

std::vector<std::complex<double>> convert_to_complex(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf_su2) {
  DataPatternLexicographical data_pattern(conf_su2.lat_dim);
  int data_size = data_pattern.get_data_size();
  std::vector<std::complex<double>> conf_complex(data_size);
  double module;
  for (int i = 0; i < data_size; i++) {
    module =
        sqrt(conf_su2[i].a0 * conf_su2[i].a0 + conf_su2[i].a3 * conf_su2[i].a3);
    conf_complex[i] =
        std::complex<double>(conf_su2[i].a0 / module, conf_su2[i].a3 / module);
  }
  return conf_complex;
}

std::vector<std::complex<double>> convert_to_complex(
    const Data::LatticeData<DataPatternLexicographical, abelian> &conf) {
  DataPatternLexicographical data_pattern(conf.lat_dim);
  int data_size = data_pattern.get_data_size();
  std::vector<std::complex<double>> conf_complex(data_size);
  for (int i = 0; i < data_size; i++) {
    conf_complex[i] = std::complex<double>(cos(conf[i].phi), sin(conf[i].phi));
  }
  return conf_complex;
}

std::vector<double> convert_complex_to_angles(
    const std::vector<std::complex<double>> &conf_complex) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  std::vector<double> conf_angles;
  conf_angles.reserve(data_size);

  double module;

  for (int i = 0; i < data_size; i++) {
    conf_angles.push_back(
        atan2(conf_complex[i].imag(), conf_complex[i].real()));
  }

  return conf_angles;
}

std::vector<abelian> convert_complex_to_abelian(
    const std::vector<std::complex<double>> &conf_complex,
    DataPatternLexicographical &data_pattern) {
  int data_size = data_pattern.get_data_size();
  std::vector<abelian> conf_angles(data_size);
  double module;
  for (int i = 0; i < data_size; i++) {
    conf_angles[i] =
        abelian(1, atan2(conf_complex[i].imag(), conf_complex[i].real()));
  }
  return conf_angles;
}

std::vector<double> generate_gauge_angles_uniform() {
  unsigned seed = time(NULL);
  // unsigned seed = 123;
  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);
  int data_size = x_size * y_size * z_size * t_size;
  std::vector<double> gauge_abelian;
  gauge_abelian.reserve(data_size);
  for (int i = 0; i < data_size; i++) {
    gauge_abelian.push_back(
        (2 * (double)random_generator() / random_generator.max() - 1) * M_PI);
  }
  return gauge_abelian;
}

std::vector<double>
generate_gauge_angles_uniform(DataPatternLexicographical &data_pattern) {
  unsigned seed = time(NULL);
  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);
  int data_size = data_pattern.get_lattice_size();
  std::vector<double> gauge_abelian(data_size);
  for (int i = 0; i < data_size; i++) {
    gauge_abelian[i] =
        (2 * (double)random_generator() / (random_generator.max() + 1) - 1) *
        M_PI;
  }
  return gauge_abelian;
}

std::vector<abelian> generate_gauge_abelian_uniform() {

  unsigned seed = time(NULL);
  // unsigned seed = 123;

  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);

  int data_size = x_size * y_size * z_size * t_size;

  std::vector<abelian> gauge_abelian;
  gauge_abelian.reserve(data_size);

  for (int i = 0; i < data_size; i++) {

    gauge_abelian.push_back(abelian(
        1,
        (2 * (double)random_generator() / random_generator.max() - 1) * M_PI));
  }

  return gauge_abelian;
}

std::vector<std::complex<double>> generate_gauge_complex_uniform() {

  unsigned seed = time(NULL);
  // unsigned seed = 123;

  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);

  int data_size = x_size * y_size * z_size * t_size;

  std::vector<std::complex<double>> gauge_complex;
  gauge_complex.reserve(data_size);

  double angle_tmp;

  for (int i = 0; i < data_size; i++) {

    angle_tmp =
        (2 * (double)random_generator() / random_generator.max() - 1) * M_PI;

    gauge_complex.push_back(
        std::complex<double>(cos(angle_tmp), sin(angle_tmp)));
  }

  return gauge_complex;
}

std::vector<std::complex<double>>
generate_gauge_complex_uniform(DataPatternLexicographical &data_pattern) {
  unsigned seed = time(NULL);
  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);
  int data_size = data_pattern.get_lattice_size();
  std::vector<std::complex<double>> gauge_complex(data_size);
  double angle_tmp;
  for (int i = 0; i < data_size; i++) {
    angle_tmp =
        (2 * (double)random_generator() / (random_generator.max() + 1) - 1) *
        M_PI;
    gauge_complex[i] = std::complex<double>(cos(angle_tmp), sin(angle_tmp));
  }
  return gauge_complex;
}

std::vector<std::complex<double>> generate_gauge_complex_unity() {

  int data_size = x_size * y_size * z_size * t_size;

  std::vector<std::complex<double>> gauge_complex;
  gauge_complex.reserve(data_size);

  for (int i = 0; i < data_size; i++) {
    gauge_complex.push_back(std::complex<double>(1, 0));
  }

  return gauge_complex;
}

std::vector<abelian>
generate_gauge_abelian_uniform(DataPatternLexicographical &data_pattern) {
  unsigned seed = time(NULL);
  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);
  int data_size = data_pattern.get_lattice_size();
  std::vector<abelian> gauge_abelian;
  gauge_abelian.reserve(data_size);
  for (int i = 0; i < data_size; i++) {
    gauge_abelian.push_back(abelian(
        1,
        (2 * (double)random_generator() / random_generator.max() - 1) * M_PI));
  }
  return gauge_abelian;
}

std::vector<std::complex<double>>
generate_gauge_complex_unity(DataPatternLexicographical &data_pattern) {
  int data_size = data_pattern.get_lattice_size();
  std::vector<std::complex<double>> gauge_complex;
  gauge_complex.reserve(data_size);
  for (int i = 0; i < data_size; i++) {
    gauge_complex.push_back(std::complex<double>(1, 0));
  }
  return gauge_complex;
}

std::vector<double> generate_random_numbers1(int vector_size) {
  unsigned seed = time(NULL);
  // unsigned seed = 123;

  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);

  std::vector<double> random_numbers;
  random_numbers.reserve(vector_size);

  double a;

  for (int i = 0; i < vector_size; i++) {
    do {
      a = (double)random_generator() / random_generator.max();

      // to aviod errors due to rounding
    } while (a == 0 || a == 1);

    random_numbers.push_back(a);
  }

  return random_numbers;
}

void heat_bath_fast(
    std::complex<double> &gauge_complex, const std::complex<double> &neighbour,
    double &temperature,
    std::subtract_with_carry_engine<unsigned, 24, 10, 24> &random_generator) {
  double temp_factor;
  // generation of new std::complex<double> variable
  temp_factor = std::norm(neighbour) / temperature;
  // generation of random variable with distribution p ~ exp(ax)
  double x;
  do {
    x = 1 + log((double)random_generator() / (random_generator.max() + 1)) /
                temp_factor;
  } while (x < -1);
  double alpha;
  if (random_generator() % 2 == 0) {
    alpha = std::acos(x) - std::atan2(neighbour.imag(), neighbour.real());
  } else {
    alpha = -std::acos(x) - std::atan2(neighbour.imag(), neighbour.real());
  }
  gauge_complex = std::complex<double>(cos(alpha), sin(alpha));
}

void heat_bath(
    std::complex<double> &gauge_complex, const std::complex<double> &neighbour,
    double &temperature,
    std::subtract_with_carry_engine<unsigned, 24, 10, 24> &random_generator) {
  double temp_factor;
  long double temp_exponent;
  // generation of new std::complex<double> variable
  temp_factor = std::norm(neighbour) / temperature;
  temp_exponent = exp(temp_factor);
  // generation of random variable with distribution p ~ exp(ax)
  double x;
  x = log(((double)(random_generator.min() + 1) / random_generator.max() +
           (double)random_generator() / (random_generator.max() + 2)) *
              (temp_exponent - 1 / temp_exponent) +
          1 / temp_exponent) /
      temp_factor;
  double alpha;
  if (random_generator() % 2 == 0) {
    alpha = std::acos(x) - std::atan2(neighbour.imag(), neighbour.real());
  } else {
    alpha = -std::acos(x) - std::atan2(neighbour.imag(), neighbour.real());
  }
  gauge_complex = std::complex<double>(cos(alpha), sin(alpha));
}

// std::complex<double>
// contribution_site(std::vector<std::complex<double>> &gauge_complex,
//                   const std::vector<std::complex<double>> &conf_complex,
//                   DataPatternLexicographical &data_pattern, int x, int y, int
//                   z, int t, int position, std::vector<int> &shift) {
//   std::complex<double> contribution = 0;
//   // mu = 0
//   if (x < data_pattern.lat_dim[0] - 1) {
//     contribution += conf_complex[position * 4] *
//                     std::conj(gauge_complex[position + shift[0]]);
//   } else
//     contribution += conf_complex[position * 4] *
//                     std::conj(gauge_complex[position + shift[0] - shift[1]]);
//   if (x > 0)
//     contribution += std::conj(conf_complex[(position - shift[0]) * 4]) *
//                     std::conj(gauge_complex[position - shift[0]]);
//   else
//     contribution +=
//         std::conj(conf_complex[(position - shift[0] + shift[1]) * 4]) *
//         std::conj(gauge_complex[position - shift[0] + shift[1]]);

//   // mu = 1
//   if (y < data_pattern.lat_dim[1] - 1)
//     contribution += conf_complex[position * 4 + 1] *
//                     std::conj(gauge_complex[position + shift[1]]);
//   else
//     contribution += conf_complex[position * 4 + 1] *
//                     std::conj(gauge_complex[position + shift[1] - shift[2]]);
//   if (y > 0)
//     contribution += std::conj(conf_complex[(position - shift[1]) * 4 + 1]) *
//                     std::conj(gauge_complex[position - shift[1]]);
//   else
//     contribution +=
//         std::conj(conf_complex[(position - shift[1] + shift[2]) * 4 + 1]) *
//         std::conj(gauge_complex[position - shift[1] + shift[2]]);

//   // mu = 2
//   if (z < data_pattern.lat_dim[2] - 1)
//     contribution += conf_complex[position * 4 + 2] *
//                     std::conj(gauge_complex[position + shift[2]]);
//   else
//     contribution += conf_complex[position * 4 + 2] *
//                     std::conj(gauge_complex[position + shift[2] - shift[3]]);
//   if (z > 0)
//     contribution += std::conj(conf_complex[(position - shift[2]) * 4 + 2]) *
//                     std::conj(gauge_complex[position - shift[2]]);
//   else
//     contribution +=
//         std::conj(conf_complex[(position - shift[2] + shift[3]) * 4 + 2]) *
//         std::conj(gauge_complex[position - shift[2] + shift[3]]);

//   // mu = 3
//   if (t < data_pattern.lat_dim[3] - 1)
//     contribution += conf_complex[position * 4 + 3] *
//                     std::conj(gauge_complex[position + shift[3]]);
//   else
//     contribution += conf_complex[position * 4 + 3] *
//                     std::conj(gauge_complex[position + shift[3] - shift[4]]);
//   if (t > 0)
//     contribution += std::conj(conf_complex[(position - shift[3]) * 4 + 3]) *
//                     std::conj(gauge_complex[position - shift[3]]);
//   else
//     contribution +=
//         std::conj(conf_complex[(position - shift[3] + shift[4]) * 4 + 3]) *
//         std::conj(gauge_complex[position - shift[3] + shift[4]]);
//   return contribution;
// }

inline void contribution_conj(std::complex<double> &contribution,
                              const std::complex<double> &a,
                              const std::complex<double> &b) {
  contribution +=
      std::complex<double>(a.real() * b.real() + a.imag() * b.imag(),
                           a.imag() * b.real() - a.real() * b.imag());
}

inline void contribution_conj_conj(std::complex<double> &contribution,
                                   const std::complex<double> &a,
                                   const std::complex<double> &b) {
  contribution +=
      std::complex<double>(a.real() * b.real() - a.imag() * b.imag(),
                           -a.imag() * b.real() - a.real() * b.imag());
}

std::complex<double>
contribution_site(std::vector<std::complex<double>> &gauge_complex,
                  const std::vector<std::complex<double>> &conf_complex,
                  DataPatternLexicographical &data_pattern, int x, int y, int z,
                  int t, int position, std::vector<int> &shift) {
  std::complex<double> contribution = 0;
  // mu = 0
  if (x < data_pattern.lat_dim[0] - 1) {
    contribution_conj(contribution, conf_complex[position * 4],
                      gauge_complex[position + shift[0]]);
  } else
    contribution_conj(contribution, conf_complex[position * 4],
                      gauge_complex[position + shift[0] - shift[1]]);
  if (x > 0)
    contribution_conj_conj(contribution,
                           conf_complex[(position - shift[0]) * 4],
                           gauge_complex[position - shift[0]]);
  else
    contribution_conj_conj(contribution,
                           conf_complex[(position - shift[0] + shift[1]) * 4],
                           gauge_complex[position - shift[0] + shift[1]]);
  // mu = 1
  if (y < data_pattern.lat_dim[1] - 1)
    contribution_conj(contribution, conf_complex[position * 4 + 1],
                      gauge_complex[position + shift[1]]);
  else
    contribution_conj(contribution, conf_complex[position * 4 + 1],
                      gauge_complex[position + shift[1] - shift[2]]);
  if (y > 0)
    contribution_conj_conj(contribution,
                           conf_complex[(position - shift[1]) * 4 + 1],
                           gauge_complex[position - shift[1]]);
  else
    contribution_conj_conj(
        contribution, conf_complex[(position - shift[1] + shift[2]) * 4 + 1],
        gauge_complex[position - shift[1] + shift[2]]);
  // mu = 2
  if (z < data_pattern.lat_dim[2] - 1)
    contribution_conj(contribution, conf_complex[position * 4 + 2],
                      gauge_complex[position + shift[2]]);
  else
    contribution_conj(contribution, conf_complex[position * 4 + 2],
                      gauge_complex[position + shift[2] - shift[3]]);
  if (z > 0)
    contribution_conj_conj(contribution,
                           conf_complex[(position - shift[2]) * 4 + 2],
                           gauge_complex[position - shift[2]]);
  else
    contribution_conj_conj(
        contribution, conf_complex[(position - shift[2] + shift[3]) * 4 + 2],
        gauge_complex[position - shift[2] + shift[3]]);
  // mu = 3
  if (t < data_pattern.lat_dim[3] - 1)
    contribution_conj(contribution, conf_complex[position * 4 + 3],
                      gauge_complex[position + shift[3]]);
  else
    contribution_conj(contribution, conf_complex[position * 4 + 3],
                      gauge_complex[position + shift[3] - shift[4]]);
  if (t > 0)
    contribution_conj_conj(contribution,
                           conf_complex[(position - shift[3]) * 4 + 3],
                           gauge_complex[position - shift[3]]);
  else
    contribution_conj_conj(
        contribution, conf_complex[(position - shift[3] + shift[4]) * 4 + 3],
        gauge_complex[position - shift[3] + shift[4]]);
  return contribution;
}

void heat_bath_update(std::vector<std::complex<double>> &gauge_angles,
                      const std::vector<std::complex<double>> &conf_angles,
                      DataPatternLexicographical &data_pattern,
                      double temperature) {
  std::vector<int> shift = {1, data_pattern.lat_dim[0],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2] *
                                data_pattern.lat_dim[3]};
  unsigned seed = time(NULL);
  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);
  std::complex<double> contribution;
  int position;
#pragma omp parallel for collapse(4) private(contribution, position)           \
    firstprivate(data_pattern, temperature, random_generator, shift)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          position = data_pattern.get_index_site();
          contribution =
              contribution_site(gauge_angles, conf_angles, data_pattern, x, y,
                                z, t, position, shift);
          heat_bath(gauge_angles[position], contribution, temperature,
                    random_generator);
        }
      }
    }
  }
}

void heat_bath_update_fast(std::vector<std::complex<double>> &gauge_angles,
                           const std::vector<std::complex<double>> &conf_angles,
                           DataPatternLexicographical &data_pattern,
                           double temperature) {
  std::vector<int> shift = {1, data_pattern.lat_dim[0],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2] *
                                data_pattern.lat_dim[3]};
  unsigned seed = time(NULL);
  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);
  std::complex<double> contribution;
  int position;
#pragma omp parallel for collapse(4) private(contribution, position)           \
    firstprivate(data_pattern, temperature, random_generator, shift)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          position = data_pattern.get_index_site();
          contribution =
              contribution_site(gauge_angles, conf_angles, data_pattern, x, y,
                                z, t, position, shift);
          heat_bath_fast(gauge_angles[position], contribution, temperature,
                         random_generator);
        }
      }
    }
  }
}

void overrelaxation_update(
    std::vector<std::complex<double>> &gauge_complex,
    const std::vector<std::complex<double>> &conf_complex,
    DataPatternLexicographical &data_pattern) {
  std::vector<int> shift = {1, data_pattern.lat_dim[0],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2] *
                                data_pattern.lat_dim[3]};
  int position;
  std::complex<double> contribution;
  double alpha, phi;
#pragma omp parallel for collapse(4) private(                                  \
        contribution, position, alpha, phi) firstprivate(data_pattern, shift)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          position = data_pattern.get_index_site();
          contribution =
              contribution_site(gauge_complex, conf_complex, data_pattern, x, y,
                                z, t, position, shift);
          phi = std::atan2(contribution.imag(), contribution.real());
          alpha = std::atan2(gauge_complex[position].imag(),
                             gauge_complex[position].real());
          alpha = -alpha - 2 * phi;
          gauge_complex[position] =
              std::complex<double>(cos(alpha), sin(alpha));
        }
      }
    }
  }
}

std::tuple<double, double>
relaxation_update(std::vector<std::complex<double>> &gauge_complex,
                  const std::vector<std::complex<double>> &conf_complex,
                  DataPatternLexicographical &data_pattern) {
  std::vector<int> shift = {1, data_pattern.lat_dim[0],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2] *
                                data_pattern.lat_dim[3]};
  int position;
  std::complex<double> contribution;
  double alpha, phi;
  double diff_average = 0;
  double diff_maximal = 0;
  double diff_tmp;
#pragma omp parallel for collapse(4) private(contribution, position, alpha,    \
                                                 diff_tmp)                     \
    firstprivate(data_pattern, shift) reduction(+ : diff_average)              \
    reduction(max : diff_maximal)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          position = data_pattern.get_index_site();
          contribution =
              contribution_site(gauge_complex, conf_complex, data_pattern, x, y,
                                z, t, position, shift);
          phi = std::atan2(gauge_complex[position].imag(),
                           gauge_complex[position].real());
          alpha = std::atan2(contribution.imag(), contribution.real());
          diff_tmp = 1 - cos(alpha + phi);
          if (diff_tmp > diff_maximal)
            diff_maximal = diff_tmp;
          diff_average += diff_tmp;
          gauge_complex[position] =
              std::complex<double>(cos(alpha), -sin(alpha));
        }
      }
    }
  }
  diff_average = diff_average / data_pattern.get_lattice_size();
  return std::tuple<double, double>(diff_maximal, diff_average);
}

void heat_bath_test1(abelian &g, const abelian &K, double temperature,
                     double &random_numbers) {
  double temp_factor;
  long double temp_exponent;
  spin spin_new;

  // generation of new spin variable
  temp_factor = K.r / temperature;
  double x;
  if (temp_factor < 700) {
    temp_exponent = exp(temp_factor);

    // generation of random variable with distribution p ~ exp(ax)
    x = log(random_numbers * (temp_exponent - 1 / temp_exponent) +
            1 / temp_exponent) /
        temp_factor;
  } else
    x = log(random_numbers) / temp_factor + 1;

  g = abelian(1, acos(x) - K.phi);
}

abelian contribution_site_test1(const std::vector<abelian> &gauge_abelian,
                                const std::vector<abelian> &conf_abelian, int x,
                                int y, int z, int t, int position,
                                std::vector<int> &shift) {

  abelian A(0, 0);

  // mu = 0
  if (x < x_size - 1)
    A = A +
        conf_abelian[position * 4] * gauge_abelian[position + shift[0]].conj();
  else
    A = A + conf_abelian[position * 4] *
                gauge_abelian[position + shift[0] - shift[1]].conj();

  if (x > 0)
    A = A + conf_abelian[(position - shift[0]) * 4].conj() *
                gauge_abelian[position - shift[0]].conj();
  else
    A = A + conf_abelian[(position - shift[0] + shift[1]) * 4].conj() *
                gauge_abelian[position - shift[0] + shift[1]].conj();

  // mu = 1
  if (y < y_size - 1)
    A = A + conf_abelian[position * 4 + 1] *
                gauge_abelian[position + shift[1]].conj();
  else
    A = A + conf_abelian[position * 4 + 1] *
                gauge_abelian[position + shift[1] - shift[2]].conj();

  if (y > 0)
    A = A + conf_abelian[(position - shift[1]) * 4 + 1].conj() *
                gauge_abelian[position - shift[1]].conj();
  else
    A = A + conf_abelian[(position - shift[1] + shift[2]) * 4 + 1].conj() *
                gauge_abelian[position - shift[1] + shift[2]].conj();

  // mu = 2
  if (z < z_size - 1)
    A = A + conf_abelian[position * 4 + 2] *
                gauge_abelian[position + shift[2]].conj();
  else
    A = A + conf_abelian[position * 4 + 2] *
                gauge_abelian[position + shift[2] - shift[3]].conj();

  if (z > 0)
    A = A + conf_abelian[(position - shift[2]) * 4 + 2].conj() *
                gauge_abelian[position - shift[2]].conj();
  else
    A = A + conf_abelian[(position - shift[2] + shift[3]) * 4 + 2].conj() *
                gauge_abelian[position - shift[2] + shift[3]].conj();

  // mu = 3

  if (t < t_size - 1)
    A = A + conf_abelian[position * 4 + 3] *
                gauge_abelian[position + shift[3]].conj();

  else
    A = A + conf_abelian[position * 4 + 3] *
                gauge_abelian[position + shift[3] - shift[4]].conj();

  if (t > 0)
    A = A + conf_abelian[(position - shift[3]) * 4 + 3].conj() *
                gauge_abelian[position - shift[3]].conj();
  else
    A = A + conf_abelian[(position - shift[3] + shift[4]) * 4 + 3].conj() *
                gauge_abelian[position - shift[3] + shift[4]].conj();

  return A;
}

void heat_bath_update_test1(std::vector<abelian> &gauge_abelian,
                            const std::vector<abelian> &conf_abelian,
                            double temperature) {

  std::vector<double> random_numbers;

  // generate random numbers for heat bath
  random_numbers = generate_random_numbers1(x_size * y_size * z_size * t_size);

  std::vector<int> shift = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  abelian A;
  int position = 0;

  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {

          A = contribution_site_test1(gauge_abelian, conf_abelian, x, y, z, t,
                                      position, shift);

          heat_bath_test1(gauge_abelian[position], A, temperature,
                          random_numbers[position]);

          position++;
        }
      }
    }
  }

  normalize_abelian(gauge_abelian);
}

void normalize_abelian(std::vector<abelian> &abelian) {
  for (int i = 0; i < abelian.size(); i++) {
    abelian[i].r = 1.;
  }
}

void normalize_complex(std::vector<std::complex<double>> &gauge_complex) {
  double norm;
  for (int i = 0; i < gauge_complex.size(); i++) {
    norm = sqrt(gauge_complex[i].real() * gauge_complex[i].real() +
                gauge_complex[i].imag() * gauge_complex[i].imag());
    gauge_complex[i] = gauge_complex[i] / norm;
  }
}

double
Landau_functional_gauge_abelian(const std::vector<abelian> &gauge_abelian,
                                const std::vector<abelian> &conf_abelian) {
  link1 link(x_size, y_size, z_size, t_size);
  double result = 0;
  abelian tmp;

  SPACE_ITER_START

  for (int mu = 0; mu < 4; mu++) {

    tmp = gauge_abelian[link.place / 4] * conf_abelian[link.place + mu];
    link.move(mu, 1);
    result += (tmp * gauge_abelian[link.place / 4].conj()).tr();
    link.move(mu, -1);
  }

  SPACE_ITER_END

  return result / (x_size * y_size * z_size * t_size * 4);
}

double
Landau_functional_gauge_angles(const std::vector<double> &gauge_angles,
                               const std::vector<abelian> &conf_abelian) {
  link1 link(x_size, y_size, z_size, t_size);
  double result = 0;
  abelian tmp;

  SPACE_ITER_START

  for (int mu = 0; mu < 4; mu++) {

    tmp.r = conf_abelian[link.place + mu].r;
    tmp.phi = conf_abelian[link.place + mu].phi + gauge_angles[link.place / 4];

    link.move(mu, 1);
    tmp.phi -= gauge_angles[link.place / 4];
    result += tmp.tr();
    link.move(mu, -1);
  }

  SPACE_ITER_END

  return result / (x_size * y_size * z_size * t_size * 4);
}

double Landau_functional_gauge(const std::vector<double> &gauge_angles,
                               const std::vector<double> &conf_angles) {
  link1 link(x_size, y_size, z_size, t_size);
  double result = 0;
  double angle_tmp;

  SPACE_ITER_START

  for (int mu = 0; mu < 4; mu++) {

    angle_tmp = conf_angles[link.place + mu] + gauge_angles[link.place / 4];

    link.move(mu, 1);
    result += cos(angle_tmp - gauge_angles[link.place / 4]);
    link.move(mu, -1);
  }

  SPACE_ITER_END

  return result / (x_size * y_size * z_size * t_size * 4);
}

double Landau_functional_gauge_complex(
    const std::vector<std::complex<double>> &gauge_complex,
    const std::vector<std::complex<double>> &conf_complex) {
  link1 link(x_size, y_size, z_size, t_size);
  double result = 0;
  std::complex<double> tmp;

  SPACE_ITER_START

  for (int mu = 0; mu < 4; mu++) {

    tmp = conf_complex[link.place + mu] * gauge_complex[link.place / 4];

    link.move(mu, 1);

    result += (tmp * std::conj(gauge_complex[link.place / 4])).real();

    link.move(mu, -1);
  }

  SPACE_ITER_END

  return result / (x_size * y_size * z_size * t_size * 4);
}

double Landau_functional(const std::vector<su2> &conf_su2) {
  double result = 0;
  for (int i = 0; i < conf_su2.size(); i++) {
    result += conf_su2[i].a0 / sqrt(conf_su2[i].a0 * conf_su2[i].a0 +
                                    conf_su2[i].a3 * conf_su2[i].a3);
  }
  return result / conf_su2.size();
}

double Landau_functional_abelian(std::vector<abelian> &conf_abelian) {
  double result = 0;
  for (int i = 0; i < conf_abelian.size(); i++) {
    result += conf_abelian[i].tr();
  }
  return result / conf_abelian.size();
}

double Landau_functional_complex(
    const std::vector<std::complex<double>> &conf_complex) {
  double result = 0;
  for (int i = 0; i < conf_complex.size(); i++) {
    result += conf_complex[i].real();
  }
  return result / conf_complex.size();
}

double Landau_functional_conf_complex(
    const std::vector<std::complex<double>> &conf_complex,
    const std::vector<std::complex<double>> &gauge_complex,
    DataPatternLexicographical &data_pattern) {
  double result = 0;
  int index_gauge;
  int index_conf;
#pragma omp parallel for collapse(4) private(index_gauge, index_conf)          \
    firstprivate(data_pattern) reduction(+ : result)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index_gauge = data_pattern.get_index_site();
          for (int mu = 0; mu < 4; mu++) {
            index_conf = data_pattern.get_index_link(mu);
            data_pattern.move_forward(1, mu);
            result += (gauge_complex[index_gauge] * conf_complex[index_conf] *
                       std::conj(gauge_complex[data_pattern.get_index_site()]))
                          .real();
            data_pattern.move_backward(1, mu);
          }
        }
      }
    }
  }
  return result / (data_pattern.get_data_size());
}

void gauge_tranformation_abelian(std::vector<abelian> &gauge_abelian,
                                 std::vector<abelian> &conf_abelian) {
  link1 link(x_size, y_size, z_size, t_size);
  abelian tmp;
  SPACE_ITER_START
  for (int mu = 0; mu < 4; mu++) {
    tmp = gauge_abelian[link.place / 4] * conf_abelian[link.place + mu];
    link.move(mu, 1);
    tmp = tmp * gauge_abelian[link.place / 4].conj();
    link.move(mu, -1);
    conf_abelian[link.place + mu] = tmp;
  }
  SPACE_ITER_END
}

void gauge_tranformation(std::vector<abelian> &gauge_abelian,
                         std::vector<su2> &conf_su2) {
  link1 link(x_size, y_size, z_size, t_size);
  su2 tmp;
  SPACE_ITER_START
  for (int mu = 0; mu < 4; mu++) {
    tmp = su2(cos(gauge_abelian[link.place / 4].phi), 0., 0.,
              sin(gauge_abelian[link.place / 4].phi));

    conf_su2[link.place + mu] = tmp * conf_su2[link.place + mu];

    link.move(mu, -1);

    conf_su2[link.place + mu] = conf_su2[link.place + mu] ^ tmp;

    link.move(mu, 1);
  }
  SPACE_ITER_END
}

void overrelaxation_update_test1(std::vector<abelian> &gauge_abelian,
                                 std::vector<abelian> &conf_abelian) {

  std::vector<int> shift = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  abelian A(0, 0);
  int position = 0;
  int count = 0;

  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {

          A = contribution_site_test1(gauge_abelian, conf_abelian, x, y, z, t,
                                      position, shift);

          gauge_abelian[position].phi =
              -gauge_abelian[position].phi - 2 * A.phi;

          position++;
          count++;
        }
      }
    }
  }

  normalize_abelian(gauge_abelian);
}

std::tuple<double, double>
relaxation_update_test1(std::vector<abelian> &gauge_abelian,
                        std::vector<abelian> &conf_abelian) {

  std::vector<int> shift = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  abelian A(0, 0);
  int position = 0;
  int count = 0;

  double diff_average = 0;
  double diff_maximal = 0;
  double diff_tmp;

  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {

          A = contribution_site_test1(gauge_abelian, conf_abelian, x, y, z, t,
                                      position, shift);

          diff_tmp = gauge_abelian[position].phi;
          gauge_abelian[position].phi = -A.phi;
          diff_tmp = 1 - cos(gauge_abelian[position].phi - diff_tmp);

          if (diff_tmp > diff_maximal)
            diff_maximal = diff_tmp;
          diff_average += diff_tmp;

          position++;
          count++;
        }
      }
    }
  }

  diff_average = diff_average / (x_size * y_size * z_size * t_size);

  normalize_abelian(gauge_abelian);

  return std::tuple<double, double>(diff_maximal, diff_average);
}

void make_simulated_annealing_test1(std::vector<abelian> &gauge_abelian,
                                    std::vector<abelian> &conf_abelian,
                                    double T_init, double T_final,
                                    double T_step, int OR_steps,
                                    int thermalization_steps) {

  for (int i = 0; i < thermalization_steps; i++) {
    heat_bath_update_test1(gauge_abelian, conf_abelian, T_init);
  }

  double T = T_init;

  while (T > T_final) {

    heat_bath_update_test1(gauge_abelian, conf_abelian, T);

    for (int i = 0; i < OR_steps; i++) {
      heat_bath_update_test1(gauge_abelian, conf_abelian, T);
    }

    if (T <= 1.4 && T >= 0.8)
      T -= T_step / 4;
    else
      T -= T_step;
  }
}

void make_maximization_approximate(std::vector<abelian> &gauge_abelian,
                                   std::vector<abelian> &conf_abelian,
                                   int OR_steps, int tolerance_digits) {

  double functional =
      Landau_functional_gauge_abelian(gauge_abelian, conf_abelian);
  double functional_old = functional;

  bool is_equal = false;

  do {

    relaxation_update_test1(gauge_abelian, conf_abelian);

    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update_test1(gauge_abelian, conf_abelian);
    }

    functional = Landau_functional_gauge_abelian(gauge_abelian, conf_abelian);

    is_equal = (trunc(powf(10., tolerance_digits) * functional) ==
                trunc(powf(10., tolerance_digits) * functional_old)) &&
               (signbit(functional) == signbit(functional_old));

    functional_old = functional;

  } while (!is_equal);
}

void make_maximization_final(std::vector<abelian> &gauge_abelian,
                             std::vector<abelian> &conf_abelian, int OR_steps,
                             double tolerance_maximal,
                             double tolerance_average) {

  double functional =
      Landau_functional_gauge_abelian(gauge_abelian, conf_abelian);
  double functional_old = functional;

  bool is_converged = false;
  std::tuple<double, double> difference;

  do {

    difference = relaxation_update_test1(gauge_abelian, conf_abelian);

    is_converged = (std::get<0>(difference) < tolerance_maximal) &&
                   (std::get<1>(difference) < tolerance_average);

    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update_test1(gauge_abelian, conf_abelian);
    }

  } while (!is_converged);
}

void heat_bath_test2(double &g, abelian &K, double temperature,
                     double &random_numbers) {
  double temp_factor;
  long double temp_exponent;
  spin spin_new;

  // generation of new spin variable
  temp_factor = K.r / temperature;
  double x;
  if (temp_factor < 700) {
    temp_exponent = exp(temp_factor);

    // generation of random variable with distribution p ~ exp(ax)
    x = log(random_numbers * (temp_exponent - 1 / temp_exponent) +
            1 / temp_exponent) /
        temp_factor;
  } else
    x = log(random_numbers) / temp_factor + 1;

  g = acos(x) - K.phi;
}

abelian contribution_site_test2(std::vector<double> &gauge_abelian,
                                std::vector<abelian> &conf_abelian, int x,
                                int y, int z, int t, int position,
                                std::vector<int> &shift) {

  std::complex<double> complex_tmp(0, 0);

  // mu = 0
  if (x < x_size - 1) {
    complex_tmp += std::complex<double>(
        conf_abelian[position * 4].r * cos(conf_abelian[position * 4].phi -
                                           gauge_abelian[position + shift[0]]),
        conf_abelian[position * 4].r * sin(conf_abelian[position * 4].phi -
                                           gauge_abelian[position + shift[0]]));
  } else {
    complex_tmp += std::complex<double>(
        conf_abelian[position * 4].r *
            cos(conf_abelian[position * 4].phi -
                gauge_abelian[position + shift[0] - shift[1]]),
        conf_abelian[position * 4].r *
            sin(conf_abelian[position * 4].phi -
                gauge_abelian[position + shift[0] - shift[1]]));
  }

  if (x > 0) {
    complex_tmp += std::complex<double>(
        conf_abelian[(position - shift[0]) * 4].r *
            cos(-conf_abelian[(position - shift[0]) * 4].phi -
                gauge_abelian[position - shift[0]]),
        conf_abelian[(position - shift[0]) * 4].r *
            sin(-conf_abelian[(position - shift[0]) * 4].phi -
                gauge_abelian[position - shift[0]]));
  } else {
    complex_tmp += std::complex<double>(
        conf_abelian[(position - shift[0] + shift[1]) * 4].r *
            cos(-conf_abelian[(position - shift[0] + shift[1]) * 4].phi -
                gauge_abelian[position - shift[0] + shift[1]]),
        conf_abelian[(position - shift[0] + shift[1]) * 4].r *
            sin(-conf_abelian[(position - shift[0] + shift[1]) * 4].phi -
                gauge_abelian[position - shift[0] + shift[1]]));
  }

  // mu = 1

  if (y < y_size - 1) {
    complex_tmp +=
        std::complex<double>(conf_abelian[position * 4 + 1].r *
                                 cos(conf_abelian[position * 4 + 1].phi -
                                     gauge_abelian[position + shift[1]]),
                             conf_abelian[position * 4 + 1].r *
                                 sin(conf_abelian[position * 4 + 1].phi -
                                     gauge_abelian[position + shift[1]]));
  } else {
    complex_tmp += std::complex<double>(
        conf_abelian[position * 4 + 1].r *
            cos(conf_abelian[position * 4 + 1].phi -
                gauge_abelian[position + shift[1] - shift[2]]),
        conf_abelian[position * 4 + 1].r *
            sin(conf_abelian[position * 4 + 1].phi -
                gauge_abelian[position + shift[1] - shift[2]]));
  }

  if (y > 0) {
    complex_tmp += std::complex<double>(
        conf_abelian[(position - shift[1]) * 4 + 1].r *
            cos(-conf_abelian[(position - shift[1]) * 4 + 1].phi -
                gauge_abelian[position - shift[1]]),
        conf_abelian[(position - shift[1]) * 4 + 1].r *
            sin(-conf_abelian[(position - shift[1]) * 4 + 1].phi -
                gauge_abelian[position - shift[1]]));
  } else {
    complex_tmp += std::complex<double>(
        conf_abelian[(position - shift[1] + shift[2]) * 4 + 1].r *
            cos(-conf_abelian[(position - shift[1] + shift[2]) * 4 + 1].phi -
                gauge_abelian[position - shift[1] + shift[2]]),
        conf_abelian[(position - shift[1] + shift[2]) * 4 + 1].r *
            sin(-conf_abelian[(position - shift[1] + shift[2]) * 4 + 1].phi -
                gauge_abelian[position - shift[1] + shift[2]]));
  }

  // mu = 2

  if (z < z_size - 1) {
    complex_tmp +=
        std::complex<double>(conf_abelian[position * 4 + 2].r *
                                 cos(conf_abelian[position * 4 + 2].phi -
                                     gauge_abelian[position + shift[2]]),
                             conf_abelian[position * 4 + 2].r *
                                 sin(conf_abelian[position * 4 + 2].phi -
                                     gauge_abelian[position + shift[2]]));
  } else {
    complex_tmp += std::complex<double>(
        conf_abelian[position * 4 + 2].r *
            cos(conf_abelian[position * 4 + 2].phi -
                gauge_abelian[position + shift[2] - shift[3]]),
        conf_abelian[position * 4 + 2].r *
            sin(conf_abelian[position * 4 + 2].phi -
                gauge_abelian[position + shift[2] - shift[3]]));
  }

  if (z > 0) {
    complex_tmp += std::complex<double>(
        conf_abelian[(position - shift[2]) * 4 + 2].r *
            cos(-conf_abelian[(position - shift[2]) * 4 + 2].phi -
                gauge_abelian[position - shift[2]]),
        conf_abelian[(position - shift[2]) * 4 + 2].r *
            sin(-conf_abelian[(position - shift[2]) * 4 + 2].phi -
                gauge_abelian[position - shift[2]]));
  } else {
    complex_tmp += std::complex<double>(
        conf_abelian[(position - shift[2] + shift[3]) * 4 + 2].r *
            cos(-conf_abelian[(position - shift[2] + shift[3]) * 4 + 2].phi -
                gauge_abelian[position - shift[2] + shift[3]]),
        conf_abelian[(position - shift[2] + shift[3]) * 4 + 2].r *
            sin(-conf_abelian[(position - shift[2] + shift[3]) * 4 + 2].phi -
                gauge_abelian[position - shift[2] + shift[3]]));
  }

  // mu = 3

  if (t < t_size - 1) {
    complex_tmp +=
        std::complex<double>(conf_abelian[position * 4 + 3].r *
                                 cos(conf_abelian[position * 4 + 3].phi -
                                     gauge_abelian[position + shift[3]]),
                             conf_abelian[position * 4 + 3].r *
                                 sin(conf_abelian[position * 4 + 3].phi -
                                     gauge_abelian[position + shift[3]]));
  } else {
    complex_tmp += std::complex<double>(
        conf_abelian[position * 4 + 3].r *
            cos(conf_abelian[position * 4 + 3].phi -
                gauge_abelian[position + shift[3] - shift[4]]),
        conf_abelian[position * 4 + 3].r *
            sin(conf_abelian[position * 4 + 3].phi -
                gauge_abelian[position + shift[3] - shift[4]]));
  }

  if (t > 0) {
    complex_tmp += std::complex<double>(
        conf_abelian[(position - shift[3]) * 4 + 3].r *
            cos(-conf_abelian[(position - shift[3]) * 4 + 3].phi -
                gauge_abelian[position - shift[3]]),
        conf_abelian[(position - shift[3]) * 4 + 3].r *
            sin(-conf_abelian[(position - shift[3]) * 4 + 3].phi -
                gauge_abelian[position - shift[3]]));
  } else {
    complex_tmp += std::complex<double>(
        conf_abelian[(position - shift[3] + shift[4]) * 4 + 3].r *
            cos(-conf_abelian[(position - shift[3] + shift[4]) * 4 + 3].phi -
                gauge_abelian[position - shift[3] + shift[4]]),
        conf_abelian[(position - shift[3] + shift[4]) * 4 + 3].r *
            sin(-conf_abelian[(position - shift[3] + shift[4]) * 4 + 3].phi -
                gauge_abelian[position - shift[3] + shift[4]]));
  }

  return abelian(sqrt(complex_tmp.real() * complex_tmp.real() +
                      complex_tmp.imag() * complex_tmp.imag()),
                 atan2(complex_tmp.imag(), complex_tmp.real()));
}

void heat_bath_update_test2(std::vector<double> &gauge_angles,
                            std::vector<abelian> &conf_abelian,
                            double temperature) {

  std::vector<double> random_numbers;

  // generate random numbers for heat bath
  random_numbers = generate_random_numbers1(x_size * y_size * z_size * t_size);

  std::vector<int> shift = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  abelian A;
  int position = 0;

  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {

          A = contribution_site_test2(gauge_angles, conf_abelian, x, y, z, t,
                                      position, shift);

          heat_bath_test2(gauge_angles[position], A, temperature,
                          random_numbers[position]);

          position++;
        }
      }
    }
  }
}

abelian contribution_site_test3(std::vector<double> &gauge_angles,
                                std::vector<double> &conf_angles, int x, int y,
                                int z, int t, int position,
                                std::vector<int> &shift) {

  std::complex<double> complex_tmp(0, 0);

  // mu = 0
  if (x < x_size - 1) {
    complex_tmp += std::complex<double>(
        cos(conf_angles[position * 4] - gauge_angles[position + shift[0]]),
        sin(conf_angles[position * 4] - gauge_angles[position + shift[0]]));
  } else {
    complex_tmp +=
        std::complex<double>(cos(conf_angles[position * 4] -
                                 gauge_angles[position + shift[0] - shift[1]]),
                             sin(conf_angles[position * 4] -
                                 gauge_angles[position + shift[0] - shift[1]]));
  }

  if (x > 0) {
    complex_tmp +=
        std::complex<double>(cos(-conf_angles[(position - shift[0]) * 4] -
                                 gauge_angles[position - shift[0]]),
                             sin(-conf_angles[(position - shift[0]) * 4] -
                                 gauge_angles[position - shift[0]]));
  }

  // mu = 1

  if (y < y_size - 1) {
    complex_tmp += std::complex<double>(
        cos(conf_angles[position * 4 + 1] - gauge_angles[position + shift[1]]),
        sin(conf_angles[position * 4 + 1] - gauge_angles[position + shift[1]]));
  } else {
    complex_tmp +=
        std::complex<double>(cos(conf_angles[position * 4 + 1] -
                                 gauge_angles[position + shift[1] - shift[2]]),
                             sin(conf_angles[position * 4 + 1] -
                                 gauge_angles[position + shift[1] - shift[2]]));
  }

  if (y > 0) {
    complex_tmp +=
        std::complex<double>(cos(-conf_angles[(position - shift[1]) * 4 + 1] -
                                 gauge_angles[position - shift[1]]),
                             sin(-conf_angles[(position - shift[1]) * 4 + 1] -
                                 gauge_angles[position - shift[1]]));

  } else {
    complex_tmp += std::complex<double>(
        cos(-conf_angles[(position - shift[1] + shift[2]) * 4 + 1] -
            gauge_angles[position - shift[1] + shift[2]]),
        sin(-conf_angles[(position - shift[1] + shift[2]) * 4 + 1] -
            gauge_angles[position - shift[1] + shift[2]]));
  }

  // mu = 2

  if (z < z_size - 1) {
    complex_tmp += std::complex<double>(
        cos(conf_angles[position * 4 + 2] - gauge_angles[position + shift[2]]),
        sin(conf_angles[position * 4 + 2] - gauge_angles[position + shift[2]]));
  } else {
    complex_tmp +=
        std::complex<double>(cos(conf_angles[position * 4 + 2] -
                                 gauge_angles[position + shift[2] - shift[3]]),
                             sin(conf_angles[position * 4 + 2] -
                                 gauge_angles[position + shift[2] - shift[3]]));
  }

  if (z > 0) {
    complex_tmp +=
        std::complex<double>(cos(-conf_angles[(position - shift[2]) * 4 + 2] -
                                 gauge_angles[position - shift[2]]),
                             sin(-conf_angles[(position - shift[2]) * 4 + 2] -
                                 gauge_angles[position - shift[2]]));
  } else {
    complex_tmp += std::complex<double>(
        cos(-conf_angles[(position - shift[2] + shift[3]) * 4 + 2] -
            gauge_angles[position - shift[2] + shift[3]]),
        sin(-conf_angles[(position - shift[2] + shift[3]) * 4 + 2] -
            gauge_angles[position - shift[2] + shift[3]]));
  }

  // mu = 3

  if (t < t_size - 1) {
    complex_tmp += std::complex<double>(
        cos(conf_angles[position * 4 + 3] - gauge_angles[position + shift[3]]),
        sin(conf_angles[position * 4 + 3] - gauge_angles[position + shift[3]]));
  } else {
    complex_tmp +=
        std::complex<double>(cos(conf_angles[position * 4 + 3] -
                                 gauge_angles[position + shift[3] - shift[4]]),
                             sin(conf_angles[position * 4 + 3] -
                                 gauge_angles[position + shift[3] - shift[4]]));
  }

  if (t > 0) {
    complex_tmp +=
        std::complex<double>(cos(-conf_angles[(position - shift[3]) * 4 + 3] -
                                 gauge_angles[position - shift[3]]),
                             sin(-conf_angles[(position - shift[3]) * 4 + 3] -
                                 gauge_angles[position - shift[3]]));
  } else {
    complex_tmp += std::complex<double>(
        cos(-conf_angles[(position - shift[3] + shift[4]) * 4 + 3] -
            gauge_angles[position - shift[3] + shift[4]]),
        sin(-conf_angles[(position - shift[3] + shift[4]) * 4 + 3] -
            gauge_angles[position - shift[3] + shift[4]]));
  }

  return abelian(sqrt(complex_tmp.real() * complex_tmp.real() +
                      complex_tmp.imag() * complex_tmp.imag()),
                 atan2(complex_tmp.imag(), complex_tmp.real()));
}

void heat_bath_update_test3(std::vector<double> &gauge_angles,
                            std::vector<double> &conf_angles,
                            double temperature) {

  std::vector<double> random_numbers;

  // generate random numbers for heat bath
  random_numbers = generate_random_numbers1(x_size * y_size * z_size * t_size);

  std::vector<int> shift = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  abelian A;
  int position = 0;

  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {

          A = contribution_site_test3(gauge_angles, conf_angles, x, y, z, t,
                                      position, shift);

          heat_bath_test2(gauge_angles[position], A, temperature,
                          random_numbers[position]);

          position++;
        }
      }
    }
  }
}

void heat_bath(std::complex<double> &g, std::complex<double> &K,
               double temperature, double *random_numbers) {
  double temp_factor;
  long double temp_exponent;
  spin spin_new;

  // generation of new spin variable
  temp_factor = sqrt(K.real() * K.real() + K.imag() * K.imag()) / temperature;
  double x;
  if (temp_factor < 700) {
    temp_exponent = exp(temp_factor);

    // generation of random variable with distribution p ~ exp(ax)
    x = log(random_numbers[0] * (temp_exponent - 1 / temp_exponent) +
            1 / temp_exponent) /
        temp_factor;
  } else
    x = log(random_numbers[0]) / temp_factor + 1;

  double angle_new;

  if (random_numbers[1] > 0)
    angle_new = acos(x) - atan2(K.imag(), K.real());
  else
    angle_new = -acos(x) - atan2(K.imag(), K.real());

  g = std::complex<double>(cos(angle_new), sin(angle_new));
}

std::complex<double>
contribution_site(std::vector<std::complex<double>> &gauge_complex,
                  std::vector<std::complex<double>> &conf_complex, int x, int y,
                  int z, int t, int position, std::vector<int> &shift) {

  std::complex<double> complex_tmp(0, 0);

  // mu = 0
  if (x < x_size - 1) {
    complex_tmp = complex_tmp + (conf_complex[position * 4] *
                                 std::conj(gauge_complex[position + shift[0]]));
  } else {
    complex_tmp = complex_tmp +
                  (conf_complex[position * 4] *
                   std::conj(gauge_complex[position + shift[0] - shift[1]]));
  }

  if (x > 0) {
    complex_tmp =
        complex_tmp + std::conj(conf_complex[(position - shift[0]) * 4] *
                                gauge_complex[position - shift[0]]);
  } else {
    complex_tmp = complex_tmp +
                  std::conj(conf_complex[(position - shift[0] + shift[1]) * 4] *
                            gauge_complex[position - shift[0] + shift[1]]);
  }

  // mu = 1

  if (y < y_size - 1) {
    complex_tmp = complex_tmp + (conf_complex[position * 4 + 1] *
                                 std::conj(gauge_complex[position + shift[1]]));
  } else {
    complex_tmp = complex_tmp +
                  (conf_complex[position * 4 + 1] *
                   std::conj(gauge_complex[position + shift[1] - shift[2]]));
  }

  if (y > 0) {
    complex_tmp =
        complex_tmp + std::conj(conf_complex[(position - shift[1]) * 4 + 1] *
                                gauge_complex[position - shift[1]]);
  } else {
    complex_tmp =
        complex_tmp +
        std::conj(conf_complex[(position - shift[1] + shift[2]) * 4 + 1] *
                  gauge_complex[position - shift[1] + shift[2]]);
  }

  // mu = 2

  if (z < z_size - 1) {
    complex_tmp = complex_tmp + (conf_complex[position * 4 + 2] *
                                 std::conj(gauge_complex[position + shift[2]]));
  } else {
    complex_tmp = complex_tmp +
                  (conf_complex[position * 4 + 2] *
                   std::conj(gauge_complex[position + shift[2] - shift[3]]));
  }

  if (z > 0) {
    complex_tmp =
        complex_tmp + std::conj(conf_complex[(position - shift[2]) * 4 + 2] *
                                gauge_complex[position - shift[2]]);
  } else {
    complex_tmp =
        complex_tmp +
        std::conj(conf_complex[(position - shift[2] + shift[3]) * 4 + 2] *
                  gauge_complex[position - shift[2] + shift[3]]);
  }

  // mu = 3

  if (t < t_size - 1) {
    complex_tmp = complex_tmp + (conf_complex[position * 4 + 3] *
                                 std::conj(gauge_complex[position + shift[3]]));
  } else {
    complex_tmp = complex_tmp +
                  (conf_complex[position * 4 + 3] *
                   std::conj(gauge_complex[position + shift[3] - shift[4]]));
  }

  if (t > 0) {
    complex_tmp =
        complex_tmp + std::conj(conf_complex[(position - shift[3]) * 4 + 3] *
                                gauge_complex[position - shift[3]]);
  } else {
    complex_tmp =
        complex_tmp +
        std::conj(conf_complex[(position - shift[3] + shift[4]) * 4 + 3] *
                  gauge_complex[position - shift[3] + shift[4]]);
  }

  return complex_tmp;
}

void heat_bath_update(std::vector<std::complex<double>> &gauge_complex,
                      std::vector<std::complex<double>> &conf_complex,
                      double temperature) {

  std::vector<double> random_numbers;

  // generate random numbers for heat bath
  random_numbers =
      generate_random_numbers1(2 * x_size * y_size * z_size * t_size);

  std::vector<int> shift = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  std::complex<double> A;
  int position = 0;

  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {

          A = contribution_site(gauge_complex, conf_complex, x, y, z, t,
                                position, shift);

          heat_bath(gauge_complex[position], A, temperature,
                    &random_numbers[2 * position]);

          position++;
        }
      }
    }
  }
}

void overrelaxation_update(std::vector<std::complex<double>> &gauge_complex,
                           std::vector<std::complex<double>> &conf_complex) {

  std::vector<int> shift = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  std::complex<double> A;
  int position = 0;
  double A_phi;
  double gauge_phi;

  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {

          A = contribution_site(gauge_complex, conf_complex, x, y, z, t,
                                position, shift);

          A_phi = atan2(A.imag(), A.real());
          gauge_phi = -atan2(gauge_complex[position].imag(),
                             gauge_complex[position].real()) -
                      2 * A_phi;

          gauge_complex[position] =
              std::complex<double>(cos(gauge_phi), sin(gauge_phi));
          position++;
        }
      }
    }
  }
}

std::tuple<double, double>
relaxation_update(std::vector<std::complex<double>> &gauge_complex,
                  std::vector<std::complex<double>> &conf_complex) {

  std::vector<int> shift = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  std::complex<double> A;
  int position = 0;

  double diff_average = 0;
  double diff_maximal = 0;
  double diff_tmp;

  double gauge_phi;
  double A_phi;

  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {

          A = contribution_site(gauge_complex, conf_complex, x, y, z, t,
                                position, shift);

          gauge_phi = atan2(gauge_complex[position].imag(),
                            gauge_complex[position].real());
          A_phi = atan2(A.imag(), A.real());

          gauge_complex[position] =
              std::complex<double>(cos(A_phi), -sin(A_phi));
          diff_tmp = 1 - cos(A_phi + gauge_phi);

          if (diff_tmp > diff_maximal)
            diff_maximal = diff_tmp;
          diff_average += diff_tmp;

          position++;
        }
      }
    }
  }

  diff_average = diff_average / (x_size * y_size * z_size * t_size);

  return std::tuple<double, double>(diff_maximal, diff_average);
}

void make_simulated_annealing(std::vector<std::complex<double>> &gauge_complex,
                              std::vector<std::complex<double>> &conf_complex,
                              double T_init, double T_final, double T_step,
                              int OR_steps, int thermalization_steps) {

  for (int i = 0; i < thermalization_steps; i++) {
    heat_bath_update(gauge_complex, conf_complex, T_init);
  }

  double T = T_init;

  while (T > T_final) {

    heat_bath_update(gauge_complex, conf_complex, T);

    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(gauge_complex, conf_complex);
    }

    T -= T_step;
  }
}

void make_maximization_final(std::vector<std::complex<double>> &gauge_complex,
                             std::vector<std::complex<double>> &conf_complex,
                             int OR_steps, double tolerance_maximal,
                             double tolerance_average) {

  bool is_converged = false;
  std::tuple<double, double> difference;

  do {

    difference = relaxation_update(gauge_complex, conf_complex);

    is_converged = (std::get<0>(difference) < tolerance_maximal) &&
                   (std::get<1>(difference) < tolerance_average);

    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(gauge_complex, conf_complex);
    }

  } while (!is_converged);
}

void apply_gauge_Landau_complex(
    std::vector<std::complex<double>> &gauge_complex,
    std::vector<std::complex<double>> &conf_complex) {
  link1 link(x_size, y_size, z_size, t_size);

  std::complex<double> complex_tmp;
  double norm;

  SPACE_ITER_START

  for (int mu = 0; mu < 4; mu++) {

    complex_tmp = gauge_complex[link.place / 4] * conf_complex[link.place + mu];

    link.move(mu, 1);

    complex_tmp = complex_tmp * std::conj(gauge_complex[link.place / 4]);

    link.move(mu, -1);

    conf_complex[link.place + mu] = complex_tmp;
  }

  SPACE_ITER_END
}

void apply_gauge_Landau(std::vector<std::complex<double>> &gauge_complex,
                        std::vector<su2> &conf_su2) {
  link1 link(x_size, y_size, z_size, t_size);

  su2 A;

  SPACE_ITER_START

  A = su2(gauge_complex[link.place / 4].real(), 0, 0,
          gauge_complex[link.place / 4].imag());

  for (int mu = 0; mu < 4; mu++) {

    conf_su2[link.place + mu] = A * conf_su2[link.place + mu];
  }

  for (int mu = 0; mu < 4; mu++) {

    link.move(mu, -1);

    conf_su2[link.place + mu] = conf_su2[link.place + mu] ^ A;

    link.move(mu, 1);
  }

  SPACE_ITER_END
}

void apply_gauge_Landau_complex(
    std::vector<std::complex<double>> &gauge_complex,
    std::vector<std::complex<double>> &conf_complex,
    DataPatternLexicographical &data_pattern) {
  std::complex<double> complex_tmp;
  double norm;
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          for (int mu = 0; mu < 4; mu++) {
            complex_tmp = gauge_complex[data_pattern.get_index_site()] *
                          conf_complex[data_pattern.get_index_link(mu)];
            data_pattern.move_forward(1, mu);
            complex_tmp =
                complex_tmp *
                std::conj(gauge_complex[data_pattern.get_index_site()]);
            data_pattern.move_backward(1, mu);
            conf_complex[data_pattern.get_index_link(mu)] = complex_tmp;
          }
        }
      }
    }
  }
}

void apply_gauge_Landau(
    std::vector<std::complex<double>> &gauge_complex,
    Data::LatticeData<DataPatternLexicographical, su2> &conf_su2) {
  DataPatternLexicographical data_pattern(conf_su2.lat_dim);
  su2 A;
  int index;
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_index_site();
          A = su2(gauge_complex[index].real(), 0, 0,
                  gauge_complex[index].imag());
          for (int mu = 0; mu < 4; mu++) {
            index = data_pattern.get_index_link(mu);
            conf_su2[index] = A * conf_su2[index];
          }
          for (int mu = 0; mu < 4; mu++) {
            data_pattern.move_backward(1, mu);
            index = data_pattern.get_index_link(mu);
            conf_su2[index] = conf_su2[index] ^ A;
            data_pattern.move_forward(1, mu);
          }
        }
      }
    }
  }
}

std::map<double, double> simulated_annealing_thermalization_test(
    const std::vector<std::complex<double>> &conf_complex,
    std::vector<std::complex<double>> &gauge_complex,
    DataPatternLexicographical &data_pattern, double T_init, double T_final,
    double T_step, int OR_steps, int thermalization_steps,
    int local_thermalization_steps) {
  double omp_time;
  std::map<double, double> result;
  for (int i = 0; i < thermalization_steps; i++) {
    heat_bath_update(gauge_complex, conf_complex, data_pattern, T_init);
  }
  double T = T_init;
  while (T > 1) {
    omp_time = omp_get_wtime();
    for (int i = 0; i < local_thermalization_steps; i++) {
      heat_bath_update(gauge_complex, conf_complex, data_pattern, T);
      for (int i = 0; i < OR_steps; i++) {
        overrelaxation_update(gauge_complex, conf_complex, data_pattern);
      }
    }
    std::cout << "heat_bath_update time: " << omp_get_wtime() - omp_time
              << std::endl;
    result[T] = Landau_functional_conf_complex(conf_complex, gauge_complex,
                                               data_pattern);
    std::cout << "T: " << T << " functional: " << result[T] << std::endl;
    if (T <= 7.5 + T_step && T >= 7.2)
      T -= T_step / 20;
    else
      T -= T_step;
  }
  while (T > T_final) {
    omp_time = omp_get_wtime();
    for (int i = 0; i < local_thermalization_steps; i++) {
      heat_bath_update_fast(gauge_complex, conf_complex, data_pattern, T);
      for (int i = 0; i < OR_steps; i++) {
        overrelaxation_update(gauge_complex, conf_complex, data_pattern);
      }
    }
    std::cout << "heat_bath_update time: " << omp_get_wtime() - omp_time
              << std::endl;
    result[T] = Landau_functional_conf_complex(conf_complex, gauge_complex,
                                               data_pattern);
    std::cout << "T: " << T << " functional: " << result[T] << std::endl;
    T -= T_step;
  }
  return result;
}

void make_simulated_annealing(
    const std::vector<std::complex<double>> &conf_complex,
    std::vector<std::complex<double>> &gauge_complex,
    DataPatternLexicographical &data_pattern, double T_init, double T_final,
    double T_step, int OR_steps, int thermalization_steps) {
  for (int i = 0; i < thermalization_steps; i++) {
    heat_bath_update(gauge_complex, conf_complex, data_pattern, T_init);
  }
  double T = T_init;
  while (T > 1) {
    heat_bath_update(gauge_complex, conf_complex, data_pattern, T);
    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(gauge_complex, conf_complex, data_pattern);
    }
    if (T <= 7.5 + T_step && T >= 7.2)
      T -= T_step / 20;
    else
      T -= T_step;
  }
  while (T > T_final) {
    heat_bath_update_fast(gauge_complex, conf_complex, data_pattern, T);
    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(gauge_complex, conf_complex, data_pattern);
    }
    T -= T_step;
  }
}

void make_maximization_final(
    const std::vector<std::complex<double>> &conf_complex,
    std::vector<std::complex<double>> &gauge_complex,
    DataPatternLexicographical &data_pattern, int OR_steps,
    double tolerance_maximal, double tolerance_average) {
  bool is_converged = false;
  std::tuple<double, double> difference;
  do {
    difference = relaxation_update(gauge_complex, conf_complex, data_pattern);
    // std::cout << "difference: " << std::get<0>(difference) << " "
    //           << std::get<1>(difference) << std::endl;
    is_converged = (std::get<0>(difference) < tolerance_maximal) &&
                   (std::get<1>(difference) < tolerance_average);
    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(gauge_complex, conf_complex, data_pattern);
    }
  } while (!is_converged);
}