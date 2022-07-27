#include "../include/Landau_U1.h"
#include "../include/link.h"
#include "../include/matrix.h"

#include <fstream>
#include <math.h>
#include <random>
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

std::vector<double> convert_to_angles(std::vector<su2> conf_su2) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  std::vector<double> conf_abelian;
  conf_abelian.reserve(data_size);

  for (int i = 0; i < data_size; i++) {
    conf_abelian.push_back(atan2(conf_su2[i].a3, conf_su2[i].a0));
  }

  return conf_abelian;
}

std::vector<abelian> convert_to_abelian(std::vector<su2> conf_su2) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  std::vector<abelian> conf_abelian;
  conf_abelian.reserve(data_size);

  for (int i = 0; i < data_size; i++) {
    conf_abelian.push_back(abelian(1, atan2(conf_su2[i].a3, conf_su2[i].a0)));
  }

  return conf_abelian;
}

std::vector<complex_t> convert_to_complex(std::vector<su2> conf_su2) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  std::vector<complex_t> conf_complex;
  conf_complex.reserve(data_size);

  double module;

  for (int i = 0; i < data_size; i++) {
    module =
        sqrt(conf_su2[i].a0 * conf_su2[i].a0 + conf_su2[i].a3 * conf_su2[i].a3);
    conf_complex.push_back(
        complex_t(conf_su2[i].a0 / module, conf_su2[i].a3 / module));
  }

  return conf_complex;
}

std::vector<double>
convert_complex_to_angles(std::vector<complex_t> conf_complex) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  std::vector<double> conf_angles;
  conf_complex.reserve(data_size);

  double module;

  for (int i = 0; i < data_size; i++) {
    conf_angles.push_back(atan2(conf_complex[i].imag, conf_complex[i].real));
  }

  return conf_angles;
}

std::vector<double> generate_gauge_angles_uniform() {

  unsigned seed = time(NULL);
  // unsigned seed = 123;

  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);

  int data_szie = x_size * y_size * z_size * t_size;

  std::vector<double> gauge_abelian;
  gauge_abelian.reserve(data_szie);

  for (int i = 0; i < data_szie; i++) {

    gauge_abelian.push_back(
        (2 * (double)random_generator() / random_generator.max() - 1) * M_PI);
  }

  return gauge_abelian;
}

std::vector<abelian> generate_gauge_abelian_uniform() {

  unsigned seed = time(NULL);
  // unsigned seed = 123;

  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);

  int data_szie = x_size * y_size * z_size * t_size;

  std::vector<abelian> gauge_abelian;
  gauge_abelian.reserve(data_szie);

  for (int i = 0; i < data_szie; i++) {

    gauge_abelian.push_back(abelian(
        1,
        (2 * (double)random_generator() / random_generator.max() - 1) * M_PI));
  }

  return gauge_abelian;
}

std::vector<complex_t> generate_gauge_complex_uniform() {

  unsigned seed = time(NULL);
  // unsigned seed = 123;

  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);

  int data_szie = x_size * y_size * z_size * t_size;

  std::vector<complex_t> gauge_complex;
  gauge_complex.reserve(data_szie);

  double angle_tmp;

  for (int i = 0; i < data_szie; i++) {

    angle_tmp =
        (2 * (double)random_generator() / random_generator.max() - 1) * M_PI;

    gauge_complex.push_back(complex_t(cos(angle_tmp), sin(angle_tmp)));
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

void heat_bath_test1(abelian &g, abelian &K, double temperature,
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

abelian contribution_site_test1(std::vector<abelian> &gauge_abelian,
                                std::vector<abelian> &conf_abelian, int x,
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

void normalize_complex(std::vector<complex_t> &gauge_complex) {
  double norm;
  for (int i = 0; i < gauge_complex.size(); i++) {
    norm = sqrt(gauge_complex[i].real * gauge_complex[i].real +
                gauge_complex[i].imag * gauge_complex[i].imag);
    gauge_complex[i] = gauge_complex[i] / norm;
  }
}

double Landau_functional_gauge_abelian(std::vector<abelian> &gauge_abelian,
                                       std::vector<abelian> &conf_abelian) {
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

double Landau_functional_gauge_angles(std::vector<double> &gauge_angles,
                                      std::vector<abelian> &conf_abelian) {
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

double Landau_functional_gauge(std::vector<double> &gauge_angles,
                               std::vector<double> &conf_angles) {
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

double Landau_functional_gauge_complex(std::vector<complex_t> &gauge_complex,
                                       std::vector<complex_t> &conf_complex) {
  link1 link(x_size, y_size, z_size, t_size);
  double result = 0;
  complex_t tmp;

  SPACE_ITER_START

  for (int mu = 0; mu < 4; mu++) {

    tmp = conf_complex[link.place + mu] * gauge_complex[link.place / 4];

    link.move(mu, 1);

    result += (tmp ^ gauge_complex[link.place / 4]).real;

    link.move(mu, -1);
  }

  SPACE_ITER_END

  return result / (x_size * y_size * z_size * t_size * 4);
}

double Landau_functional(std::vector<su2> &conf_su2) {
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

double Landau_functional_complex(std::vector<complex_t> &conf_complex) {
  double result = 0;
  for (int i = 0; i < conf_complex.size(); i++) {
    result += conf_complex[i].real;
  }
  return result / conf_complex.size();
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

  complex_t complex_tmp(0, 0);

  // mu = 0
  if (x < x_size - 1) {
    complex_tmp.real +=
        conf_abelian[position * 4].r * cos(conf_abelian[position * 4].phi -
                                           gauge_abelian[position + shift[0]]);
    complex_tmp.imag +=
        conf_abelian[position * 4].r * sin(conf_abelian[position * 4].phi -
                                           gauge_abelian[position + shift[0]]);
  } else {
    complex_tmp.real += conf_abelian[position * 4].r *
                        cos(conf_abelian[position * 4].phi -
                            gauge_abelian[position + shift[0] - shift[1]]);
    complex_tmp.imag += conf_abelian[position * 4].r *
                        sin(conf_abelian[position * 4].phi -
                            gauge_abelian[position + shift[0] - shift[1]]);
  }

  if (x > 0) {
    complex_tmp.real += conf_abelian[(position - shift[0]) * 4].r *
                        cos(-conf_abelian[(position - shift[0]) * 4].phi -
                            gauge_abelian[position - shift[0]]);
    complex_tmp.imag += conf_abelian[(position - shift[0]) * 4].r *
                        sin(-conf_abelian[(position - shift[0]) * 4].phi -
                            gauge_abelian[position - shift[0]]);
  } else {
    complex_tmp.real +=
        conf_abelian[(position - shift[0] + shift[1]) * 4].r *
        cos(-conf_abelian[(position - shift[0] + shift[1]) * 4].phi -
            gauge_abelian[position - shift[0] + shift[1]]);
    complex_tmp.imag +=
        conf_abelian[(position - shift[0] + shift[1]) * 4].r *
        sin(-conf_abelian[(position - shift[0] + shift[1]) * 4].phi -
            gauge_abelian[position - shift[0] + shift[1]]);
  }

  // mu = 1

  if (y < y_size - 1) {
    complex_tmp.real += conf_abelian[position * 4 + 1].r *
                        cos(conf_abelian[position * 4 + 1].phi -
                            gauge_abelian[position + shift[1]]);
    complex_tmp.imag += conf_abelian[position * 4 + 1].r *
                        sin(conf_abelian[position * 4 + 1].phi -
                            gauge_abelian[position + shift[1]]);
  } else {
    complex_tmp.real += conf_abelian[position * 4 + 1].r *
                        cos(conf_abelian[position * 4 + 1].phi -
                            gauge_abelian[position + shift[1] - shift[2]]);
    complex_tmp.imag += conf_abelian[position * 4 + 1].r *
                        sin(conf_abelian[position * 4 + 1].phi -
                            gauge_abelian[position + shift[1] - shift[2]]);
  }

  if (y > 0) {
    complex_tmp.real += conf_abelian[(position - shift[1]) * 4 + 1].r *
                        cos(-conf_abelian[(position - shift[1]) * 4 + 1].phi -
                            gauge_abelian[position - shift[1]]);
    complex_tmp.imag += conf_abelian[(position - shift[1]) * 4 + 1].r *
                        sin(-conf_abelian[(position - shift[1]) * 4 + 1].phi -
                            gauge_abelian[position - shift[1]]);
  } else {
    complex_tmp.real +=
        conf_abelian[(position - shift[1] + shift[2]) * 4 + 1].r *
        cos(-conf_abelian[(position - shift[1] + shift[2]) * 4 + 1].phi -
            gauge_abelian[position - shift[1] + shift[2]]);
    complex_tmp.imag +=
        conf_abelian[(position - shift[1] + shift[2]) * 4 + 1].r *
        sin(-conf_abelian[(position - shift[1] + shift[2]) * 4 + 1].phi -
            gauge_abelian[position - shift[1] + shift[2]]);
  }

  // mu = 2

  if (z < z_size - 1) {
    complex_tmp.real += conf_abelian[position * 4 + 2].r *
                        cos(conf_abelian[position * 4 + 2].phi -
                            gauge_abelian[position + shift[2]]);
    complex_tmp.imag += conf_abelian[position * 4 + 2].r *
                        sin(conf_abelian[position * 4 + 2].phi -
                            gauge_abelian[position + shift[2]]);
  } else {
    complex_tmp.real += conf_abelian[position * 4 + 2].r *
                        cos(conf_abelian[position * 4 + 2].phi -
                            gauge_abelian[position + shift[2] - shift[3]]);
    complex_tmp.imag += conf_abelian[position * 4 + 2].r *
                        sin(conf_abelian[position * 4 + 2].phi -
                            gauge_abelian[position + shift[2] - shift[3]]);
  }

  if (z > 0) {
    complex_tmp.real += conf_abelian[(position - shift[2]) * 4 + 2].r *
                        cos(-conf_abelian[(position - shift[2]) * 4 + 2].phi -
                            gauge_abelian[position - shift[2]]);
    complex_tmp.imag += conf_abelian[(position - shift[2]) * 4 + 2].r *
                        sin(-conf_abelian[(position - shift[2]) * 4 + 2].phi -
                            gauge_abelian[position - shift[2]]);
  } else {
    complex_tmp.real +=
        conf_abelian[(position - shift[2] + shift[3]) * 4 + 2].r *
        cos(-conf_abelian[(position - shift[2] + shift[3]) * 4 + 2].phi -
            gauge_abelian[position - shift[2] + shift[3]]);
    complex_tmp.imag +=
        conf_abelian[(position - shift[2] + shift[3]) * 4 + 2].r *
        sin(-conf_abelian[(position - shift[2] + shift[3]) * 4 + 2].phi -
            gauge_abelian[position - shift[2] + shift[3]]);
  }

  // mu = 3

  if (t < t_size - 1) {
    complex_tmp.real += conf_abelian[position * 4 + 3].r *
                        cos(conf_abelian[position * 4 + 3].phi -
                            gauge_abelian[position + shift[3]]);
    complex_tmp.imag += conf_abelian[position * 4 + 3].r *
                        sin(conf_abelian[position * 4 + 3].phi -
                            gauge_abelian[position + shift[3]]);
  } else {
    complex_tmp.real += conf_abelian[position * 4 + 3].r *
                        cos(conf_abelian[position * 4 + 3].phi -
                            gauge_abelian[position + shift[3] - shift[4]]);
    complex_tmp.imag += conf_abelian[position * 4 + 3].r *
                        sin(conf_abelian[position * 4 + 3].phi -
                            gauge_abelian[position + shift[3] - shift[4]]);
  }

  if (t > 0) {
    complex_tmp.real += conf_abelian[(position - shift[3]) * 4 + 3].r *
                        cos(-conf_abelian[(position - shift[3]) * 4 + 3].phi -
                            gauge_abelian[position - shift[3]]);
    complex_tmp.imag += conf_abelian[(position - shift[3]) * 4 + 3].r *
                        sin(-conf_abelian[(position - shift[3]) * 4 + 3].phi -
                            gauge_abelian[position - shift[3]]);
  } else {
    complex_tmp.real +=
        conf_abelian[(position - shift[3] + shift[4]) * 4 + 3].r *
        cos(-conf_abelian[(position - shift[3] + shift[4]) * 4 + 3].phi -
            gauge_abelian[position - shift[3] + shift[4]]);
    complex_tmp.imag +=
        conf_abelian[(position - shift[3] + shift[4]) * 4 + 3].r *
        sin(-conf_abelian[(position - shift[3] + shift[4]) * 4 + 3].phi -
            gauge_abelian[position - shift[3] + shift[4]]);
  }

  return abelian(sqrt(complex_tmp.real * complex_tmp.real +
                      complex_tmp.imag * complex_tmp.imag),
                 atan2(complex_tmp.imag, complex_tmp.real));
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

  complex_t complex_tmp(0, 0);

  // mu = 0
  if (x < x_size - 1) {
    complex_tmp.real +=
        cos(conf_angles[position * 4] - gauge_angles[position + shift[0]]);
    complex_tmp.imag +=
        sin(conf_angles[position * 4] - gauge_angles[position + shift[0]]);
  } else {
    complex_tmp.real += cos(conf_angles[position * 4] -
                            gauge_angles[position + shift[0] - shift[1]]);
    complex_tmp.imag += sin(conf_angles[position * 4] -
                            gauge_angles[position + shift[0] - shift[1]]);
  }

  if (x > 0) {
    complex_tmp.real += cos(-conf_angles[(position - shift[0]) * 4] -
                            gauge_angles[position - shift[0]]);
    complex_tmp.imag += sin(-conf_angles[(position - shift[0]) * 4] -
                            gauge_angles[position - shift[0]]);
  } else {
    complex_tmp.real += cos(-conf_angles[(position - shift[0] + shift[1]) * 4] -
                            gauge_angles[position - shift[0] + shift[1]]);
    complex_tmp.imag += sin(-conf_angles[(position - shift[0] + shift[1]) * 4] -
                            gauge_angles[position - shift[0] + shift[1]]);
  }

  // mu = 1

  if (y < y_size - 1) {
    complex_tmp.real +=
        cos(conf_angles[position * 4 + 1] - gauge_angles[position + shift[1]]);
    complex_tmp.imag +=
        sin(conf_angles[position * 4 + 1] - gauge_angles[position + shift[1]]);
  } else {
    complex_tmp.real += cos(conf_angles[position * 4 + 1] -
                            gauge_angles[position + shift[1] - shift[2]]);
    complex_tmp.imag += sin(conf_angles[position * 4 + 1] -
                            gauge_angles[position + shift[1] - shift[2]]);
  }

  if (y > 0) {
    complex_tmp.real += cos(-conf_angles[(position - shift[1]) * 4 + 1] -
                            gauge_angles[position - shift[1]]);
    complex_tmp.imag += sin(-conf_angles[(position - shift[1]) * 4 + 1] -
                            gauge_angles[position - shift[1]]);

  } else {
    complex_tmp.real +=
        cos(-conf_angles[(position - shift[1] + shift[2]) * 4 + 1] -
            gauge_angles[position - shift[1] + shift[2]]);
    complex_tmp.imag +=
        sin(-conf_angles[(position - shift[1] + shift[2]) * 4 + 1] -
            gauge_angles[position - shift[1] + shift[2]]);
  }

  // mu = 2

  if (z < z_size - 1) {
    complex_tmp.real +=
        cos(conf_angles[position * 4 + 2] - gauge_angles[position + shift[2]]);
    complex_tmp.imag +=
        sin(conf_angles[position * 4 + 2] - gauge_angles[position + shift[2]]);
  } else {
    complex_tmp.real += cos(conf_angles[position * 4 + 2] -
                            gauge_angles[position + shift[2] - shift[3]]);
    complex_tmp.imag += sin(conf_angles[position * 4 + 2] -
                            gauge_angles[position + shift[2] - shift[3]]);
  }

  if (z > 0) {
    complex_tmp.real += cos(-conf_angles[(position - shift[2]) * 4 + 2] -
                            gauge_angles[position - shift[2]]);
    complex_tmp.imag += sin(-conf_angles[(position - shift[2]) * 4 + 2] -
                            gauge_angles[position - shift[2]]);
  } else {
    complex_tmp.real +=
        cos(-conf_angles[(position - shift[2] + shift[3]) * 4 + 2] -
            gauge_angles[position - shift[2] + shift[3]]);
    complex_tmp.imag +=
        sin(-conf_angles[(position - shift[2] + shift[3]) * 4 + 2] -
            gauge_angles[position - shift[2] + shift[3]]);
  }

  // mu = 3

  if (t < t_size - 1) {
    complex_tmp.real +=
        cos(conf_angles[position * 4 + 3] - gauge_angles[position + shift[3]]);
    complex_tmp.imag +=
        sin(conf_angles[position * 4 + 3] - gauge_angles[position + shift[3]]);
  } else {
    complex_tmp.real += cos(conf_angles[position * 4 + 3] -
                            gauge_angles[position + shift[3] - shift[4]]);
    complex_tmp.imag += sin(conf_angles[position * 4 + 3] -
                            gauge_angles[position + shift[3] - shift[4]]);
  }

  if (t > 0) {
    complex_tmp.real += cos(-conf_angles[(position - shift[3]) * 4 + 3] -
                            gauge_angles[position - shift[3]]);
    complex_tmp.imag += sin(-conf_angles[(position - shift[3]) * 4 + 3] -
                            gauge_angles[position - shift[3]]);
  } else {
    complex_tmp.real +=
        cos(-conf_angles[(position - shift[3] + shift[4]) * 4 + 3] -
            gauge_angles[position - shift[3] + shift[4]]);
    complex_tmp.imag +=
        sin(-conf_angles[(position - shift[3] + shift[4]) * 4 + 3] -
            gauge_angles[position - shift[3] + shift[4]]);
  }

  return abelian(sqrt(complex_tmp.real * complex_tmp.real +
                      complex_tmp.imag * complex_tmp.imag),
                 atan2(complex_tmp.imag, complex_tmp.real));
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

void heat_bath(complex_t &g, complex_t &K, double temperature,
               double *random_numbers) {
  double temp_factor;
  long double temp_exponent;
  spin spin_new;

  // generation of new spin variable
  temp_factor = sqrt(K.real * K.real + K.imag * K.imag) / temperature;
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
    angle_new = acos(x) - atan2(K.imag, K.real);
  else
    angle_new = -acos(x) - atan2(K.imag, K.real);

  g = complex_t(cos(angle_new), sin(angle_new));
}

complex_t contribution_site(std::vector<complex_t> &gauge_complex,
                            std::vector<complex_t> &conf_complex, int x, int y,
                            int z, int t, int position,
                            std::vector<int> &shift) {

  complex_t complex_tmp(0, 0);

  // mu = 0
  if (x < x_size - 1) {
    complex_tmp = complex_tmp + (conf_complex[position * 4] ^
                                 gauge_complex[position + shift[0]]);
  } else {
    complex_tmp = complex_tmp + (conf_complex[position * 4] ^
                                 gauge_complex[position + shift[0] - shift[1]]);
  }

  if (x > 0) {
    complex_tmp = complex_tmp + (conf_complex[(position - shift[0]) * 4] &
                                 gauge_complex[position - shift[0]]);
  } else {
    complex_tmp =
        complex_tmp + (conf_complex[(position - shift[0] + shift[1]) * 4] &
                       gauge_complex[position - shift[0] + shift[1]]);
  }

  // mu = 1

  if (y < y_size - 1) {
    complex_tmp = complex_tmp + (conf_complex[position * 4 + 1] ^
                                 gauge_complex[position + shift[1]]);
  } else {
    complex_tmp = complex_tmp + (conf_complex[position * 4 + 1] ^
                                 gauge_complex[position + shift[1] - shift[2]]);
  }

  if (y > 0) {
    complex_tmp = complex_tmp + (conf_complex[(position - shift[1]) * 4 + 1] &
                                 gauge_complex[position - shift[1]]);
  } else {
    complex_tmp =
        complex_tmp + (conf_complex[(position - shift[1] + shift[2]) * 4 + 1] &
                       gauge_complex[position - shift[1] + shift[2]]);
  }

  // mu = 2

  if (z < z_size - 1) {
    complex_tmp = complex_tmp + (conf_complex[position * 4 + 2] ^
                                 gauge_complex[position + shift[2]]);
  } else {
    complex_tmp = complex_tmp + (conf_complex[position * 4 + 2] ^
                                 gauge_complex[position + shift[2] - shift[3]]);
  }

  if (z > 0) {
    complex_tmp = complex_tmp + (conf_complex[(position - shift[2]) * 4 + 2] &
                                 gauge_complex[position - shift[2]]);
  } else {
    complex_tmp =
        complex_tmp + (conf_complex[(position - shift[2] + shift[3]) * 4 + 2] &
                       gauge_complex[position - shift[2] + shift[3]]);
  }

  // mu = 3

  if (t < t_size - 1) {
    complex_tmp = complex_tmp + (conf_complex[position * 4 + 3] ^
                                 gauge_complex[position + shift[3]]);
  } else {
    complex_tmp = complex_tmp + (conf_complex[position * 4 + 3] ^
                                 gauge_complex[position + shift[3] - shift[4]]);
  }

  if (t > 0) {
    complex_tmp = complex_tmp + (conf_complex[(position - shift[3]) * 4 + 3] &
                                 gauge_complex[position - shift[3]]);
  } else {
    complex_tmp =
        complex_tmp + (conf_complex[(position - shift[3] + shift[4]) * 4 + 3] &
                       gauge_complex[position - shift[3] + shift[4]]);
  }

  return complex_tmp;
}

void heat_bath_update(std::vector<complex_t> &gauge_complex,
                      std::vector<complex_t> &conf_complex,
                      double temperature) {

  std::vector<double> random_numbers;

  // generate random numbers for heat bath
  random_numbers =
      generate_random_numbers1(2 * x_size * y_size * z_size * t_size);

  std::vector<int> shift = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  complex_t A;
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

void overrelaxation_update(std::vector<complex_t> &gauge_complex,
                           std::vector<complex_t> &conf_complex) {

  std::vector<int> shift = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  complex_t A;
  int position = 0;
  double A_phi;
  double gauge_phi;

  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {

          A = contribution_site(gauge_complex, conf_complex, x, y, z, t,
                                position, shift);

          A_phi = atan2(A.imag, A.real);
          gauge_phi = -atan2(gauge_complex[position].imag,
                             gauge_complex[position].real) -
                      2 * A_phi;

          gauge_complex[position].real = cos(gauge_phi);
          gauge_complex[position].imag = sin(gauge_phi);

          position++;
        }
      }
    }
  }
}

std::tuple<double, double>
relaxation_update(std::vector<complex_t> &gauge_complex,
                  std::vector<complex_t> &conf_complex) {

  std::vector<int> shift = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  complex_t A;
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

          gauge_phi =
              atan2(gauge_complex[position].imag, gauge_complex[position].real);
          A_phi = atan2(A.imag, A.real);

          gauge_complex[position].real = cos(A_phi);
          gauge_complex[position].imag = -sin(A_phi);

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

void make_simulated_annealing(std::vector<complex_t> &gauge_complex,
                              std::vector<complex_t> &conf_complex,
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

void make_maximization_final(std::vector<complex_t> &gauge_complex,
                             std::vector<complex_t> &conf_complex, int OR_steps,
                             double tolerance_maximal,
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

// int get_dirac_plaket_number(std::vector<complex_t> &gauge_complex,
//                             std::vector<complex_t> &conf_complex) {
//   int data_size = x_size * y_size * z_size * t_size;
//   std::vector<std::vector<int>> singular(6, std::vector<int>(data_size));
//   link1 link(x_size, y_size, z_size, t_size);

//   int number = 0;

//   SPACE_ITER_START
//   for (int mu = 0; mu < 4; mu++) {
//     link.move_dir(mu);
//     for (int nu = mu + 1; nu < 4; nu++) {

//       if (link.monopole_plaket_singular_mu(angles, nu) != 0)
//         number++;
//     }
//   }
//   SPACE_ITER_END
//   return number;
// }

void apply_gauge_Landau_complex(std::vector<complex_t> &gauge_complex,
                                std::vector<complex_t> &conf_complex) {
  link1 link(x_size, y_size, z_size, t_size);

  complex_t complex_tmp;
  double norm;

  SPACE_ITER_START

  for (int mu = 0; mu < 4; mu++) {

    complex_tmp = gauge_complex[link.place / 4] * conf_complex[link.place + mu];

    link.move(mu, 1);

    complex_tmp = complex_tmp ^ gauge_complex[link.place / 4];

    link.move(mu, -1);

    conf_complex[link.place + mu] = complex_tmp;
  }

  SPACE_ITER_END
}

void apply_gauge_Landau(std::vector<complex_t> &gauge_complex,
                        std::vector<su2> &conf_su2) {
  link1 link(x_size, y_size, z_size, t_size);

  su2 A;

  SPACE_ITER_START

  A = su2(gauge_complex[link.place / 4].real, 0, 0,
          gauge_complex[link.place / 4].imag);

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