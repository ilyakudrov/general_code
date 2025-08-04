#include "../include/mag.h"
#include "../include/data.h"
#include "../include/link.h"
#include "../include/matrix.h"

#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
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

std::vector<Eigen::Matrix3d> make_conf_contribution(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf_su2) {
  std::vector<Eigen::Matrix3d> conf_contribution(conf_su2.array.size());
  for (int i = 0; i < conf_su2.array.size(); i++) {
    double q1 = conf_su2[i].a1 * conf_su2[i].a1;
    double q2 = conf_su2[i].a2 * conf_su2[i].a2;
    double q3 = conf_su2[i].a3 * conf_su2[i].a3;
    double q12 = conf_su2[i].a1 * conf_su2[i].a2;
    double q13 = conf_su2[i].a1 * conf_su2[i].a3;
    double q23 = conf_su2[i].a2 * conf_su2[i].a3;
    double q01 = conf_su2[i].a0 * conf_su2[i].a1;
    double q02 = conf_su2[i].a0 * conf_su2[i].a2;
    double q03 = conf_su2[i].a0 * conf_su2[i].a3;
    double A5 = conf_su2[i].a0 * conf_su2[i].a0 - q1 - q2 - q3;
    conf_contribution[i] << A5 + 2 * q1, 2 * (q12 + q03), 2 * (q13 - q02),
        2 * (q12 - q03), A5 + 2 * q2, 2 * (q23 + q01), 2 * (q13 + q02),
        2 * (q23 - q01), A5 + 2 * q3;
  }
  return conf_contribution;
}

std::vector<spin> read_spins(std::string spins_path) {

  int data_size = x_size * y_size * z_size * t_size;

  std::vector<spin> spins;
  spins.reserve(x_size * y_size * z_size * t_size);

  std::ifstream stream(spins_path);
  std::vector<double> v(3 * data_size);

  if (!stream.read((char *)&v[0], 3 * data_size * sizeof(double)))
    std::cout << "read_spins error: " << spins_path << std::endl;

  for (int i = 0; i < data_size; i++) {
    spins.push_back(spin(v[3 * i], v[3 * i + 1], v[3 * i + 2]));
  }
  stream.close();

  return spins;
}

std::vector<spin> read_spins(std::string spins_path,
                             DataPatternLexicographical &data_pattern) {
  int data_size = data_pattern.get_lattice_size();
  std::vector<spin> spins;
  spins.reserve(data_size);
  std::ifstream stream(spins_path);
  std::vector<double> v(3 * data_size);
  if (!stream.read((char *)&v[0], 3 * data_size * sizeof(double)))
    std::cout << "read_spins error: " << spins_path << std::endl;
  for (int i = 0; i < data_size; i++) {
    spins.push_back(spin(v[3 * i], v[3 * i + 1], v[3 * i + 2]));
  }
  stream.close();
  return spins;
}

std::vector<Eigen::Vector3d>
read_vectors(std::string spins_path, DataPatternLexicographical &data_pattern) {
  int data_size = data_pattern.get_lattice_size();
  std::vector<Eigen::Vector3d> spins;
  spins.reserve(data_size);
  std::ifstream stream(spins_path);
  std::vector<double> v(3 * data_size);
  if (!stream.read((char *)&v[0], 3 * data_size * sizeof(double)))
    std::cout << "read_spins error: " << spins_path << std::endl;
  for (int i = 0; i < data_size; i++) {
    spins.push_back(Eigen::Vector3d(v[3 * i], v[3 * i + 1], v[3 * i + 2]));
  }
  stream.close();
  return spins;
}

void write_spins(std::string output_path, std::vector<spin> spins) {

  int data_size = x_size * y_size * z_size * t_size;

  std::ofstream stream(output_path);
  std::vector<double> v;
  v.reserve(3 * data_size);

  for (int i = 0; i < data_size; i++) {
    v.push_back(spins[i].a1);
    v.push_back(spins[i].a2);
    v.push_back(spins[i].a3);
  }

  if (!stream.write((char *)&v[0], 3 * data_size * sizeof(double)))
    std::cout << "write_spins error: " << output_path << std::endl;

  stream.close();
}

void write_spins(std::string output_path, std::vector<Eigen::Vector3d> spins,
                 DataPatternLexicographical &data_pattern) {
  int data_size = data_pattern.get_lattice_size();
  std::ofstream stream(output_path);
  std::vector<double> v;
  v.reserve(3 * data_size);
  for (int i = 0; i < data_size; i++) {
    v.push_back(spins[i][0]);
    v.push_back(spins[i][1]);
    v.push_back(spins[i][2]);
  }
  if (!stream.write((char *)&v[0], 3 * data_size * sizeof(double)))
    std::cout << "write_spins error: " << output_path << std::endl;
  stream.close();
}

std::vector<double> generate_random_numbers_sphere(int vector_size) {
  // using time to set a seed
  unsigned seed = time(NULL);
  //   unsigned seed = 123;
  // Since C rand() algorythm seems to have small cycle length it's important to
  // use more sophisticated method. Subtract_with_carry_engine is defined in
  // C++11. It realises subtract with carry lagged Fibonacci method. It seems to
  // meet the requirements
  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);
  std::vector<double> random_numbers;
  random_numbers.reserve(2 * vector_size);

  double x1, x2;

  for (int i = 0; i < vector_size; i++) {
    do {
      x1 = 2 * (double)random_generator() / random_generator.max() - 1;
      x2 = 2 * (double)random_generator() / random_generator.max() - 1;

      // to aviod errors due to rounding
    } while (x1 * x1 + x2 * x2 >= 1);

    random_numbers.push_back(x1);
    random_numbers.push_back(x2);
  }

  return random_numbers;
}

std::vector<spin> generate_spins_uniform() {
  int data_size = x_size * y_size * z_size * t_size;

  std::vector<double> random_numbers =
      generate_random_numbers_sphere(data_size);

  std::vector<spin> spins;
  spins.reserve(data_size);

  spin spin_tmp;
  double a;
  for (int i = 0; i < data_size; i++) {

    a = 1 - random_numbers[2 * i] * random_numbers[2 * i] -
        random_numbers[2 * i + 1] * random_numbers[2 * i + 1];

    spin_tmp =
        spin(2 * random_numbers[2 * i] * sqrt(a),
             2 * random_numbers[2 * i + 1] * sqrt(a),
             1 - 2 * (random_numbers[2 * i] * random_numbers[2 * i] +
                      random_numbers[2 * i + 1] * random_numbers[2 * i + 1]));

    spins.push_back(spin_tmp);
  }

  return spins;
}

std::vector<spin>
generate_spins_uniform(DataPatternLexicographical &data_pattern) {
  int data_size = data_pattern.get_lattice_size();
  std::vector<double> random_numbers =
      generate_random_numbers_sphere(data_size);
  std::vector<spin> spins;
  spins.reserve(data_size);
  spin spin_tmp;
  double a;
  for (int i = 0; i < data_size; i++) {
    a = 1 - random_numbers[2 * i] * random_numbers[2 * i] -
        random_numbers[2 * i + 1] * random_numbers[2 * i + 1];
    spin_tmp =
        spin(2 * random_numbers[2 * i] * sqrt(a),
             2 * random_numbers[2 * i + 1] * sqrt(a),
             1 - 2 * (random_numbers[2 * i] * random_numbers[2 * i] +
                      random_numbers[2 * i + 1] * random_numbers[2 * i + 1]));
    spins.push_back(spin_tmp);
  }
  return spins;
}

std::vector<Eigen::Vector3d>
generate_vectors_uniform(DataPatternLexicographical &data_pattern) {
  int data_size = data_pattern.get_lattice_size();
  std::vector<double> random_numbers =
      generate_random_numbers_sphere(data_size);
  std::vector<Eigen::Vector3d> spins;
  spins.reserve(data_size);
  Eigen::Vector3d spin_tmp;
  double a;
  for (int i = 0; i < data_size; i++) {
    a = 1 - random_numbers[2 * i] * random_numbers[2 * i] -
        random_numbers[2 * i + 1] * random_numbers[2 * i + 1];
    spin_tmp = Eigen::Vector3d(
        2 * random_numbers[2 * i] * sqrt(a),
        2 * random_numbers[2 * i + 1] * sqrt(a),
        1 - 2 * (random_numbers[2 * i] * random_numbers[2 * i] +
                 random_numbers[2 * i + 1] * random_numbers[2 * i + 1]));
    spins.push_back(spin_tmp);
  }
  return spins;
}

std::vector<double> generate_random_numbers(int vector_size) {
  // using time to set a seed
  unsigned seed = time(NULL);
  //   unsigned seed = 123;

  // Since C rand() algorythm seems to have small cycle length it's important to
  // use more sophisticated method. Subtract_with_carry_engine is defined in
  // C++11. It realises subtract with carry lagged Fibonacci method. It seems to
  // meet the requirements
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

void heat_bath(spin &spins, spin &neighbour, double temperature,
               double *random_numbers) {
  double temp_factor;
  long double temp_exponent;
  spin spin_new;

  // generation of new spin variable
  double neighbour_norm = neighbour.norm();
  temp_factor = neighbour_norm / temperature;
  double R;
  if (temp_factor < 700) {
    temp_exponent = exp(temp_factor);

    // generation of random variable with distribution p ~ exp(ax)
    spin_new.a3 = log(random_numbers[0] * (temp_exponent - 1 / temp_exponent) +
                      1 / temp_exponent) /
                  temp_factor;
  } else
    spin_new.a3 = log(random_numbers[0]) / temp_factor + 1;

  R = 1 - spin_new.a3 * spin_new.a3;

  spin_new.a1 = random_numbers[1] - 0.5;
  spin_new.a2 = random_numbers[2] - 0.5;

  R = sqrt(R / (spin_new.a1 * spin_new.a1 + spin_new.a2 * spin_new.a2));

  spin_new.a1 = spin_new.a1 * R;
  spin_new.a2 = spin_new.a2 * R;

  double neighbour3 = neighbour.a3 / neighbour_norm;

  double c, c1;
  double f1, f2;
  // there's no O31
  double O11, O12, O13, O21, O22, O23, O32, O33;

  c = sqrt(1 - neighbour3 * neighbour3);

  if (c == 0.) {
    // result of heat bath update
    spins = spin(spin_new.a1, spin_new.a2, spin_new.a3);

  } else {
    c1 = 1. / (neighbour_norm * c);
    f1 = neighbour.a1 * c1;
    f2 = neighbour.a2 * c1;
    O11 = f2;
    O12 = f1 * neighbour3;
    O13 = f1 * c;
    O21 = -f1;
    O22 = f2 * neighbour3;
    O23 = f2 * c;
    O32 = -c;
    O33 = neighbour3;

    // result of heat bath update
    spins = spin(O11 * spin_new.a1 + O12 * spin_new.a2 + O13 * spin_new.a3,
                 O21 * spin_new.a1 + O22 * spin_new.a2 + O23 * spin_new.a3,
                 O32 * spin_new.a2 + O33 * spin_new.a3);
  }
}

void heat_bath_test1(spin &spins, const spin &neighbour, double &temperature,
                     double *random_numbers) {
  double temp_factor;
  long double temp_exponent;
  spin spin_new;

  // generation of new spin variable
  double neighbour_norm = neighbour.norm();
  temp_factor = neighbour_norm / temperature;
  double R;
  // generation of random variable with distribution p ~ exp(ax)
  temp_exponent = exp(temp_factor);
  spin_new.a3 = log(random_numbers[0] * (temp_exponent - 1 / temp_exponent) +
                    1 / temp_exponent) /
                temp_factor;

  R = sqrt(1 - spin_new.a3 * spin_new.a3);
  double phi = 2 * M_PI * random_numbers[1];
  spin_new.a1 = R * cos(phi);
  spin_new.a2 = R * sin(phi);
  double neighbour3 = neighbour.a3 / neighbour_norm;

  double c, c1;
  double f1, f2;
  c = sqrt(1 - neighbour3 * neighbour3);
  c1 = 1. / (neighbour_norm * c);
  f1 = neighbour.a1 * c1;
  f2 = neighbour.a2 * c1;

  // result of heat bath update
  spins = spin(
      f2 * spin_new.a1 + f1 * neighbour3 * spin_new.a2 + f1 * c * spin_new.a3,
      -f1 * spin_new.a1 + f2 * neighbour3 * spin_new.a2 + f2 * c * spin_new.a3,
      -c * spin_new.a2 + neighbour3 * spin_new.a3);
}

void heat_bath_test2(
    spin &spins, const spin &neighbour, double &temperature,
    std::subtract_with_carry_engine<unsigned, 24, 10, 24> &random_generator) {
  double temp_factor;
  long double temp_exponent;
  spin spin_new;
  // generation of new spin variable
  double neighbour_norm = neighbour.norm();
  temp_factor = neighbour_norm / temperature;
  // generation of random variable with distribution p ~ exp(ax)
  temp_exponent = exp(temp_factor);
  spin_new.a3 =
      log(((double)(random_generator.min() + 1) / random_generator.max() +
           (double)random_generator() / (random_generator.max() + 2)) *
              (temp_exponent - 1 / temp_exponent) +
          1 / temp_exponent) /
      temp_factor;

  double R = sqrt(1 - spin_new.a3 * spin_new.a3);
  double phi =
      2 * M_PI * (double)random_generator() / (random_generator.max() + 1);
  spin_new.a1 = R * cos(phi);
  spin_new.a2 = R * sin(phi);
  double neighbour3 = neighbour.a3 / neighbour_norm;
  double c = sqrt(1 - neighbour3 * neighbour3);
  double c1 = 1. / (neighbour_norm * c);
  double f1 = neighbour.a1 * c1;
  double f2 = neighbour.a2 * c1;

  // result of heat bath update
  spins = spin(
      f2 * spin_new.a1 + f1 * neighbour3 * spin_new.a2 + f1 * c * spin_new.a3,
      -f1 * spin_new.a1 + f2 * neighbour3 * spin_new.a2 + f2 * c * spin_new.a3,
      -c * spin_new.a2 + neighbour3 * spin_new.a3);
}

void heat_bath_test3(
    spin &spins, const spin &neighbour, double &temperature,
    std::subtract_with_carry_engine<unsigned, 24, 10, 24> &random_generator) {
  double temp_factor;
  spin spin_new;
  // generation of new spin variable
  double neighbour_norm = neighbour.norm();
  temp_factor = neighbour_norm / temperature;
  // generation of random variable with distribution p ~ exp(ax)
  do {
    spin_new.a3 =
        1 + log((double)random_generator() / (random_generator.max() + 1)) /
                temp_factor;
  } while (spin_new.a3 < -1);

  double R = sqrt(1 - spin_new.a3 * spin_new.a3);
  double phi =
      2 * M_PI * (double)random_generator() / (random_generator.max() + 1);
  spin_new.a1 = R * cos(phi);
  spin_new.a2 = R * sin(phi);
  double neighbour3 = neighbour.a3 / neighbour_norm;
  double c = sqrt(1 - neighbour3 * neighbour3);
  double c1 = 1. / (neighbour_norm * c);
  double f1 = neighbour.a1 * c1;
  double f2 = neighbour.a2 * c1;

  // result of heat bath update
  spins = spin(
      f2 * spin_new.a1 + f1 * neighbour3 * spin_new.a2 + f1 * c * spin_new.a3,
      -f1 * spin_new.a1 + f2 * neighbour3 * spin_new.a2 + f2 * c * spin_new.a3,
      -c * spin_new.a2 + neighbour3 * spin_new.a3);
}

void heat_bath_test4(
    Eigen::Vector3d &spins, const Eigen::Vector3d &neighbour,
    double &temperature,
    std::subtract_with_carry_engine<unsigned, 24, 10, 24> &random_generator) {
  double temp_factor;
  long double temp_exponent;
  spin spin_new;

  // generation of new spin variable
  double neighbour_norm = neighbour.norm();
  temp_factor = neighbour_norm / temperature;
  // generation of random variable with distribution p ~ exp(ax)
  temp_exponent = exp(temp_factor);
  spin_new.a3 =
      log(((double)(random_generator.min() + 1) / random_generator.max() +
           (double)random_generator() / (random_generator.max() + 2)) *
              (temp_exponent - 1 / temp_exponent) +
          1 / temp_exponent) /
      temp_factor;

  double R = sqrt(1 - spin_new.a3 * spin_new.a3);

  double phi =
      2 * M_PI * (double)random_generator() / (random_generator.max() + 1);
  spin_new.a1 = R * cos(phi);
  spin_new.a2 = R * sin(phi);
  double neighbour3 = neighbour[2] / neighbour_norm;
  double c = sqrt(1 - neighbour3 * neighbour3);
  double c1 = 1. / (neighbour_norm * c);
  double f1 = neighbour[0] * c1;
  double f2 = neighbour[1] * c1;

  // result of heat bath update
  spins = Eigen::Vector3d(
      f2 * spin_new.a1 + f1 * neighbour3 * spin_new.a2 + f1 * c * spin_new.a3,
      -f1 * spin_new.a1 + f2 * neighbour3 * spin_new.a2 + f2 * c * spin_new.a3,
      -c * spin_new.a2 + neighbour3 * spin_new.a3);
}

void heat_bath_test5(
    Eigen::Vector3d &spins, const Eigen::Vector3d &neighbour,
    double &temperature,
    std::subtract_with_carry_engine<unsigned, 24, 10, 24> &random_generator) {
  double temp_factor;
  long double temp_exponent;
  Eigen::Vector3d spin_new;

  // generation of new spin variable
  double neighbour_norm = neighbour.norm();
  temp_factor = neighbour_norm / temperature;
  // generation of random variable with distribution p ~ exp(ax)
  temp_exponent = exp(temp_factor);
  do {
    spin_new[2] =
        1 + log((double)random_generator() / (random_generator.max() + 1)) /
                temp_factor;
  } while (spin_new[2] < -1);

  double R = sqrt(1 - spin_new[2] * spin_new[2]);
  double phi =
      2 * M_PI * (double)random_generator() / (random_generator.max() + 1);
  spin_new[0] = R * cos(phi);
  spin_new[1] = R * sin(phi);
  double neighbour3 = neighbour[2] / neighbour_norm;
  double c = sqrt(1 - neighbour3 * neighbour3);
  double c1 = 1. / (neighbour_norm * c);
  double f1 = neighbour[0] * c1;
  double f2 = neighbour[1] * c1;

  // result of heat bath update
  spins = Eigen::Vector3d(
      f2 * spin_new[0] + f1 * neighbour3 * spin_new[1] + f1 * c * spin_new[2],
      -f1 * spin_new[0] + f2 * neighbour3 * spin_new[1] + f2 * c * spin_new[2],
      -c * spin_new[1] + neighbour3 * spin_new[2]);
}

spin contribution_site(std::vector<spin> &spins,
                       const std::vector<su2> &conf_su2, int x, int y, int z,
                       int t, int position, std::vector<int> &shift) {
  spin A(0, 0, 0);
  // mu = 0
  if (x < x_size - 1) {
    A.contribution1(conf_su2[position * 4], spins[position + shift[0]]);
  } else
    A.contribution1(conf_su2[position * 4],
                    spins[position + shift[0] - shift[1]]);
  if (x > 0)
    A.contribution1_conj(conf_su2[(position - shift[0]) * 4],
                         spins[position - shift[0]]);
  else
    A.contribution1_conj(conf_su2[(position - shift[0] + shift[1]) * 4],
                         spins[position - shift[0] + shift[1]]);

  // mu = 1
  if (y < y_size - 1)
    A.contribution1(conf_su2[position * 4 + 1], spins[position + shift[1]]);
  else
    A.contribution1(conf_su2[position * 4 + 1],
                    spins[position + shift[1] - shift[2]]);
  if (y > 0)
    A.contribution1_conj(conf_su2[(position - shift[1]) * 4 + 1],
                         spins[position - shift[1]]);
  else
    A.contribution1_conj(conf_su2[(position - shift[1] + shift[2]) * 4 + 1],
                         spins[position - shift[1] + shift[2]]);

  // mu = 2
  if (z < z_size - 1)
    A.contribution1(conf_su2[position * 4 + 2], spins[position + shift[2]]);
  else
    A.contribution1(conf_su2[position * 4 + 2],
                    spins[position + shift[2] - shift[3]]);
  if (z > 0)
    A.contribution1_conj(conf_su2[(position - shift[2]) * 4 + 2],
                         spins[position - shift[2]]);
  else
    A.contribution1_conj(conf_su2[(position - shift[2] + shift[3]) * 4 + 2],
                         spins[position - shift[2] + shift[3]]);

  // mu = 3
  if (t < t_size - 1)
    A.contribution1(conf_su2[position * 4 + 3], spins[position + shift[3]]);
  else
    A.contribution1(conf_su2[position * 4 + 3],
                    spins[position + shift[3] - shift[4]]);
  if (t > 0)
    A.contribution1_conj(conf_su2[(position - shift[3]) * 4 + 3],
                         spins[position - shift[3]]);
  else
    A.contribution1_conj(conf_su2[(position - shift[3] + shift[4]) * 4 + 3],
                         spins[position - shift[3] + shift[4]]);
  return A;
}

spin contribution_site(std::vector<spin> &spins,
                       const std::vector<su2> &conf_su2,
                       DataPatternLexicographical &data_pattern, int x, int y,
                       int z, int t, int position, std::vector<int> &shift) {
  spin A(0, 0, 0);
  // mu = 0
  if (x < data_pattern.lat_dim[0] - 1) {
    A.contribution1(conf_su2[position * 4], spins[position + shift[0]]);
  } else
    A.contribution1(conf_su2[position * 4],
                    spins[position + shift[0] - shift[1]]);
  if (x > 0)
    A.contribution1_conj(conf_su2[(position - shift[0]) * 4],
                         spins[position - shift[0]]);
  else
    A.contribution1_conj(conf_su2[(position - shift[0] + shift[1]) * 4],
                         spins[position - shift[0] + shift[1]]);

  // mu = 1
  if (y < data_pattern.lat_dim[1] - 1)
    A.contribution1(conf_su2[position * 4 + 1], spins[position + shift[1]]);
  else
    A.contribution1(conf_su2[position * 4 + 1],
                    spins[position + shift[1] - shift[2]]);
  if (y > 0)
    A.contribution1_conj(conf_su2[(position - shift[1]) * 4 + 1],
                         spins[position - shift[1]]);
  else
    A.contribution1_conj(conf_su2[(position - shift[1] + shift[2]) * 4 + 1],
                         spins[position - shift[1] + shift[2]]);

  // mu = 2
  if (z < data_pattern.lat_dim[2] - 1)
    A.contribution1(conf_su2[position * 4 + 2], spins[position + shift[2]]);
  else
    A.contribution1(conf_su2[position * 4 + 2],
                    spins[position + shift[2] - shift[3]]);
  if (z > 0)
    A.contribution1_conj(conf_su2[(position - shift[2]) * 4 + 2],
                         spins[position - shift[2]]);
  else
    A.contribution1_conj(conf_su2[(position - shift[2] + shift[3]) * 4 + 2],
                         spins[position - shift[2] + shift[3]]);

  // mu = 3
  if (t < data_pattern.lat_dim[3] - 1)
    A.contribution1(conf_su2[position * 4 + 3], spins[position + shift[3]]);
  else
    A.contribution1(conf_su2[position * 4 + 3],
                    spins[position + shift[3] - shift[4]]);
  if (t > 0)
    A.contribution1_conj(conf_su2[(position - shift[3]) * 4 + 3],
                         spins[position - shift[3]]);
  else
    A.contribution1_conj(conf_su2[(position - shift[3] + shift[4]) * 4 + 3],
                         spins[position - shift[3] + shift[4]]);
  return A;
}

Eigen::Vector3d
contribution_site(std::vector<Eigen::Vector3d> &spins,
                  const std::vector<Eigen::Matrix3d> &conf_contribution,
                  DataPatternLexicographical &data_pattern, int x, int y, int z,
                  int t, int position, std::vector<int> &shift) {
  Eigen::Vector3d A(0, 0, 0);
  // mu = 0
  if (x < data_pattern.lat_dim[0] - 1) {
    A += conf_contribution[position * 4] * spins[position + shift[0]];
  } else
    A +=
        conf_contribution[position * 4] * spins[position + shift[0] - shift[1]];
  if (x > 0)
    A += conf_contribution[(position - shift[0]) * 4].transpose() *
         spins[position - shift[0]];
  else
    A += conf_contribution[(position - shift[0] + shift[1]) * 4].transpose() *
         spins[position - shift[0] + shift[1]];

  // mu = 1
  if (y < data_pattern.lat_dim[1] - 1)
    A += conf_contribution[position * 4 + 1] * spins[position + shift[1]];
  else
    A += conf_contribution[position * 4 + 1] *
         spins[position + shift[1] - shift[2]];
  if (y > 0)
    A += conf_contribution[(position - shift[1]) * 4 + 1].transpose() *
         spins[position - shift[1]];
  else
    A += conf_contribution[(position - shift[1] + shift[2]) * 4 + 1]
             .transpose() *
         spins[position - shift[1] + shift[2]];

  // mu = 2
  if (z < data_pattern.lat_dim[2] - 1)
    A += conf_contribution[position * 4 + 2] * spins[position + shift[2]];
  else
    A += conf_contribution[position * 4 + 2] *
         spins[position + shift[2] - shift[3]];
  if (z > 0)
    A += conf_contribution[(position - shift[2]) * 4 + 2].transpose() *
         spins[position - shift[2]];
  else
    A += conf_contribution[(position - shift[2] + shift[3]) * 4 + 2]
             .transpose() *
         spins[position - shift[2] + shift[3]];

  // mu = 3
  if (t < data_pattern.lat_dim[3] - 1)
    A += conf_contribution[position * 4 + 3] * spins[position + shift[3]];
  else
    A += conf_contribution[position * 4 + 3] *
         spins[position + shift[3] - shift[4]];
  if (t > 0)
    A += conf_contribution[(position - shift[3]) * 4 + 3].transpose() *
         spins[position - shift[3]];
  else
    A += conf_contribution[(position - shift[3] + shift[4]) * 4 + 3]
             .transpose() *
         spins[position - shift[3] + shift[4]];
  return A;
}

void heat_bath_update(std::vector<spin> &spins, std::vector<su2> &conf_su2,
                      double temperature) {
  std::vector<double> random_numbers;
  double time_tmp;
  double contribution_time = 0;
  double heat_bath_time = 0;

  // generate random numbers for heat bath
  random_numbers =
      generate_random_numbers(x_size * y_size * z_size * t_size * 3);

  //   std::subtract_with_carry_engine<unsigned, 24, 10, 24>
  //   random_generator(seed);

  std::vector<int> shift = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  spin A(0, 0, 0);
  int position = 0;
  int count = 0;

  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
          time_tmp = omp_get_wtime();
          A = contribution_site(spins, conf_su2, x, y, z, t, position, shift);
          contribution_time += omp_get_wtime() - time_tmp;
          time_tmp = omp_get_wtime();
          heat_bath(spins[position], A, temperature,
                    &random_numbers[count * 3]);
          position++;
          count++;
          heat_bath_time += omp_get_wtime() - time_tmp;
        }
      }
    }
  }
  std::cout << "contribution_time: " << contribution_time << std::endl;
  std::cout << "heat_bath_time: " << heat_bath_time << std::endl;
  normalize_spin(spins);
}

spin contribution_site(
    const std::vector<spin> &spins,
    const Data::LatticeData<DataPatternLexicographical, su2> &conf_su2,
    DataPatternLexicographical &data_pattern) {
  spin A(0, 0, 0);
  int index1 = data_pattern.get_index_site() * 4;
  int index2;
  for (int mu = 0; mu < 4; mu++) {
    data_pattern.move_forward(1, mu);
    index2 = data_pattern.get_index_site();
    data_pattern.move_backward(2, mu);
    A.contribution1(conf_su2[index1 + mu], spins[index2]);
    index2 = data_pattern.get_index_site();
    A.contribution1_conj(conf_su2[index2 * 4 + mu], spins[index2]);
    data_pattern.move_forward(1, mu);
  }
  return A;
}

void heat_bath_update(
    std::vector<spin> &spins,
    const Data::LatticeData<DataPatternLexicographical, su2> &conf_su2,
    double temperature) {
  std::vector<int> shift = {
      1, conf_su2.lat_dim[0], conf_su2.lat_dim[0] * conf_su2.lat_dim[1],
      conf_su2.lat_dim[0] * conf_su2.lat_dim[1] * conf_su2.lat_dim[2],
      conf_su2.lat_dim[0] * conf_su2.lat_dim[1] * conf_su2.lat_dim[2] *
          conf_su2.lat_dim[3]};
  DataPatternLexicographical data_pattern(conf_su2.lat_dim);
  unsigned seed = time(NULL);
  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);
  spin A;
  int position;
#pragma omp parallel for collapse(4) private(A, position)                      \
    firstprivate(data_pattern, temperature, random_generator, shift)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          position = data_pattern.get_index_site();
          A = contribution_site(spins, conf_su2.array, data_pattern, x, y, z, t,
                                position, shift);
          heat_bath_test2(spins[position], A, temperature, random_generator);
        }
      }
    }
  }
}

void heat_bath_update_fast(
    std::vector<spin> &spins,
    const Data::LatticeData<DataPatternLexicographical, su2> &conf_su2,
    double temperature) {
  std::vector<int> shift = {
      1, conf_su2.lat_dim[0], conf_su2.lat_dim[0] * conf_su2.lat_dim[1],
      conf_su2.lat_dim[0] * conf_su2.lat_dim[1] * conf_su2.lat_dim[2],
      conf_su2.lat_dim[0] * conf_su2.lat_dim[1] * conf_su2.lat_dim[2] *
          conf_su2.lat_dim[3]};
  int position;
  DataPatternLexicographical data_pattern(conf_su2.lat_dim);
  unsigned seed = time(NULL);
  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);
  spin A;
#pragma omp parallel for collapse(4) private(A, position)                      \
    firstprivate(data_pattern, temperature, random_generator, shift)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          position = data_pattern.get_index_site();
          A = contribution_site(spins, conf_su2.array, data_pattern, x, y, z, t,
                                position, shift);
          heat_bath_test3(spins[position], A, temperature, random_generator);
        }
      }
    }
  }
}

Eigen::Vector3d
contribution_site(const std::vector<Eigen::Vector3d> &spins,
                  const std::vector<Eigen::Matrix3d> &conf_contribution,
                  DataPatternLexicographical &data_pattern) {
  Eigen::Vector3d A(0, 0, 0);
  int index1 = data_pattern.get_index_site() * 4;
  int index2;
  for (int mu = 0; mu < 4; mu++) {
    data_pattern.move_forward(1, mu);
    index2 = data_pattern.get_index_site();
    data_pattern.move_backward(2, mu);
    A += conf_contribution[index1 + mu] * spins[index2];
    index2 = data_pattern.get_index_site();
    A += conf_contribution[index2 * 4 + mu].transpose() * spins[index2];
    data_pattern.move_forward(1, mu);
  }
  return A;
}

void heat_bath_update(std::vector<Eigen::Vector3d> &spins,
                      const std::vector<Eigen::Matrix3d> &conf_contribution,
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
  Eigen::Vector3d A;
  int position;
#pragma omp parallel for collapse(4) private(A, position)                      \
    firstprivate(data_pattern, temperature, random_generator, shift)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          position = data_pattern.get_index_site();
          A = contribution_site(spins, conf_contribution, data_pattern, x, y, z,
                                t, position, shift);
          heat_bath_test4(spins[position], A, temperature, random_generator);
        }
      }
    }
  }
}

void heat_bath_update_fast(
    std::vector<Eigen::Vector3d> &spins,
    const std::vector<Eigen::Matrix3d> &conf_contribution,
    DataPatternLexicographical &data_pattern, double temperature) {
  std::vector<int> shift = {1, data_pattern.lat_dim[0],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2] *
                                data_pattern.lat_dim[3]};
  unsigned seed = time(NULL);
  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);
  Eigen::Vector3d A;
  int position;
#pragma omp parallel for collapse(4) private(A, position)                      \
    firstprivate(data_pattern, temperature, random_generator, shift)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          position = data_pattern.get_index_site();
          A = contribution_site(spins, conf_contribution, data_pattern, x, y, z,
                                t, position, shift);
          heat_bath_test5(spins[position], A, temperature, random_generator);
        }
      }
    }
  }
}

void normalize_spin(std::vector<spin> &spins) {
  for (int i = 0; i < spins.size(); i++) {
    spins[i].normalize();
  }
}

double MAG_functional_su2_spin(std::vector<su2> &conf_su2,
                               std::vector<spin> &spins) {
  link1 link(x_size, y_size, z_size, t_size);
  double result = 0;

  SPACE_ITER_START

  for (int mu = 0; mu < 4; mu++) {

    link.move_dir(mu);
    result += link.get_consecutive_spin(spins, mu)->contribution(
                  *link.get_matrix(conf_su2)) *
              *link.get_spin(spins);
  }

  SPACE_ITER_END

  return result / (x_size * y_size * z_size * t_size * 4);
}

double MAG_functional_su2_spin(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf_su2,
    std::vector<spin> &spins) {
  DataPatternLexicographical data_pattern(conf_su2.lat_dim);
  double result = 0;
  int index_spin;
  int index_su2;
#pragma omp parallel for collapse(4) private(index_spin, index_su2)            \
    firstprivate(data_pattern) reduction(+ : result)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index_spin = data_pattern.get_index_site();
          for (int mu = 0; mu < 4; mu++) {
            index_su2 = data_pattern.get_index_link(mu);
            data_pattern.move_forward(1, mu);
            result += spins[data_pattern.get_index_site()].contribution(
                          conf_su2[index_su2]) *
                      spins[index_spin];
            data_pattern.move_backward(1, mu);
          }
        }
      }
    }
  }
  return result / (data_pattern.get_data_size());
}

double
MAG_functional_su2_spin(const std::vector<Eigen::Matrix3d> &conf_contribution,
                        std::vector<Eigen::Vector3d> &spins,
                        DataPatternLexicographical &data_pattern) {
  double result = 0;
  int index_spin;
  int index_su2;
#pragma omp parallel for collapse(4) private(index_spin, index_su2)            \
    firstprivate(data_pattern) reduction(+ : result)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index_spin = data_pattern.get_index_site();
          for (int mu = 0; mu < 4; mu++) {
            index_su2 = data_pattern.get_index_link(mu);
            data_pattern.move_forward(1, mu);
            result += spins[index_spin].transpose() *
                      conf_contribution[index_su2] *
                      spins[data_pattern.get_index_site()];
            data_pattern.move_backward(1, mu);
          }
        }
      }
    }
  }
  return result / (data_pattern.get_data_size());
}

double MAG_functional_su2(const std::vector<su2> &array) {
  double result = 0;
  for (int i = 0; i < x_size * y_size * z_size * t_size * 4; i++) {
    result += (array[i].sigma3_mult() * array[i].conj()).tr();
  }
  return result / (x_size * y_size * z_size * t_size * 4);
}

double MAG_functional_su2(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf_su2) {
  DataPatternLexicographical data_pattern(conf_su2.lat_dim);
  double result = 0;
  for (int i = 0; i < data_pattern.get_data_size(); i++) {
    result += (conf_su2[i].sigma3_mult() * conf_su2[i].conj()).tr();
  }
  return result / (data_pattern.get_data_size());
}

std::vector<su2> make_gauge(std::vector<spin> &spins) {
  std::vector<su2> gauge(x_size * y_size * z_size * t_size);
  for (int i = 0; i < spins.size(); i++) {
    gauge[i] = spins[i].GetGaugeMatrix();
  }
  return gauge;
}

std::vector<su2> gauge_tranformation(std::vector<su2> &conf_su2,
                                     std::vector<su2> &gauge) {
  std::vector<su2> conf_new(x_size * y_size * z_size * t_size * 4);
  link1 link(x_size, y_size, z_size, t_size);
  su2 A;
  SPACE_ITER_START
  for (int mu = 0; mu < 4; mu++) {
    link.move_dir(mu);
    A = gauge[link.place / 4] * *link.get_matrix(conf_su2);
    link.move(mu, 1);
    A = A * gauge[link.place / 4].conj();
    link.move(mu, -1);
    conf_new[link.place + mu] = A;
  }
  SPACE_ITER_END
  return conf_new;
}

void gauge_tranformation_spins(std::vector<su2> &conf_su2,
                               std::vector<spin> &spins) {
  link1 link(x_size, y_size, z_size, t_size);
  su2 A;
  SPACE_ITER_START
  for (int mu = 0; mu < 4; mu++) {
    link.move_dir(mu);
    A = spins[link.place / 4].GetGaugeMatrix() * *link.get_matrix(conf_su2);
    link.move(mu, 1);
    A = A * spins[link.place / 4].GetGaugeMatrix().conj();
    link.move(mu, -1);
    conf_su2[link.place + mu] = A;
  }
  SPACE_ITER_END
}

void gauge_tranformation_spins(
    Data::LatticeData<DataPatternLexicographical, su2> &conf_su2,
    std::vector<spin> &spins) {
  DataPatternLexicographical data_pattern(conf_su2.lat_dim);
  su2 A;
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          for (int mu = 0; mu < 4; mu++) {
            A = spins[data_pattern.get_index_site()].GetGaugeMatrix() *
                conf_su2[data_pattern.get_index_link(mu)];
            data_pattern.move_forward(1, mu);
            A = A *
                spins[data_pattern.get_index_site()].GetGaugeMatrix().conj();
            data_pattern.move_backward(1, mu);
            conf_su2[data_pattern.get_index_link(mu)] = A;
          }
        }
      }
    }
  }
}

void gauge_tranformation_spins(
    Data::LatticeData<DataPatternLexicographical, su2> &conf_su2,
    std::vector<Eigen::Vector3d> &spins) {
  DataPatternLexicographical data_pattern(conf_su2.lat_dim);
  su2 A;
  spin s;
  int index;
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          for (int mu = 0; mu < 4; mu++) {
            index = data_pattern.get_index_site();
            s.a1 = spins[index][0];
            s.a2 = spins[index][1];
            s.a3 = spins[index][2];
            A = s.GetGaugeMatrix() * conf_su2[data_pattern.get_index_link(mu)];
            data_pattern.move_forward(1, mu);
            index = data_pattern.get_index_site();
            s.a1 = spins[index][0];
            s.a2 = spins[index][1];
            s.a3 = spins[index][2];
            A = A * s.GetGaugeMatrix().conj();
            data_pattern.move_backward(1, mu);
            conf_su2[data_pattern.get_index_link(mu)] = A;
          }
        }
      }
    }
  }
}

std::vector<int> make_indices_qube(int qube_size) {
  link1 link(x_size, y_size, z_size, t_size);

  std::vector<int> indices;
  indices.reserve(x_size * y_size * z_size * t_size * 9);

  for (int t = 0; t < t_size; t += qube_size) {
    for (int z = 0; z < z_size; z += qube_size) {
      for (int y = 0; y < y_size; y += qube_size) {
        for (int x = 0; x < x_size; x += qube_size) {

          for (int t1 = t; t1 < t + qube_size; t1++) {
            for (int z1 = z; z1 < z + qube_size; z1++) {
              for (int y1 = y; y1 < y + qube_size; y1++) {
                for (int x1 = x; x1 < x + qube_size; x1++) {
                  link.go_update(x1, y1, z1, t1);

                  for (int mu = 0; mu < 4; mu++) {
                    link.move(mu, 1);

                    indices.push_back(link.place / 4);

                    link.move(mu, -1);
                  }

                  for (int mu = 0; mu < 4; mu++) {
                    link.move(mu, -1);

                    indices.push_back(link.place / 4);

                    link.move(mu, 1);
                  }

                  indices.push_back(link.place / 4);
                }
              }
            }
          }
        }
      }
    }
  }

  return indices;
}

spin contribution_site2(std::vector<spin> &spins, std::vector<su2> &conf_su2,
                        std::vector<int> &indices, int index) {

  spin A(0, 0, 0);

  int index_su2 = indices[index + 8];

  for (int mu = 0; mu < 4; mu++) {
    A.contribution1(conf_su2[index_su2 * 4 + mu], spins[indices[index + mu]]);
  }

  for (int mu = 0; mu < 4; mu++) {
    A.contribution1_conj(conf_su2[indices[index + 4 + mu] * 4 + mu],
                         spins[indices[index + 4 + mu]]);
  }

  return A;
}

void heat_bath_update_tets2(std::vector<spin> &spins,
                            std::vector<su2> &conf_su2,
                            std::vector<int> indices, double temperature) {

  std::vector<double> random_numbers;

  // generate random numbers for heat bath
  random_numbers =
      generate_random_numbers(x_size * y_size * z_size * t_size * 3);

  //   std::subtract_with_carry_engine<unsigned, 24, 10, 24>
  //   random_generator(seed);

  spin A(0, 0, 0);
  int count = 0;

  for (int index = 0; index < x_size * y_size * z_size * t_size; index++) {

    A = contribution_site2(spins, conf_su2, indices, index);

    heat_bath(spins[indices[index * 9 + 8]], A, temperature,
              &random_numbers[count * 3]);

    count++;
  }

  normalize_spin(spins);
}

void heat_bath_update_tets3(std::vector<spin> &spins,
                            std::vector<su2> &conf_su2, double temperature) {

  std::vector<double> random_numbers;

  // generate random numbers for heat bath
  random_numbers =
      generate_random_numbers(x_size * y_size * z_size * t_size * 3);

  //   std::subtract_with_carry_engine<unsigned, 24, 10, 24>
  //   random_generator(seed);

  std::vector<int> shift = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  spin A(0, 0, 0);
  int count = 0;
  int position;

  int qube_size = 2;

  for (int t = 0; t < t_size; t += qube_size) {
    for (int z = 0; z < z_size; z += qube_size) {
      for (int y = 0; y < y_size; y += qube_size) {
        for (int x = 0; x < x_size; x += qube_size) {

          for (int t1 = t; t1 < t + qube_size; t1++) {
            for (int z1 = z; z1 < z + qube_size; z1++) {
              for (int y1 = y; y1 < y + qube_size; y1++) {
                for (int x1 = x; x1 < x + qube_size; x1++) {

                  position = t1 * x_size * y_size * z_size +
                             z1 * x_size * y_size + y1 * x_size + x1;

                  A = contribution_site(spins, conf_su2, x1, y1, z1, t1,
                                        position, shift);

                  heat_bath(spins[position], A, temperature,
                            &random_numbers[count * 3]);

                  count++;
                }
              }
            }
          }
        }
      }
    }
  }

  normalize_spin(spins);
}

void make_indices_qube1(std::vector<int> &indices,
                        std::vector<char> &coordinates, int qube_size) {
  link1 link(x_size, y_size, z_size, t_size);

  for (int t = 0; t < t_size; t += qube_size) {
    for (int z = 0; z < z_size; z += qube_size) {
      for (int y = 0; y < y_size; y += qube_size) {
        for (int x = 0; x < x_size; x += qube_size) {

          for (int t1 = t; t1 < t + qube_size; t1++) {
            for (int z1 = z; z1 < z + qube_size; z1++) {
              for (int y1 = y; y1 < y + qube_size; y1++) {
                for (int x1 = x; x1 < x + qube_size; x1++) {
                  link.go_update(x1, y1, z1, t1);

                  indices.push_back(link.place / 4);

                  coordinates.push_back(x1);
                  coordinates.push_back(y1);
                  coordinates.push_back(z1);
                  coordinates.push_back(t1);
                }
              }
            }
          }
        }
      }
    }
  }
}

spin contribution_site4(std::vector<spin> &spins, std::vector<su2> &conf_su2,
                        char *coordinate, int position, std::vector<int> &shift,
                        std::vector<char> &lattice_size) {

  spin A(0, 0, 0);

  for (int mu = 0; mu < 4; mu++) {
    if (coordinate[mu] < lattice_size[mu] - 1)
      A.contribution1(conf_su2[position * 4 + mu], spins[position + shift[mu]]);
    else
      A.contribution1(conf_su2[position * 4 + mu],
                      spins[position + shift[mu] - shift[mu + 1]]);

    if (coordinate[mu] > 0)
      A.contribution1_conj(conf_su2[(position - shift[mu]) * 4 + mu],
                           spins[position - shift[mu]]);
    else
      A.contribution1_conj(
          conf_su2[(position - shift[mu] + shift[mu + 1]) * 4 + mu],
          spins[position - shift[mu] + shift[mu + 1]]);
  }

  return A;
}

void heat_bath_update_tets4(std::vector<spin> &spins,
                            std::vector<su2> &conf_su2,
                            std::vector<int> &indices,
                            std::vector<char> &coordinates,
                            double temperature) {

  std::vector<double> random_numbers;

  // generate random numbers for heat bath
  random_numbers =
      generate_random_numbers(x_size * y_size * z_size * t_size * 3);

  //   std::subtract_with_carry_engine<unsigned, 24, 10, 24>
  //   random_generator(seed);

  std::vector<int> shift = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  std::vector<char> lattice_size = {(char)x_size, (char)y_size, (char)z_size,
                                    (char)t_size};

  spin A(0, 0, 0);
  int count = 0;

  for (int index = 0; index < x_size * y_size * z_size * t_size; index++) {

    A = contribution_site4(spins, conf_su2, &coordinates[index * 4],
                           indices[index], shift, lattice_size);

    heat_bath(spins[indices[index]], A, temperature,
              &random_numbers[count * 3]);

    count++;
  }

  normalize_spin(spins);
}

void overrelaxation_update(std::vector<spin> &spins,
                           std::vector<su2> &conf_su2) {

  std::vector<int> shift = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  spin A(0, 0, 0);
  int position = 0;
  int count = 0;

  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {

          A = contribution_site(spins, conf_su2, x, y, z, t, position, shift);

          spins[position].reflect_fast(A);

          position++;
          count++;
        }
      }
    }
  }
  normalize_spin(spins);
}

void overrelaxation_update(
    std::vector<spin> &spins,
    const Data::LatticeData<DataPatternLexicographical, su2> &conf_su2) {
  std::vector<int> shift = {
      1, conf_su2.lat_dim[0], conf_su2.lat_dim[0] * conf_su2.lat_dim[1],
      conf_su2.lat_dim[0] * conf_su2.lat_dim[1] * conf_su2.lat_dim[2],
      conf_su2.lat_dim[0] * conf_su2.lat_dim[1] * conf_su2.lat_dim[2] *
          conf_su2.lat_dim[3]};
  int position;
  DataPatternLexicographical data_pattern(conf_su2.lat_dim);
  spin A;
#pragma omp parallel for collapse(4) private(A, position)                      \
    firstprivate(data_pattern, shift)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          position = data_pattern.get_index_site();
          A = contribution_site(spins, conf_su2.array, data_pattern, x, y, z, t,
                                position, shift);
          spins[position].reflect_fast(A);
        }
      }
    }
  }
}

void reflect(Eigen::Vector3d &s, Eigen::Vector3d &V) {
  double tmp1 = (V[0] * V[0] - V[1] * V[1] - V[2] * V[2]) * s[0] +
                2 * V[0] * (V[1] * s[1] + V[2] * s[2]);
  double tmp2 = (-V[0] * V[0] + V[1] * V[1] - V[2] * V[2]) * s[1] +
                2 * V[1] * (V[0] * s[0] + V[2] * s[2]);
  double tmp3 = (-V[0] * V[0] - V[1] * V[1] + V[2] * V[2]) * s[2] +
                2 * V[2] * (V[0] * s[0] + V[1] * s[1]);
  double norm2 = V.norm() * V.norm();
  s[0] = tmp1 / norm2;
  s[1] = tmp2 / norm2;
  s[2] = tmp3 / norm2;
}

void overrelaxation_update(
    std::vector<Eigen::Vector3d> &spins,
    const std::vector<Eigen::Matrix3d> &conf_contribution,
    DataPatternLexicographical &data_pattern) {
  std::vector<int> shift = {1, data_pattern.lat_dim[0],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2] *
                                data_pattern.lat_dim[3]};
  Eigen::Vector3d A;
  int position;
#pragma omp parallel for collapse(4) private(A, position)                      \
    firstprivate(data_pattern, shift)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          position = data_pattern.get_index_site();
          A = contribution_site(spins, conf_contribution, data_pattern, x, y, z,
                                t, position, shift);
          reflect(spins[position], A);
        }
      }
    }
  }
}

std::tuple<double, double> relaxation_update(std::vector<spin> &spins,
                                             std::vector<su2> &conf_su2) {
  std::vector<int> shift = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};
  spin A(0, 0, 0);
  int position = 0;
  int count = 0;
  double diff_average = 0;
  double diff_maximal = 0;
  double diff_tmp;
  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
          A = contribution_site(spins, conf_su2, x, y, z, t, position, shift);
          diff_tmp = spins[position].parallel(A);
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
  normalize_spin(spins);
  return std::tuple<double, double>(diff_maximal, diff_average);
}

std::tuple<double, double> relaxation_update(
    std::vector<spin> &spins,
    const Data::LatticeData<DataPatternLexicographical, su2> &conf_su2) {
  std::vector<int> shift = {
      1, conf_su2.lat_dim[0], conf_su2.lat_dim[0] * conf_su2.lat_dim[1],
      conf_su2.lat_dim[0] * conf_su2.lat_dim[1] * conf_su2.lat_dim[2],
      conf_su2.lat_dim[0] * conf_su2.lat_dim[1] * conf_su2.lat_dim[2] *
          conf_su2.lat_dim[3]};
  int position;
  DataPatternLexicographical data_pattern(conf_su2.lat_dim);
  spin A(0, 0, 0);
  double diff_average = 0;
  double diff_maximal = 0;
  double diff_tmp;
#pragma omp parallel for collapse(4) private(A, position, diff_tmp)            \
    firstprivate(data_pattern, shift) reduction(+ : diff_average)              \
    reduction(max : diff_maximal)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          position = data_pattern.get_index_site();
          A = contribution_site(spins, conf_su2.array, data_pattern, x, y, z, t,
                                position, shift);
          diff_tmp = spins[position].parallel(A);
          if (diff_tmp > diff_maximal)
            diff_maximal = diff_tmp;
          diff_average += diff_tmp;
        }
      }
    }
  }
  diff_average = diff_average / data_pattern.get_lattice_size();
  return std::tuple<double, double>(diff_maximal, diff_average);
}

double parallel(Eigen::Vector3d &s, const Eigen::Vector3d &V) {
  double vNorm = V.norm();
  double b1 = V[0] / vNorm;
  double b2 = V[1] / vNorm;
  double b3 = V[2] / vNorm;
  double c = b1 * s[0] + b2 * s[1] + b3 * s[2];
  s[0] = b1;
  s[1] = b2;
  s[2] = b3;
  return 1 - c;
}

std::tuple<double, double>
relaxation_update(std::vector<Eigen::Vector3d> &spins,
                  const std::vector<Eigen::Matrix3d> &conf_contribution,
                  DataPatternLexicographical &data_pattern) {
  std::vector<int> shift = {1, data_pattern.lat_dim[0],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2] *
                                data_pattern.lat_dim[3]};
  Eigen::Vector3d A;
  double diff_average = 0;
  double diff_maximal = 0;
  double diff_tmp;
  int position;
#pragma omp parallel for collapse(4) private(A, position)                      \
    firstprivate(data_pattern, shift) reduction(+ : diff_average)              \
    reduction(max : diff_maximal)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          position = data_pattern.get_index_site();
          A = contribution_site(spins, conf_contribution, data_pattern, x, y, z,
                                t, position, shift);
          diff_tmp = parallel(spins[position], A);
          if (diff_tmp > diff_maximal)
            diff_maximal = diff_tmp;
          diff_average += diff_tmp;
        }
      }
    }
  }
  diff_average = diff_average / data_pattern.get_lattice_size();
  return std::tuple<double, double>(diff_maximal, diff_average);
}

void make_simulated_annealing(std::vector<su2> &conf_su2,
                              std::vector<spin> &spins, double T_init,
                              double T_final, double T_step, int OR_steps,
                              int thermalization_steps) {
  for (int i = 0; i < thermalization_steps; i++) {
    heat_bath_update(spins, conf_su2, T_init);
  }
  double T = T_init;
  while (T > T_final) {
    heat_bath_update(spins, conf_su2, T);
    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(spins, conf_su2);
    }
    if (T <= 1.6 && T >= 1.2)
      T -= T_step / 4;
    else
      T -= T_step;
  }
}

void make_simulated_annealing(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf_su2,
    std::vector<spin> &spins, double T_init, double T_final, double T_step,
    int OR_steps, int thermalization_steps) {
  for (int i = 0; i < thermalization_steps; i++) {
    heat_bath_update_fast(spins, conf_su2, T_init);
  }
  double T = T_init;
  while (T > T_final) {
    heat_bath_update_fast(spins, conf_su2, T);
    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(spins, conf_su2);
    }
    if (T <= 1.6 && T >= 1.2)
      T -= T_step / 4;
    else
      T -= T_step;
  }
}

void make_simulated_annealing(
    const std::vector<Eigen::Matrix3d> &conf_contribution,
    std::vector<Eigen::Vector3d> &spins,
    DataPatternLexicographical &data_pattern, double T_init, double T_final,
    double T_step, int OR_steps, int thermalization_steps) {
  for (int i = 0; i < thermalization_steps; i++) {
    heat_bath_update(spins, conf_contribution, data_pattern, T_init);
  }
  double T = T_init;
  while (T >= 1.6) {
    heat_bath_update(spins, conf_contribution, data_pattern, T);
    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(spins, conf_contribution, data_pattern);
    }
    T -= T_step;
  }
  while (T >= 1.2) {
    heat_bath_update_fast(spins, conf_contribution, data_pattern, T);
    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(spins, conf_contribution, data_pattern);
    }
    T -= T_step / 4;
  }
  while (T > T_final) {
    heat_bath_update_fast(spins, conf_contribution, data_pattern, T);
    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(spins, conf_contribution, data_pattern);
    }
    T -= T_step;
  }
}

std::map<double, double> simulated_annealing_thermalization_test(
    std::vector<su2> &conf_su2, std::vector<spin> &spins, double T_init,
    double T_final, double T_step, int OR_steps, int thermalization_steps,
    int local_thermalization_steps) {
  std::map<double, double> result;
  for (int i = 0; i < thermalization_steps; i++) {
    heat_bath_update(spins, conf_su2, T_init);
  }
  double T = T_init;
  while (T > T_final) {
    for (int i = 0; i < local_thermalization_steps; i++) {
      heat_bath_update(spins, conf_su2, T);
      for (int i = 0; i < OR_steps; i++) {
        overrelaxation_update(spins, conf_su2);
      }
    }
    result[T] = MAG_functional_su2_spin(conf_su2, spins);
    T -= T_step;
  }
  return result;
}

std::map<double, double> simulated_annealing_thermalization_test(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf_su2,
    std::vector<spin> &spins, double T_init, double T_final, double T_step,
    int OR_steps, int thermalization_steps, int local_thermalization_steps) {
  std::map<double, double> result;
  for (int i = 0; i < thermalization_steps; i++) {
    heat_bath_update_fast(spins, conf_su2, T_init);
  }
  double T = T_init;
  while (T > T_final) {
    for (int i = 0; i < local_thermalization_steps; i++) {
      heat_bath_update_fast(spins, conf_su2, T);
      for (int i = 0; i < OR_steps; i++) {
        overrelaxation_update(spins, conf_su2);
      }
    }
    result[T] = MAG_functional_su2_spin(conf_su2, spins);
    T -= T_step;
  }
  return result;
}

std::map<double, double> simulated_annealing_thermalization_test(
    const std::vector<Eigen::Matrix3d> &conf_contribution,
    std::vector<Eigen::Vector3d> &spins,
    DataPatternLexicographical &data_pattern, double T_init, double T_final,
    double T_step, int OR_steps, int thermalization_steps,
    int local_thermalization_steps) {
  std::map<double, double> result;
  for (int i = 0; i < thermalization_steps; i++) {
    heat_bath_update_fast(spins, conf_contribution, data_pattern, T_init);
  }
  double T = T_init;
  while (T >= T_final) {
    for (int i = 0; i < local_thermalization_steps; i++) {
      heat_bath_update_fast(spins, conf_contribution, data_pattern, T);
      for (int j = 0; j < OR_steps; j++) {
        overrelaxation_update(spins, conf_contribution, data_pattern);
      }
    }
    result[T] = MAG_functional_su2_spin(conf_contribution, spins, data_pattern);
    T -= T_step;
  }
  return result;
}

void make_maximization_approximate(std::vector<su2> &conf_su2,
                                   std::vector<spin> &spins, int OR_steps,
                                   int tolerance_digits) {

  double functional = MAG_functional_su2_spin(conf_su2, spins);
  double functional_old = functional;

  bool is_equal = false;

  do {

    relaxation_update(spins, conf_su2);

    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(spins, conf_su2);
    }

    functional = MAG_functional_su2_spin(conf_su2, spins);

    is_equal = (trunc(powf(10., tolerance_digits) * functional) ==
                trunc(powf(10., tolerance_digits) * functional_old)) &&
               (signbit(functional) == signbit(functional_old));

    functional_old = functional;

  } while (!is_equal);
}

void make_maximization_approximate(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf_su2,
    std::vector<spin> &spins, int OR_steps, int tolerance_digits) {
  double functional = MAG_functional_su2_spin(conf_su2, spins);
  double functional_old = functional;
  bool is_equal = false;
  do {
    relaxation_update(spins, conf_su2);
    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(spins, conf_su2);
    }
    functional = MAG_functional_su2_spin(conf_su2, spins);
    is_equal = (trunc(powf(10., tolerance_digits) * functional) ==
                trunc(powf(10., tolerance_digits) * functional_old)) &&
               (signbit(functional) == signbit(functional_old));
    functional_old = functional;
  } while (!is_equal);
}

void make_maximization_approximate(
    const std::vector<Eigen::Matrix3d> &conf_contribution,
    std::vector<Eigen::Vector3d> &spins,
    DataPatternLexicographical &data_pattern, int OR_steps,
    int tolerance_digits) {
  double functional =
      MAG_functional_su2_spin(conf_contribution, spins, data_pattern);
  double functional_old = functional;
  bool is_equal = false;
  do {
    relaxation_update(spins, conf_contribution, data_pattern);
    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(spins, conf_contribution, data_pattern);
    }
    functional =
        MAG_functional_su2_spin(conf_contribution, spins, data_pattern);
    is_equal = (trunc(powf(10., tolerance_digits) * functional) ==
                trunc(powf(10., tolerance_digits) * functional_old)) &&
               (signbit(functional) == signbit(functional_old));
    functional_old = functional;
  } while (!is_equal);
}

void make_maximization_final(std::vector<su2> &conf_su2,
                             std::vector<spin> &spins, int OR_steps,
                             double tolerance_maximal,
                             double tolerance_average) {

  double functional = MAG_functional_su2_spin(conf_su2, spins);
  double functional_old = functional;

  bool is_converged = false;
  std::tuple<double, double> difference;

  do {

    difference = relaxation_update(spins, conf_su2);

    is_converged = (std::get<0>(difference) < tolerance_maximal) &&
                   (std::get<1>(difference) < tolerance_average);

    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(spins, conf_su2);
    }

  } while (!is_converged);
}

void make_maximization_final(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf_su2,
    std::vector<spin> &spins, int OR_steps, double tolerance_maximal,
    double tolerance_average) {
  bool is_converged = false;
  std::tuple<double, double> difference;
  do {
    difference = relaxation_update(spins, conf_su2);
    // std::cout << "difference: " << std::get<0>(difference) << " "
    //           << std::get<1>(difference) << std::endl;
    is_converged = (std::get<0>(difference) < tolerance_maximal) &&
                   (std::get<1>(difference) < tolerance_average);
    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(spins, conf_su2);
    }
  } while (!is_converged);
}

void make_maximization_final(
    const std::vector<Eigen::Matrix3d> &conf_contribution,
    std::vector<Eigen::Vector3d> &spins,
    DataPatternLexicographical &data_pattern, int OR_steps,
    double tolerance_maximal, double tolerance_average) {
  bool is_converged = false;
  std::tuple<double, double> difference;
  do {
    difference = relaxation_update(spins, conf_contribution, data_pattern);
    std::cout << "difference: " << std::get<0>(difference) << " "
              << std::get<1>(difference) << std::endl;
    is_converged = (std::get<0>(difference) < tolerance_maximal) &&
                   (std::get<1>(difference) < tolerance_average);
    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(spins, conf_contribution, data_pattern);
    }
  } while (!is_converged);
}

double mag_functional_su3(std::vector<su3> &conf_su3) {
  double functional = 0;
  std::vector<su3> generators_su3 = get_generators_su3();
  for (int i = 0; i < conf_su3.size(); i++) {
    for (int j = 0; j < 3; j++) {
      functional +=
          conf_su3[i].matrix(j, j).real() * conf_su3[i].matrix(j, j).real() +
          conf_su3[i].matrix(j, j).imag() * conf_su3[i].matrix(j, j).imag();
    }
  }
  return functional / (x_size * y_size * z_size * t_size * 3 * 4);
}

double mag_functional_su3(
    const Data::LatticeData<DataPatternLexicographical, su3> &conf_su3) {
  double functional = 0;
  std::vector<su3> generators_su3 = get_generators_su3();
  for (int i = 0; i < conf_su3.array.size(); i++) {
    for (int j = 0; j < 3; j++) {
      functional +=
          conf_su3[i].matrix(j, j).real() * conf_su3[i].matrix(j, j).real() +
          conf_su3[i].matrix(j, j).imag() * conf_su3[i].matrix(j, j).imag();
    }
  }
  return functional / (x_size * y_size * z_size * t_size * 3 * 4);
}