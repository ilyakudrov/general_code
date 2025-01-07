#include "../include/mag.h"
#include "../include/link.h"
#include "../include/matrix.h"

#include <fstream>
#include <iostream>
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

spin contribution_site(std::vector<spin> &spins, std::vector<su2> &conf_su2,
                       int x, int y, int z, int t, int position,
                       std::vector<int> &shift) {

  spin A(0, 0, 0);

  // mu = 0
  if (x < x_size - 1)
    A.contribution1(conf_su2[position * 4], spins[position + shift[0]]);
  else
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

void heat_bath_update(std::vector<spin> &spins, std::vector<su2> &conf_su2,
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

  spin A(0, 0, 0);
  int position = 0;
  int count = 0;

  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {

          A = contribution_site(spins, conf_su2, x, y, z, t, position, shift);

          heat_bath(spins[position], A, temperature,
                    &random_numbers[count * 3]);

          position++;
          count++;
        }
      }
    }
  }

  normalize_spin(spins);
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

double MAG_functional_su2(const std::vector<su2> &array) {
  double result = 0;
  for (int i = 0; i < x_size * y_size * z_size * t_size * 4; i++) {
    result += (array[i].sigma3_mult() * array[i].conj()).tr();
  }
  return result / (x_size * y_size * z_size * t_size * 4);
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

double mag_functional_su3(std::vector<su3> &conf_su3) {
  double functional = 0;
  std::vector<su3> generators_su3 = get_generators_su3();
  std::cout << generators_su3[2] << std::endl;
  std::cout << generators_su3[7] << std::endl;
  for (int i = 0; i < conf_su3.size(); i++) {
    // functional +=
    //     (((conf_su3[i] * generators_su3[7]) ^ conf_su3[i]) *
    //     generators_su3[7])
    //         .tr();
    // (((conf_su3[i] * generators_su3[2]) ^ conf_su3[i]) * generators_su3[2])
    //     .tr() +
    // (((conf_su3[i] * generators_su3[7]) ^ conf_su3[i]) * generators_su3[7])
    //     .tr();
    for (int j = 0; j < 3; j++) {
      functional +=
          conf_su3[i].matrix(j, j).real() * conf_su3[i].matrix(j, j).real() +
          conf_su3[i].matrix(j, j).imag() * conf_su3[i].matrix(j, j).imag();
    }
  }
  return functional / (x_size * y_size * z_size * t_size * 3 * 4);
}