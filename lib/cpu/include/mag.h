#pragma once

#include "../include/matrix.h"

#include <vector>

std::vector<spin> read_spins(std::string spins_path);

void write_spins(std::string output_path, std::vector<spin> spins);

std::vector<double> generate_random_numbers_sphere(int vector_size);

std::vector<spin> generate_spins_uniform();

std::vector<double> generate_random_numbers(int vector_size);

void heat_bath(spin &spins, spin &neighbour, double temperature,
               double *random_numbers);

spin contribution_site(std::vector<spin> &spins, std::vector<su2> &conf_su2,
                       int x, int y, int z, int t, int position,
                       std::vector<int> &shift);

void heat_bath_update(std::vector<spin> &spins, std::vector<su2> &conf_su2,
                      double temperature);

void normalize_spin(std::vector<spin> &spins);

double MAG_functional_su2_spin(std::vector<su2> &conf_su2,
                               std::vector<spin> &spins);

double MAG_functional_su2(const std::vector<su2> &array);

std::vector<su2> make_gauge(std::vector<spin> &spins);

std::vector<su2> gauge_tranformation(std::vector<su2> &conf_su2,
                                     std::vector<su2> &gauge);

void gauge_tranformation_spins(std::vector<su2> &conf_su2,
                               std::vector<spin> &spins);

std::vector<int> make_indices_qube(int qube_size);

spin contribution_site2(std::vector<spin> &spins, std::vector<su2> &conf_su2,
                        std::vector<int> &indices, int index);

void heat_bath_update_tets2(std::vector<spin> &spins,
                            std::vector<su2> &conf_su2,
                            std::vector<int> indices, double temperature);

void heat_bath_update_tets3(std::vector<spin> &spins,
                            std::vector<su2> &conf_su2, double temperature);

void make_indices_qube1(std::vector<int> &indices,
                        std::vector<char> &coordinates, int qube_size);

void heat_bath_update_tets4(std::vector<spin> &spins,
                            std::vector<su2> &conf_su2,
                            std::vector<int> &indices,
                            std::vector<char> &coordinates, double temperature);

void overrelaxation_update(std::vector<spin> &spins,
                           std::vector<su2> &conf_su2);

std::tuple<double, double> relaxation_update(std::vector<spin> &spins,
                                             std::vector<su2> &conf_su2);

void make_simulated_annealing(std::vector<su2> &conf_su2,
                              std::vector<spin> &spins, double T_init,
                              double T_final, double T_step, int OR_steps,
                              int thermalization_steps);

void make_maximization_approximate(std::vector<su2> &conf_su2,
                                   std::vector<spin> &spins, int OR_steps,
                                   int tolerance_digits);

void make_maximization_final(std::vector<su2> &conf_su2,
                             std::vector<spin> &spins, int OR_steps,
                             double tolerance_maximal,
                             double tolerance_average);

double mag_functional_su3(std::vector<su3> &conf_su3);

template <class DataPattern>
std::vector<spin> generate_spins_uniform(DataPattern &data_pattern) {
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

// spin contribution_site(std::vector<spin> &spins, std::vector<su2> &conf_su2,
//                        int x, int y, int z, int t, int position,
//                        std::vector<int> &shift) {

//   spin A(0, 0, 0);

//   // mu = 0
//   if (x < x_size - 1)
//     A.contribution1(conf_su2[position * 4], spins[position + shift[0]]);
//   else
//     A.contribution1(conf_su2[position * 4],
//                     spins[position + shift[0] - shift[1]]);

//   if (x > 0)
//     A.contribution1_conj(conf_su2[(position - shift[0]) * 4],
//                          spins[position - shift[0]]);
//   else
//     A.contribution1_conj(conf_su2[(position - shift[0] + shift[1]) * 4],
//                          spins[position - shift[0] + shift[1]]);

//   // mu = 1
//   if (y < y_size - 1)
//     A.contribution1(conf_su2[position * 4 + 1], spins[position + shift[1]]);
//   else
//     A.contribution1(conf_su2[position * 4 + 1],
//                     spins[position + shift[1] - shift[2]]);

//   if (y > 0)
//     A.contribution1_conj(conf_su2[(position - shift[1]) * 4 + 1],
//                          spins[position - shift[1]]);
//   else
//     A.contribution1_conj(conf_su2[(position - shift[1] + shift[2]) * 4 + 1],
//                          spins[position - shift[1] + shift[2]]);

//   // mu = 2
//   if (z < z_size - 1)
//     A.contribution1(conf_su2[position * 4 + 2], spins[position + shift[2]]);
//   else
//     A.contribution1(conf_su2[position * 4 + 2],
//                     spins[position + shift[2] - shift[3]]);

//   if (z > 0)
//     A.contribution1_conj(conf_su2[(position - shift[2]) * 4 + 2],
//                          spins[position - shift[2]]);
//   else
//     A.contribution1_conj(conf_su2[(position - shift[2] + shift[3]) * 4 + 2],
//                          spins[position - shift[2] + shift[3]]);

//   // mu = 3
//   if (t < t_size - 1)
//     A.contribution1(conf_su2[position * 4 + 3], spins[position + shift[3]]);
//   else
//     A.contribution1(conf_su2[position * 4 + 3],
//                     spins[position + shift[3] - shift[4]]);

//   if (t > 0)
//     A.contribution1_conj(conf_su2[(position - shift[3]) * 4 + 3],
//                          spins[position - shift[3]]);
//   else
//     A.contribution1_conj(conf_su2[(position - shift[3] + shift[4]) * 4 + 3],
//                          spins[position - shift[3] + shift[4]]);

//   return A;
// }

// template <class DataPattern>
// void heat_bath_update(std::vector<spin> &spins, std::vector<su2> &conf_su2,
// DataPattern &data_pattern,
//                       double temperature) {
//   std::vector<double> random_numbers;
//   random_numbers =
//       generate_random_numbers(data_pattern.get_lattice_size * 3);
//   spin A(0, 0, 0);
//   int position = 0;
//   int count = 0;
// #pragma omp parallel for collapse(4) firstprivate(data_pattern)
// private(place)
//   for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
//     for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
//       for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
//         for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
//           data_pattern.lat_coord = {x, y, z, t};

//           A = contribution_site(spins, conf_su2, x, y, z, t, position,
//           shift);

//           heat_bath(spins[position], A, temperature,
//                     &random_numbers[count * 3]);

//           position++;
//           count++;
//         }
//       }
//     }
//   }

//   normalize_spin(spins);
// }