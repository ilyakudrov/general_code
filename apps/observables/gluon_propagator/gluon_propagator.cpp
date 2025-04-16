#include "../../../lib/cpu/include/gluon_propagator.h"
#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/matrix.h"

#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <sstream>

#ifndef MATRIX
#define MATRIX su2
#endif

using namespace std;

int x_size;
int y_size;
int z_size;
int t_size;
int size1;
int size2;

std::vector<su2> read_conf_test(std::string &file_name, int bytes_skip) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  std::vector<su2> array;
  array.reserve(4 * x_size * y_size * z_size * t_size);
  std::ifstream stream(file_name);
  std::vector<float> v(data_size * 4);
  stream.ignore(bytes_skip);
  if (!stream.read((char *)&v[0], (data_size * 4) * sizeof(float)))
    std::cout << "read_float<su2> error: " << file_name << std::endl;
  su2 A;
  std::array<int, 4> lat_dim = {x_size, y_size, z_size, t_size};
  std::array<int, 3> multipliers = {
      lat_dim[0] * lat_dim[1], lat_dim[0] * lat_dim[1] * lat_dim[2],
      lat_dim[0] * lat_dim[1] * lat_dim[2] * lat_dim[3]};
  int index;
  for (int t = 0; t < lat_dim[3]; t++) {
    for (int z = 0; z < lat_dim[2]; z++) {
      for (int y = 0; y < lat_dim[1]; y++) {
        for (int x = 0; x < lat_dim[0]; x++) {
          for (int mu = 0; mu < 4; mu++) {
            index = multipliers[1] * t + multipliers[0] * z + lat_dim[0] * y +
                    x + mu * multipliers[2];
            A.a0 = (double)v[index * 4];
            A.a1 = (double)v[index * 4 + 1];
            A.a2 = (double)v[index * 4 + 2];
            A.a3 = (double)v[index * 4 + 3];
            array.push_back(A);
          }
        }
      }
    }
  }
  stream.close();
  return array;
}

void write_float_test(std::string &file_name, std::vector<su2> &data) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  std::ofstream stream(file_name);
  std::vector<float> v(data_size1 * 4);
  std::array<int, 4> lat_dim = {x_size, y_size, z_size, t_size};
  std::array<int, 3> multipliers = {
      lat_dim[0] * lat_dim[1], lat_dim[0] * lat_dim[1] * lat_dim[2],
      lat_dim[0] * lat_dim[1] * lat_dim[2] * lat_dim[3]};
  std::array<int, 4> lat_coord;
  int index_write;
  int index_data;
  for (int t = 0; t < lat_dim[3]; t++) {
    for (int z = 0; z < lat_dim[2]; z++) {
      for (int y = 0; y < lat_dim[1]; y++) {
        for (int x = 0; x < lat_dim[0]; x++) {
          lat_coord = {x, y, z, t};
          for (int mu = 0; mu < 4; mu++) {
            index_write = multipliers[1] * t + multipliers[0] * z +
                          lat_dim[0] * y + x + mu * multipliers[2];
            index_data = get_index_matrix(lat_coord, mu);
            v[index_write * 4] = (float)data[index_data].a0;
            v[index_write * 4 + 1] = (float)data[index_data].a1;
            v[index_write * 4 + 2] = (float)data[index_data].a2;
            v[index_write * 4 + 3] = (float)data[index_data].a3;
          }
        }
      }
    }
  }
  if (!stream.write((char *)&v[0], data_size1 * 4 * sizeof(float)))
    std::cout << "write_float<su2> error: " << file_name << std::endl;
  stream.close();
}

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double observables_time;

  string conf_format;
  string conf_path;
  string path;
  int L_spat, L_time;
  int bytes_skip = 0;
  bool convert = 0;
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "-bytes_skip") {
      bytes_skip = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-conf_path") {
      conf_path = argv[++i];
    } else if (string(argv[i]) == "-convert") {
      istringstream(string(argv[++i])) >> convert;
    } else if (string(argv[i]) == "-L_spat") {
      L_spat = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-L_time") {
      L_time = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-path") {
      path = argv[++i];
    }
  }

  x_size = L_spat;
  y_size = L_spat;
  z_size = L_spat;
  t_size = L_time;
  size1 = x_size * y_size;
  size2 = x_size * y_size * z_size;

  cout << "conf_format " << conf_format << endl;
  cout << "conf_path " << conf_path << endl;
  cout << "bytes_skip " << bytes_skip << endl;
  cout << "convert " << convert << endl;
  cout << "L_spat " << L_spat << endl;
  cout << "L_time " << L_time << endl;
  cout << "path " << path << endl;
  cout << endl;

  cout.precision(17);

  // Data::data<MATRIX> conf;

  // get_data(conf, conf_path, conf_format, bytes_skip, convert);
  std::vector<su2> array = read_conf_test(conf_path, 0);
  std::cout << "plaket " << plaket(array) << std::endl;
  std::cout << "plaket " << plaket_time(array) << std::endl;
  std::cout << "plaket " << plaket_space(array) << std::endl;

  // std::string write_path =
  //     "/home/ilya/soft/lattice/general_code/tests/confs/su2/gluodynamics/"
  //     "32^3x8/beta2.779/CONF0001_test";
  // write_float_test(write_path, conf.array);

  ofstream stream;
  stream.precision(17);
  // open file
  //   stream.open(path);
  start_time = omp_get_wtime();

  std::vector<std::array<double, 12>> vector_potential =
      get_vector_potential(array);
  end_time = omp_get_wtime();
  observables_time = end_time - start_time;
  cout << "get_vector_potential time: " << observables_time << endl;

  double multiplyer_t = 2 * M_PI / t_size;
  double multiplyer_s = 2 * M_PI / x_size;
  std::array<double, 4> momenta = {multiplyer_t, 2 * multiplyer_s,
                                   3 * multiplyer_s, 4 * multiplyer_s};

  start_time = omp_get_wtime();
  std::vector<std::complex<double>> furier_coefficients =
      get_furier_coefficients(momenta);
  end_time = omp_get_wtime();
  observables_time = end_time - start_time;
  cout << "get_furier_coefficients time: " << observables_time << endl;

  start_time = omp_get_wtime();

  // std::array<std::complex<double>, 144> gluon_propagator =
  //     calculate_gluon_propagator(vector_potential, furier_coefficients);

  // end_time = omp_get_wtime();
  // observables_time = end_time - start_time;
  // cout << "gluon propagator time: " << observables_time << endl;
  // for (int mu = 0; mu < 4; mu++) {
  //   for (int a = 0; a < 3; a++) {
  //     for (int nu = 0; nu < 4; nu++) {
  //       for (int b = 0; b < 3; b++) {
  //         std::cout << "mu: " << mu << ", nu: " << nu << ", a: " << a
  //                   << ", b: " << b << ", D:"
  //                   << gluon_propagator[(mu * 4 + a) * 12 + nu * 4 + b]
  //                   << std::endl;
  //       }
  //     }
  //   }
  // }
  std::array<std::complex<double>, 12> gluon_propagator_diagonal =
      calculate_gluon_propagator_diagonal(vector_potential,
                                          furier_coefficients);

  end_time = omp_get_wtime();
  observables_time = end_time - start_time;
  cout << "gluon propagator time: " << observables_time << endl;
  for (int mu = 0; mu < 4; mu++) {
    std::cout << "mu = " << mu << ": " << gluon_propagator_diagonal[mu * 3]
              << " " << gluon_propagator_diagonal[mu * 3 + 1] << " "
              << gluon_propagator_diagonal[mu * 3 + 2] << std::endl;
  }

  // stream.close();
}