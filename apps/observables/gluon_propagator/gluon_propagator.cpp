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

  Data::data<MATRIX> conf;

  get_data(conf, conf_path, conf_format, bytes_skip, convert);

  std::cout << "plaket: " << plaket(conf.array) << std::endl;

  ofstream stream;
  stream.precision(17);
  // open file
  //   stream.open(path);
  start_time = omp_get_wtime();

  std::vector<std::array<double, 12>> vector_potential =
      get_vector_potential(conf.array);
  end_time = omp_get_wtime();
  observables_time = end_time - start_time;
  cout << "get_vector_potential time: " << observables_time << endl;

  for (int i = 0; i < 5; i++) {
    std::cout << vector_potential[0][i] << std::endl;
  }

  double multiplyer_t = 2 * M_PI / t_size;
  double multiplyer_s = 2 * M_PI / x_size;
  std::array<double, 4> momenta = {multiplyer_t, 2 * multiplyer_s,
                                   3 * multiplyer_s, 4 * multiplyer_s};

  start_time = omp_get_wtime();
  std::vector<std::complex<double>> furier_coefficients =
      get_furier_coefficients(momenta);
  //   for (int i = 0; i < 5; i++) {
  //     std::cout << furier_coefficients[i] << std::endl;
  //   }
  end_time = omp_get_wtime();
  observables_time = end_time - start_time;
  cout << "get_furier_coefficients time: " << observables_time << endl;

  start_time = omp_get_wtime();

  std::array<std::complex<double>, 144> gluon_propagator =
      calculate_gluon_propagator(vector_potential, furier_coefficients);

  end_time = omp_get_wtime();
  observables_time = end_time - start_time;
  cout << "gluon propagator time: " << observables_time << endl;
  for (int mu = 0; mu < 4; mu++) {
    for (int a = 0; a < 3; a++) {
      for (int nu = 0; nu < 4; nu++) {
        for (int b = 0; b < 3; b++) {
          std::cout << "mu: " << mu << ", nu: " << nu << ", a: " << a
                    << ", b: " << b << ", D:"
                    << gluon_propagator[(mu * 4 + a) * 12 + nu * 4 + b]
                    << std::endl;
        }
      }
    }
  }

  stream.close();
}