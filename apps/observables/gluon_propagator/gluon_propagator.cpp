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

  double a_inv = (542.6 * 8 / 1000);
  double beta = 2.701;
  double g = 2 / sqrt(beta);
  double multiplier = 2 / g * a_inv;

  get_data(conf, conf_path, conf_format, bytes_skip, convert);
  // std::vector<su2> array = read_conf_test(conf_path, 0);
  std::cout << "plaket " << plaket(conf.array) << std::endl;
  std::cout << "plaket " << plaket_time(conf.array) << std::endl;
  std::cout << "plaket " << plaket_space(conf.array) << std::endl;

  ofstream stream;
  stream.precision(17);
  // open file
  //   stream.open(path);

  std::vector<std::array<double, 12>> vector_potential =
      get_vector_potential(conf.array);

  std::vector<std::complex<double>> furier_coefficients;

  double multiplyer_t = 2 * M_PI / t_size;
  double multiplyer_s = 2 * M_PI / x_size;
  // std::array<double, 4> momenta = {multiplyer_t, 2 * multiplyer_s,
  //                                  3 * multiplyer_s, 4 * multiplyer_s};
  std::array<double, 4> momenta = {multiplyer_s, 0, 0, 0};
  // std::array<double, 4> momenta = {multiplyer_s, multiplyer_s, multiplyer_s,
  // 0}; std::array<double, 4> momenta = {0, 0, 0, 0};

  double propagator_transversal_aver_momentum = 0;
  double propagator_longitudinal_aver_momentum = 0;
  for (int i = 0; i < 3; i++) {
    momenta = {0, 0, 0, 0};
    momenta[i] = multiplyer_s;
    furier_coefficients = get_furier_coefficients(momenta);

    start_time = omp_get_wtime();

    std::array<std::complex<double>, 144> gluon_propagator =
        calculate_gluon_propagator(vector_potential, furier_coefficients,
                                   beta / a_inv / a_inv);

    end_time = omp_get_wtime();
    observables_time = end_time - start_time;
    cout << "gluon propagator time: " << observables_time << endl;
    for (int mu = 0; mu < 4; mu++) {
      for (int a = 0; a < 3; a++) {
        for (int nu = 0; nu < 4; nu++) {
          for (int b = 0; b < 3; b++) {
            std::cout << "mu: " << mu << ", nu: " << nu << ", a: " << a
                      << ", b: " << b << ", D:"
                      << gluon_propagator[(mu * 3 + a) * 12 + nu * 3 + b] << " "
                      << (mu * 3 + a) * 12 + nu * 3 + b << std::endl;
          }
        }
      }
    }

    for (int mu = 0; mu < 3; mu++) {
      for (int a = 0; a < 3; a++) {
        propagator_transversal_aver_momentum +=
            gluon_propagator[(mu * 3 + a) * 12 + mu * 3 + a].real();
      }
    }
    for (int mu = 3; mu < 4; mu++) {
      for (int a = 0; a < 3; a++) {
        propagator_longitudinal_aver_momentum +=
            gluon_propagator[(mu * 3 + a) * 12 + mu * 3 + a].real();
      }
    }
  }

  std::cout << "propagator_transversal_aver_momentum "
            << propagator_transversal_aver_momentum / 6 / 3 << std::endl;
  std::cout << "propagator_longitudinal_aver_momentum "
            << propagator_longitudinal_aver_momentum / 3 / 3 << std::endl;

  int m = 12;
  int n = 12;

  std::cout << "propagator single: "
            << calculate_gluon_propagator_single(
                   vector_potential, furier_coefficients, multiplier, m, n)
            << std::endl;

  propagator_transversal_aver_momentum = 0;
  propagator_longitudinal_aver_momentum = 0;
  for (int i = 0; i < 3; i++) {

    momenta = {0, 0, 0, 0};
    momenta[i] = multiplyer_s;
    furier_coefficients = get_furier_coefficients(momenta);

    start_time = omp_get_wtime();

    std::array<std::complex<double>, 12> gluon_propagator_diagonal =
        calculate_gluon_propagator_diagonal(
            vector_potential, furier_coefficients, beta / a_inv / a_inv);

    end_time = omp_get_wtime();
    observables_time = end_time - start_time;
    cout << "gluon propagator diagonal time: " << observables_time << endl;
    for (int mu = 0; mu < 4; mu++) {
      std::cout << "mu = " << mu << ": " << gluon_propagator_diagonal[mu * 3]
                << " " << gluon_propagator_diagonal[mu * 3 + 1] << " "
                << gluon_propagator_diagonal[mu * 3 + 2] << std::endl;
    }

    double propagator_transversal_aver = 0;
    for (int i = 0; i < 9; i++) {
      propagator_transversal_aver += gluon_propagator_diagonal[i].real();
    }
    propagator_transversal_aver /= 6;
    propagator_transversal_aver_momentum += propagator_transversal_aver;
    std::cout << "propagator_transversal: " << propagator_transversal_aver
              << std::endl;

    double propagator_longitudinal_aver = 0;
    for (int i = 9; i < 12; i++) {
      propagator_longitudinal_aver += gluon_propagator_diagonal[i].real();
    }
    propagator_longitudinal_aver /= 3;
    propagator_longitudinal_aver_momentum += propagator_longitudinal_aver;
    std::cout << "propagator_longitudinal: " << propagator_longitudinal_aver
              << std::endl;
  }

  propagator_transversal_aver_momentum /= 3;
  propagator_longitudinal_aver_momentum /= 3;
  std::cout << "propagator_transversal_aver_momentum: "
            << propagator_transversal_aver_momentum << std::endl;
  std::cout << "propagator_longitudinal_aver_momentum: "
            << propagator_longitudinal_aver_momentum << std::endl;

  // stream.close();
}