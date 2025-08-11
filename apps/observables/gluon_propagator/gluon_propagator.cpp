#include "../../../lib/cpu/include/gluon_propagator.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/plaket.h"

#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
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
  string file_precision;
  string conf_path;
  string output_path;
  int L_spat, L_time;
  int bytes_skip = 0;
  bool convert = 0;
  double beta = 0;
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "--conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "--file_precision") {
      file_precision = argv[++i];
    } else if (string(argv[i]) == "--bytes_skip") {
      bytes_skip = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--conf_path") {
      conf_path = argv[++i];
    } else if (string(argv[i]) == "--convert") {
      istringstream(string(argv[++i])) >> convert;
    } else if (string(argv[i]) == "--L_spat") {
      L_spat = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--L_time") {
      L_time = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--output_path") {
      output_path = argv[++i];
    } else if (string(argv[i]) == "--beta") {
      beta = stod(string(argv[++i]));
    }
  }

  int x_size1 = L_spat;
  int y_size1 = L_spat;
  int z_size1 = L_spat;
  int t_size1 = L_time;

  cout << "conf_format " << conf_format << endl;
  cout << "file_precision " << file_precision << endl;
  cout << "conf_path " << conf_path << endl;
  cout << "bytes_skip " << bytes_skip << endl;
  cout << "convert " << convert << endl;
  cout << "L_spat " << L_spat << endl;
  cout << "L_time " << L_time << endl;
  cout << "output_path " << output_path << endl;
  cout << "beta " << beta << endl;
  cout << endl;

  cout.precision(17);

  Data::LatticeData<DataPatternLexicographical, MATRIX> conf(
      {x_size1, y_size1, z_size1, t_size1});
  Data::read_data_convert(conf, conf_path, conf_format, bytes_skip,
                          file_precision, convert);
  DataPatternLexicographical data_pattern(conf.lat_dim);

  double a_inv = (542.6 * 8 / 1000);
  double g = 2 / sqrt(beta);
  double multiplier = 2 / g * a_inv;

  std::cout << "plaket " << plaket(conf) << std::endl;
  std::cout << "plaket " << plaket_time(conf) << std::endl;
  std::cout << "plaket " << plaket_space(conf) << std::endl;

  std::vector<std::vector<std::array<double, 4>>> momenta =
      generate_momenta(x_size1, t_size1);

  std::vector<std::array<double, 12>> vector_potential =
      get_vector_potential(conf);

  // std::vector<std::complex<double>> furier_coefficients;

  std::map<std::tuple<double, double, double, double, int, int, int, int>,
           std::complex<double>>
      gluon_propagator_map;

  for (int i = 0; i < 100; i++) {
    start_time = omp_get_wtime();
    std::array<std::complex<double>, 144> gluon_propagator =
        calculate_gluon_propagator_group(vector_potential, momenta[i],
                                         beta / a_inv / a_inv, data_pattern);
    for (int j = 0; j < momenta[i].size(); j++) {
      for (int mu = 0; mu < 4; mu++) {
        for (int a = 0; a < 3; a++) {
          for (int nu = 0; nu < 4; nu++) {
            for (int b = 0; b < 3; b++) {
              gluon_propagator_map[std::make_tuple(
                  momenta[i][0][0], momenta[i][0][1], momenta[i][0][2],
                  momenta[i][0][3], mu, nu, a, b)] =
                  gluon_propagator[(mu * 3 + a) * 12 + nu * 3 + b];
            }
          }
        }
      }
    }
    end_time = omp_get_wtime();
    std::cout << "time: " << end_time - start_time << std::endl;
  }

  ofstream stream;
  stream.precision(17);
  stream.open(output_path);
  stream << "p1,p2,p3,p4,mu,nu,a,b,Dr,Di" << std::endl;
  for (const auto &pair : gluon_propagator_map) {
    stream << get<0>(pair.first) << "," << get<1>(pair.first) << ","
           << get<2>(pair.first) << "," << get<3>(pair.first) << ","
           << get<4>(pair.first) << "," << get<5>(pair.first) << ","
           << get<6>(pair.first) << "," << get<7>(pair.first) << ","
           << pair.second.real() << "," << pair.second.imag() << endl;
  }
  stream.close();
}