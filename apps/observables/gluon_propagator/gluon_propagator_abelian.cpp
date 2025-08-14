#include "../../../lib/cpu/include/Landau_U1.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/decomposition.h"
#include "../../../lib/cpu/include/gluon_propagator.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/monopoles.h"
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
  double omp_time;

  string conf_format;
  string file_precision;
  string conf_path;
  string path_inverse_laplacian;
  string output_path_abelian;
  string output_path_monopole;
  string output_path_photon;
  int L_spat, L_time;
  int bytes_skip = 0;
  bool convert = 0;
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
    } else if (string(argv[i]) == "--path_inverse_laplacian") {
      path_inverse_laplacian = argv[++i];
    } else if (string(argv[i]) == "--output_path_abelian") {
      output_path_abelian = argv[++i];
    } else if (string(argv[i]) == "--output_path_monopole") {
      output_path_monopole = argv[++i];
    } else if (string(argv[i]) == "--output_path_photon") {
      output_path_photon = argv[++i];
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
  cout << "path_inverse_laplacian " << path_inverse_laplacian << endl;
  cout << "output_path_abelian " << output_path_abelian << endl;
  cout << "output_path_monopole " << output_path_monopole << endl;
  cout << "output_path_photon " << output_path_photon << endl;
  cout << endl;

  cout.precision(17);

  Data::LatticeData<DataPatternLexicographical, abelian> conf(
      {x_size1, y_size1, z_size1, t_size1});
  Data::read_data_convert(conf, conf_path, conf_format, bytes_skip,
                          file_precision, convert);
  DataPatternLexicographical data_pattern(conf.lat_dim);
  DataPatternLexicographical data_pattern_laplacian(
      {conf.lat_dim[0] / 2 + 1, conf.lat_dim[1] / 2 + 1,
       conf.lat_dim[2] / 2 + 1, conf.lat_dim[3] / 2 + 1});

  std::cout << "plaket " << plaket(conf) << std::endl;
  std::cout << "plaket " << plaket_time(conf) << std::endl;
  std::cout << "plaket " << plaket_space(conf) << std::endl;

  std::vector<std::complex<double>> conf_complex = convert_to_complex(conf);
  std::vector<std::complex<double>> gauge_complex =
      generate_gauge_complex_uniform(data_pattern);
  std::cout << "Landau functional before Landau gauge: "
            << Landau_functional_conf_complex(conf_complex, gauge_complex,
                                              data_pattern)
            << std::endl;

  omp_time = omp_get_wtime();
  make_simulated_annealing(conf_complex, gauge_complex, data_pattern, 10, 0.1,
                           0.1, 4, 20);
  make_maximization_final(conf_complex, gauge_complex, data_pattern, 4, 1e-14,
                          1e-15);
  std::cout << "Landau gauge time: " << omp_get_wtime() - omp_time << std::endl;
  std::cout << "Landau functional after Landau gauge: "
            << Landau_functional_conf_complex(conf_complex, gauge_complex,
                                              data_pattern)
            << std::endl;

  conf.array = convert_complex_to_abelian(conf_complex, data_pattern);
  std::vector<std::vector<std::array<double, 4>>> momenta =
      generate_momenta(x_size1, t_size1);
  std::vector<std::array<double, 4>> vector_potential =
      get_vector_potential(conf);
  std::map<std::tuple<double, double, double, double, int>,
           std::complex<double>>
      gluon_propagator_map_abelian;
  std::array<std::complex<double>, 4> gluon_propagator;

  for (int i = 0; i < momenta.size(); i++) {
    gluon_propagator = calculate_gluon_propagator_diagonal_group(
        vector_potential, momenta[i], 1, data_pattern);
    for (int j = 0; j < momenta[i].size(); j++) {
      for (int mu = 0; mu < 4; mu++) {
        gluon_propagator_map_abelian[std::make_tuple(
            momenta[i][0][0], momenta[i][0][1], momenta[i][0][2],
            momenta[i][0][3], mu)] = gluon_propagator[mu];
      }
    }
  }

  omp_time = omp_get_wtime();
  vector<double> inverse_laplacian =
      read_inverse_laplacian(path_inverse_laplacian, data_pattern_laplacian);
  vector<vector<int>> dirac_plakets = calculate_monopole_plaket_singular(conf);
  vector<double> monopole_angles = make_monopole_angles(
      dirac_plakets, inverse_laplacian, data_pattern, data_pattern_laplacian);
  Data::LatticeData<DataPatternLexicographical, abelian> conf_monopole(
      {x_size1, y_size1, z_size1, t_size1});
  Data::LatticeData<DataPatternLexicographical, abelian> conf_photon(
      {x_size1, y_size1, z_size1, t_size1});
  for (int i = 0; i < monopole_angles.size(); i++) {
    conf_monopole.array[i] = abelian(1, monopole_angles[i]);
    conf_photon.array[i] = abelian(1, conf[i].phi + monopole_angles[i]);
  }
  monopole_angles.clear();
  monopole_angles.shrink_to_fit();
  dirac_plakets.clear();
  dirac_plakets.shrink_to_fit();
  std::cout << "decomposition time: " << omp_get_wtime() - omp_time
            << std::endl;

  vector_potential = get_vector_potential(conf_monopole);
  std::map<std::tuple<double, double, double, double, int>,
           std::complex<double>>
      gluon_propagator_map_monopole;
  for (int i = 0; i < momenta.size(); i++) {
    gluon_propagator = calculate_gluon_propagator_diagonal_group(
        vector_potential, momenta[i], 1, data_pattern);
    for (int j = 0; j < momenta[i].size(); j++) {
      for (int mu = 0; mu < 4; mu++) {
        gluon_propagator_map_monopole[std::make_tuple(
            momenta[i][0][0], momenta[i][0][1], momenta[i][0][2],
            momenta[i][0][3], mu)] = gluon_propagator[mu];
      }
    }
  }

  vector_potential = get_vector_potential(conf_photon);
  std::map<std::tuple<double, double, double, double, int>,
           std::complex<double>>
      gluon_propagator_map_photon;
  for (int i = 0; i < momenta.size(); i++) {
    gluon_propagator = calculate_gluon_propagator_diagonal_group(
        vector_potential, momenta[i], 1, data_pattern);
    for (int j = 0; j < momenta[i].size(); j++) {
      for (int mu = 0; mu < 4; mu++) {
        gluon_propagator_map_photon[std::make_tuple(
            momenta[i][0][0], momenta[i][0][1], momenta[i][0][2],
            momenta[i][0][3], mu)] = gluon_propagator[mu];
      }
    }
  }

  ofstream stream_abelian;
  stream_abelian.precision(17);
  stream_abelian.open(output_path_abelian);
  stream_abelian << "p1,p2,p3,p4,mu,Dr,Di" << std::endl;
  for (const auto &pair : gluon_propagator_map_abelian) {
    stream_abelian << get<0>(pair.first) << "," << get<1>(pair.first) << ","
                   << get<2>(pair.first) << "," << get<3>(pair.first) << ","
                   << get<4>(pair.first) << "," << pair.second.real() << ","
                   << pair.second.imag() << endl;
  }
  stream_abelian.close();
  ofstream stream_monopole;
  stream_monopole.precision(17);
  stream_monopole.open(output_path_monopole);
  stream_monopole << "p1,p2,p3,p4,mu,Dr,Di" << std::endl;
  for (const auto &pair : gluon_propagator_map_monopole) {
    stream_monopole << get<0>(pair.first) << "," << get<1>(pair.first) << ","
                    << get<2>(pair.first) << "," << get<3>(pair.first) << ","
                    << get<4>(pair.first) << "," << pair.second.real() << ","
                    << pair.second.imag() << endl;
  }
  stream_monopole.close();
  ofstream stream_photon;
  stream_photon.precision(17);
  stream_photon.open(output_path_photon);
  stream_photon << "p1,p2,p3,p4,mu,Dr,Di" << std::endl;
  for (const auto &pair : gluon_propagator_map_photon) {
    stream_photon << get<0>(pair.first) << "," << get<1>(pair.first) << ","
                  << get<2>(pair.first) << "," << get<3>(pair.first) << ","
                  << get<4>(pair.first) << "," << pair.second.real() << ","
                  << pair.second.imag() << endl;
  }
  stream_photon.close();
}