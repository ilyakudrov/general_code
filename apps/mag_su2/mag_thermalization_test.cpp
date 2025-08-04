#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/mag.h"
#include "../../lib/cpu/include/matrix.h"
#include "../../lib/cpu/include/plaket.h"

#include <fstream>
#include <iostream>
#include <map>
#include <omp.h>

using namespace std;

vector<Eigen::Vector3d> spins_to_vectors(vector<spin> &spins) {
  vector<Eigen::Vector3d> vectors(spins.size());
  for (int i = 0; i < spins.size(); i++) {
    vectors[i] = Eigen::Vector3d(spins[i].a1, spins[i].a2, spins[i].a3);
  }
  return vectors;
}

// global variables for lattice size
int x_size;
int y_size;
int z_size;
int t_size;
int size1;
int size2;

int main(int argc, char **argv) {
  double omp_time;

  int x_size1;
  int y_size1;
  int z_size1;
  int t_size1;

  string path_conf;
  string conf_format;
  string file_precision;
  int bytes_skip = 0;
  string path_functional_output;

  double T_step;
  double T_init;
  double T_final;
  int OR_steps;
  int thermalization_steps;
  int local_thermalization_steps;

  // read parameters
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "--conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "--file_precision") {
      file_precision = argv[++i];
    } else if (string(argv[i]) == "--path_conf") {
      path_conf = argv[++i];
    } else if (string(argv[i]) == "--path_functional_output") {
      path_functional_output = argv[++i];
    } else if (string(argv[i]) == "--bytes_skip") {
      bytes_skip = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--T_step") {
      T_step = stod(string(argv[++i]));
    } else if (string(argv[i]) == "--T_init") {
      T_init = stod(string(argv[++i]));
    } else if (string(argv[i]) == "--T_final") {
      T_final = stod(string(argv[++i]));
    } else if (string(argv[i]) == "--OR_steps") {
      OR_steps = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--thermalization_steps") {
      thermalization_steps = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--local_thermalization_steps") {
      local_thermalization_steps = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--x_size") {
      x_size1 = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--y_size") {
      y_size1 = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--z_size") {
      z_size1 = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--t_size") {
      t_size1 = stoi(string(argv[++i]));
    } else
      cout << "unknown parameter " << argv[i] << endl;
  }

  cout << "path_conf " << path_conf << endl;
  cout << "conf_format " << conf_format << endl;
  cout << "file_precision " << file_precision << endl;
  cout << "path_functional_output " << path_functional_output << endl;
  cout << "bytes_skip " << bytes_skip << endl;
  cout << "OR_steps " << OR_steps << endl;
  cout << "T_step " << T_step << endl;
  cout << "T_init " << T_init << endl;
  cout << "T_final " << T_final << endl;
  cout << "thermalization_steps " << thermalization_steps << endl;
  cout << "local_thermalization_steps " << local_thermalization_steps << endl;
  cout << "x_size " << x_size1 << endl;
  cout << "y_size " << y_size1 << endl;
  cout << "z_size " << z_size1 << endl;
  cout << "t_size " << t_size1 << endl;

  // x_size = x_size1;
  // y_size = y_size1;
  // z_size = z_size1;
  // t_size = t_size1;

  Data::LatticeData<DataPatternLexicographical, su2> conf_su2(
      {x_size1, y_size1, z_size1, t_size1});
  Data::read_data_convert(conf_su2, path_conf, conf_format, bytes_skip,
                          file_precision, 0);
  DataPatternLexicographical data_pattern(conf_su2.lat_dim);
  cout.precision(17);
  vector<spin> spins = generate_spins_uniform(data_pattern);
  // vector<spin> spins = generate_spins_uniform();
  // vector<Eigen::Vector3d> vectors = spins_to_vectors(spins);
  // vector<Eigen::Vector3d> spins = generate_vectors_uniform(data_pattern);
  std::vector<Eigen::Matrix3d> conf_contribution =
      make_conf_contribution(conf_su2);
  std::cout << plaket(conf_su2) << std::endl;

  omp_time = omp_get_wtime();

  // std::map<double, double> functional =
  // simulated_annealing_thermalization_test(
  //     conf_contribution, spins, data_pattern, T_init, T_final, T_step,
  //     OR_steps, thermalization_steps, local_thermalization_steps);
  std::map<double, double> functional = simulated_annealing_thermalization_test(
      conf_su2, spins, T_init, T_final, T_step, OR_steps, thermalization_steps,
      local_thermalization_steps);

  std::cout << "thermalization test time: " << omp_get_wtime() - omp_time
            << std::endl;

  ofstream functional_stream;
  functional_stream.open(path_functional_output);
  functional_stream.precision(17);
  functional_stream << "temperature,functional" << endl;
  for (auto pair : functional) {
    functional_stream << pair.first << "," << pair.second << endl;
  }
  functional_stream.close();
}
