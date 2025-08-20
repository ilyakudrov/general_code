#include "../../lib/cpu/include/Landau_su2.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/matrix.h"
#include "../../lib/cpu/include/plaket.h"

#include <iostream>
#include <omp.h>

using namespace std;

// global variables for lattice size
int x_size;
int y_size;
int z_size;
int t_size;
int size1;
int size2;

int main(int argc, char **argv) {
  double omp_time;
  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;
  int x_size1;
  int y_size1;
  int z_size1;
  int t_size1;
  string path_conf;
  string conf_format;
  string file_precision;
  int bytes_skip = 0;
  double T_step;
  double T_init;
  double T_final;
  int OR_steps;
  int thermalization_steps;
  double tolerance_maximal;
  double tolerance_average;
  // read parameters
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "--conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "--file_precision") {
      file_precision = argv[++i];
    } else if (string(argv[i]) == "--path_conf") {
      path_conf = argv[++i];
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
    } else if (string(argv[i]) == "--tolerance_maximal") {
      tolerance_maximal = stod(string(argv[++i]));
    } else if (string(argv[i]) == "--tolerance_average") {
      tolerance_average = stod(string(argv[++i]));
    } else if (string(argv[i]) == "--thermalization_steps") {
      thermalization_steps = stoi(string(argv[++i]));
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
  cout << "bytes_skip " << bytes_skip << endl;
  cout << "OR_steps " << OR_steps << endl;
  cout << "tolerance_maximal " << tolerance_maximal << endl;
  cout << "tolerance_average " << tolerance_average << endl;
  cout << "T_step " << T_step << endl;
  cout << "T_init " << T_init << endl;
  cout << "T_final " << T_final << endl;
  cout << "thermalization_steps " << thermalization_steps << endl;
  cout << "x_size " << x_size1 << endl;
  cout << "y_size " << y_size1 << endl;
  cout << "z_size " << z_size1 << endl;
  cout << "t_size " << t_size1 << endl;

  Data::LatticeData<DataPatternLexicographical, su2> conf_su2(
      {x_size1, y_size1, z_size1, t_size1});
  Data::read_data_convert(conf_su2, path_conf, conf_format, bytes_skip,
                          file_precision, 0);
  DataPatternLexicographical data_pattern(conf_su2.lat_dim);
  cout.precision(17);
  std::vector<su2> gauge = generate_gauge_su2_uniform(data_pattern);
  std::cout << "plaket: " << plaket(conf_su2) << std::endl;
  std::cout << "Landau su2 functional: "
            << Landau_functional_conf(conf_su2, gauge) << std::endl;

  omp_time = omp_get_wtime();
  make_simulated_annealing(conf_su2, gauge, T_init, T_final, T_step, OR_steps,
                           thermalization_steps);
  make_maximization_final(conf_su2, gauge, OR_steps, tolerance_maximal,
                          tolerance_average);
  std::cout << "gauge fixing time: " << omp_get_wtime() - omp_time << std::endl;

  apply_gauge_Landau(gauge, conf_su2);
  std::cout << "Landau su2 functional: " << Landau_su2_functional(conf_su2)
            << std::endl;
  std::cout << "plaket: " << plaket(conf_su2) << std::endl;
}
