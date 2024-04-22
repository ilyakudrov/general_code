#include "../../lib/cpu/include/basic_observables.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/matrix.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <omp.h>

#ifndef MATRIX
#define MATRIX su2
#endif

using namespace std;

// global variables for lattice size
int x_size;
int y_size;
int z_size;
int t_size;

int main(int argc, char **argv) {

  double start_time;
  double end_time;
  double search_time;

  string path_conf;
  string conf_format;
  int bytes_skip = 0;
  string path_output_correlator;
  string correlator_type;
  bool convert = 0;

  int D_max;

  // read parameters
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "-path_conf") {
      path_conf = argv[++i];
    } else if (string(argv[i]) == "-path_output_correlator") {
      path_output_correlator = argv[++i];
    } else if (string(argv[i]) == "-correlator_type") {
      correlator_type = argv[++i];
    } else if (string(argv[i]) == "-bytes_skip") {
      bytes_skip = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-convert") {
      istringstream(string(argv[++i])) >> convert;
    } else if (string(argv[i]) == "-D_max") {
      D_max = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-x_size") {
      x_size = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-y_size") {
      y_size = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-z_size") {
      z_size = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-t_size") {
      t_size = stoi(string(argv[++i]));
    } else
      cout << "unknown parameter " << argv[i] << endl;
  }

  cout << "path_conf " << path_conf << endl;
  cout << "conf_format " << conf_format << endl;
  cout << "bytes_skip " << bytes_skip << endl;
  cout << "D_max " << D_max << endl;
  cout << "convert " << convert << endl;

  cout << "path_output_correlator " << path_output_correlator << endl;

  cout << "x_size " << x_size << endl;
  cout << "y_size " << y_size << endl;
  cout << "z_size " << z_size << endl;
  cout << "t_size " << t_size << endl;

  data<MATRIX> conf;

  // read configuration
  get_data(conf, path_conf, conf_format, bytes_skip, convert);

  cout.precision(17);

  std::ofstream output_correlator(path_output_correlator);
  output_correlator.precision(17);

  output_correlator << "distance,correlator" << endl;

  std::vector<double> polyakov_correlator_vec;

  vector<vector<MATRIX>> conf_separated = separate_wilson(conf.array);
  conf.array.clear();
  conf.array.shrink_to_fit();

  start_time = omp_get_wtime();

  if (correlator_type == "singlet") {
    polyakov_correlator_vec =
        polyakov_loop_correlator_singlet(conf_separated, D_max);
  } else if (correlator_type == "color_average") {
    polyakov_correlator_vec = polyakov_loop_correlator(conf_separated, D_max);
  } else {
    cout << "invalid correlator_type" << endl;
  }

  std ::map<double, double> polyakov_correlator =
      polyakov_average_directions(polyakov_correlator_vec, D_max);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "polyakov_loop_correlator time: " << search_time << std::endl;

  for (auto it = polyakov_correlator.begin(); it != polyakov_correlator.end();
       it++) {
    output_correlator << it->first << "," << it->second << std::endl;
  }

  output_correlator.close();
}
