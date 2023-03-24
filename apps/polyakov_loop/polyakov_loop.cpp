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
  string path_output;
  string correlator_type;
  bool convert = 0;

  // read parameters
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "-path_conf") {
      path_conf = argv[++i];
    } else if (string(argv[i]) == "-path_output") {
      path_output = argv[++i];
    } else if (string(argv[i]) == "-correlator_type") {
      correlator_type = argv[++i];
    } else if (string(argv[i]) == "-bytes_skip") {
      bytes_skip = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-convert") {
      istringstream(string(argv[++i])) >> convert;
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
  cout << "convert " << convert << endl;

  cout << "path_output " << path_output << endl;

  cout << "x_size " << x_size << endl;
  cout << "y_size " << y_size << endl;
  cout << "z_size " << z_size << endl;
  cout << "t_size " << t_size << endl;

  data<MATRIX> conf;

  // read configuration
  get_data(conf, path_conf, conf_format, bytes_skip, convert);

  cout.precision(17);

  std::ofstream output_polyakov_loop(path_output);
  output_polyakov_loop.precision(17);

  output_polyakov_loop << "polyakov_loop" << endl;

  std::vector<double> polyakov_polyakov_loop_vec;

  start_time = omp_get_wtime();

  double ploop = polyakov_loop(conf.array);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "polyakov_loop time: " << search_time << std::endl;

  output_polyakov_loop << ploop << std::endl;

  output_polyakov_loop.close();
}