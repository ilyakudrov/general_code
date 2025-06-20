#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/matrix.h"
#include "../../lib/cpu/include/polyakov_loops.h"

#include <fstream>
#include <map>
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

int main(int argc, char **argv) {
  double start_time;
  double end_time;
  double search_time;

  int x_size1;
  int y_size1;
  int z_size1;
  int t_size1;
  string path_conf;
  string conf_format;
  string file_precision;
  int bytes_skip = 0;
  string path_output_correlator;
  string correlator_type;
  bool convert = 0;
  int D_max;

  // read parameters
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "--conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "--file_precision") {
      file_precision = argv[++i];
    } else if (string(argv[i]) == "--path_conf") {
      path_conf = argv[++i];
    } else if (string(argv[i]) == "--path_output_correlator") {
      path_output_correlator = argv[++i];
    } else if (string(argv[i]) == "--correlator_type") {
      correlator_type = argv[++i];
    } else if (string(argv[i]) == "--bytes_skip") {
      bytes_skip = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--convert") {
      istringstream(string(argv[++i])) >> convert;
    } else if (string(argv[i]) == "--D_max") {
      D_max = stoi(string(argv[++i]));
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
  cout << "file_precision " << file_precision << endl;
  cout << "conf_format " << conf_format << endl;
  cout << "bytes_skip " << bytes_skip << endl;
  cout << "D_max " << D_max << endl;
  cout << "convert " << convert << endl;
  cout << "path_output_correlator " << path_output_correlator << endl;
  cout << "x_size " << x_size1 << endl;
  cout << "y_size " << y_size1 << endl;
  cout << "z_size " << z_size1 << endl;
  cout << "t_size " << t_size1 << endl;

  Data::LatticeData<DataPatternLexicographical, MATRIX> conf(
      {x_size1, y_size1, z_size1, t_size1});
  Data::read_data_convert(conf, path_conf, conf_format, bytes_skip,
                          file_precision, convert);
  cout.precision(17);

  std::cout << polyakov_loop(conf) << std::endl;

  std::vector<double> polyakov_correlator_vec;

  start_time = omp_get_wtime();
  if (correlator_type == "singlet") {
    polyakov_correlator_vec = polyakov_loop_correlator_singlet(conf, D_max);
  } else if (correlator_type == "color_average") {
    polyakov_correlator_vec =
        polyakov_loop_correlator_color_average(conf, D_max);
  } else {
    cout << "invalid correlator_type" << endl;
  }
  std ::map<double, double> polyakov_correlator =
      polyakov_average_directions(polyakov_correlator_vec, D_max);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "polyakov_loop_correlator time: " << search_time << std::endl;

  std::ofstream output_correlator(path_output_correlator);
  output_correlator.precision(17);
  output_correlator << "distance,correlator" << endl;

  for (auto it = polyakov_correlator.begin(); it != polyakov_correlator.end();
       it++) {
    output_correlator << it->first << "," << it->second << std::endl;
  }

  output_correlator.close();
}
