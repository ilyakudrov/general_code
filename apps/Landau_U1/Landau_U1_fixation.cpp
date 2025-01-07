#include "../../lib/cpu/include/Landau_U1.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/decomposition.h"
#include "../../lib/cpu/include/matrix.h"

#include <fstream>
#include <iostream>
#include <omp.h>

using namespace std;

// global variables for lattice size
int x_size;
int y_size;
int z_size;
int t_size;

int main(int argc, char **argv) {

  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  string path_conf;
  string conf_format;
  int bytes_skip = 0;
  string path_conf_output;
  string path_functional_output;

  // read parameters
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "-path_conf") {
      path_conf = argv[++i];
    } else if (string(argv[i]) == "-path_conf_output") {
      path_conf_output = argv[++i];
    } else if (string(argv[i]) == "-path_functional_output") {
      path_functional_output = argv[++i];
    } else if (string(argv[i]) == "-bytes_skip") {
      bytes_skip = stoi(string(argv[++i]));
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
  cout << "path_conf_output " << path_conf_output << endl;
  cout << "path_functional_output " << path_functional_output << endl;
  cout << "bytes_skip " << bytes_skip << endl;
  cout << "x_size " << x_size << endl;
  cout << "y_size " << y_size << endl;
  cout << "z_size " << z_size << endl;
  cout << "t_size " << t_size << endl;

  Data::data<su2> conf_su2;

  // read configuration
  get_data(conf_su2, path_conf, conf_format, bytes_skip, 0);

  double tolerance_maximal = 1e-8;
  double tolerance_average = 1e-12;
  int OR_steps = 4;

  cout.precision(17);

  vector<std::complex<double>> conf_complex =
      convert_to_complex(conf_su2.array);
  vector<std::complex<double>> gauge_complex = generate_gauge_complex_uniform();

  cout << "initial Landau U1 functional "
       << Landau_functional_complex(conf_complex) << endl;

  start_time = omp_get_wtime();

  make_maximization_final(gauge_complex, conf_complex, OR_steps,
                          tolerance_maximal, tolerance_average);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "make_maximization_final time: " << search_time << endl;

  normalize_complex(gauge_complex);

  apply_gauge_Landau_complex(gauge_complex, conf_complex);
  apply_gauge_Landau(gauge_complex, conf_su2.array);

  double Landau_functional = Landau_functional_complex(conf_complex);

  cout << "Landau functional " << Landau_functional << endl;

  ofstream functional_stream;
  functional_stream.open(path_functional_output);

  functional_stream.precision(17);

  functional_stream << "functional" << endl;
  functional_stream << Landau_functional << endl;

  write_double_su2(path_conf_output, conf_su2.array);
}
