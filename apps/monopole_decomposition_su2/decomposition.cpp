#include "../../lib/cpu/include/decomposition.h"
#include "../../lib/cpu/include/Landau_U1.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/matrix.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <map>

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
  string path_conf_monopole;
  string path_conf_monopoless;
  string path_inverse_laplacian;

  int ml5_conf_num = 0;

  // read parameters
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "-path_conf") {
      path_conf = argv[++i];
    } else if (string(argv[i]) == "-path_conf_monopole") {
      path_conf_monopole = argv[++i];
    } else if (string(argv[i]) == "-path_conf_monopoless") {
      path_conf_monopoless = argv[++i];
    } else if (string(argv[i]) == "-path_inverse_laplacian") {
      path_inverse_laplacian = argv[++i];
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
    } else if (string(argv[i]) == "-ml5_conf_num") {
      ml5_conf_num = stoi(string(argv[++i]));
    } else
      cout << "unknown parameter " << argv[i] << endl;
  }

  cout << "path_conf " << path_conf << endl;
  cout << "conf_format " << conf_format << endl;
  cout << "path_conf_monopole " << path_conf_monopole << endl;
  cout << "path_conf_monopoless " << path_conf_monopoless << endl;
  cout << "path_inverse_laplacian " << path_inverse_laplacian << endl;
  cout << "bytes_skip " << bytes_skip << endl;

  cout << "x_size " << x_size << endl;
  cout << "y_size " << y_size << endl;
  cout << "z_size " << z_size << endl;
  cout << "t_size " << t_size << endl;
  cout << "ml5_conf_num " << ml5_conf_num << endl;

  data<su2> conf_su2;

  // for ml5 configuration
  vector<float> ml5_data;

  // read configuration
  if (std::string(conf_format) == "float") {
    conf_su2.read_float(path_conf, bytes_skip);
  } else if (std::string(conf_format) == "double") {
    conf_su2.read_double(path_conf, bytes_skip);
  } else if (std::string(conf_format) == "double_qc2dstag") {
    conf_su2.read_double_qc2dstag(path_conf);
  } else if (conf_format == "ml5") {
    ml5_data = read_full_ml5(path_conf, ml5_conf_num);
    conf_su2.read_float_ml5(ml5_data, ml5_conf_num);
  } else {
    cout << "wrong conf format: " << conf_format << endl;
    return 0;
  }

  cout.precision(17);

  double tolerance_maximal = 1e-5;
  double tolerance_average = 1e-7;
  int OR_steps = 4;

  vector<complex_t> conf_complex = convert_to_complex(conf_su2.array);
  vector<complex_t> gauge_complex = generate_gauge_complex_uniform();

  cout << "initial Landau U1 functional "
       << Landau_functional_complex(conf_complex) << endl;

  start_time = clock();

  make_maximization_final(gauge_complex, conf_complex, OR_steps,
                          tolerance_maximal, tolerance_average);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "make_maximization_final time: " << search_time * 1. / CLOCKS_PER_SEC
       << endl;

  normalize_complex(gauge_complex);

  apply_gauge_Landau_complex(gauge_complex, conf_complex);
  apply_gauge_Landau(gauge_complex, conf_su2.array);

  cout << "functional after applying gauge "
       << Landau_functional_complex(conf_complex) << endl;

  vector<double> conf_angles_U1 = convert_complex_to_angles(conf_complex);

  gauge_complex.clear();
  gauge_complex.shrink_to_fit();
  conf_complex.clear();
  conf_complex.shrink_to_fit();

  vector<double> inverse_laplacian =
      read_inverse_laplacian(path_inverse_laplacian);

  start_time = clock();

  vector<double> monopole_angles =
      make_monopole_angles(conf_angles_U1, inverse_laplacian);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "make_monopole_angles time: " << search_time * 1. / CLOCKS_PER_SEC
       << endl;

  write_double_angles(path_conf_monopole, monopole_angles);

  data<su2> conf_monopoless;
  conf_monopoless.array = get_monopoless(conf_su2.array, monopole_angles);

  conf_monopoless.write_double(path_conf_monopoless);
}
