#include "../../lib/cpu/include/abelian_projection_su3.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/decomposition.h"
#include "../../lib/cpu/include/matrix.h"

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

  double start_time;
  double end_time;
  double search_time;

  string path_conf;
  string conf_format;
  int bytes_skip = 0;
  string path_conf_monopole;
  string path_conf_monopoless;

  bool parallel = false;

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
  cout << "path_conf_monopole " << path_conf_monopole << endl;
  cout << "path_conf_monopoless " << path_conf_monopoless << endl;
  cout << "bytes_skip " << bytes_skip << endl;

  cout << "x_size " << x_size << endl;
  cout << "y_size " << y_size << endl;
  cout << "z_size " << z_size << endl;
  cout << "t_size " << t_size << endl;
  cout << "ml5_conf_num " << ml5_conf_num << endl;

  size1 = x_size * y_size;
  size2 = x_size * y_size * z_size;

  Data::data<su3> conf_su3;

  // read configuration
  if (std::string(conf_format) == "float") {
    conf_su3.read_float(path_conf, bytes_skip);
  } else if (std::string(conf_format) == "double") {
    conf_su3.read_double(path_conf, bytes_skip);
  } else if (std::string(conf_format) == "double_qc2dstag") {
    conf_su3.read_double_qc2dstag(path_conf);
  } else if (std::string(conf_format) == "ildg") {
    conf_su3.read_ildg(path_conf);
  } else {
    cout << "wrong conf format: " << conf_format << endl;
    return 0;
  }

  cout.precision(17);

  std::vector<std::vector<double>> angles_monopoole =
      read_double_angles_su3(path_conf_monopole);

  make_unitary(angles_monopoole);

  write_double_angles_su3(path_conf_monopole, angles_monopoole);

  get_monopoless_optimized_su3(conf_su3.array, angles_monopoole);

  write_double_su3(path_conf_monopoless, conf_su3.array);
}
