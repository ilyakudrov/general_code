#include "../../lib/cpu/include/Landau_U1.h"
#include "../../lib/cpu/include/abelian_projection_su3.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/decomposition.h"
#include "../../lib/cpu/include/matrix.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <omp.h>

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
  string path_conf_monopole;
  string path_conf_monopoless;
  string path_inverse_laplacian;

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
    } else if (string(argv[i]) == "-parallel") {
      istringstream(string(argv[++i])) >> parallel;
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

  data<su3> conf_su3;

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

  // std::vector<std::vector<double>> angles_su3 =
  // get_angles_su3(conf_su3.array);
  std::vector<std::vector<double>> angles_su3 = make_angles_SU3(conf_su3.array);

  vector<double> inverse_laplacian =
      read_inverse_laplacian(path_inverse_laplacian);

  double cos_sum = 0;
  for (int i = 0; i < angles_su3[0].size(); i++) {
    cos_sum += cos(angles_su3[0][i]);
  }
  cos_sum = cos_sum / angles_su3[0].size();

  cout << "cos sum = " << cos_sum << endl;

  vector<vector<double>> angles_monopoole(3);

  start_time = omp_get_wtime();

  for (int i = 0; i < 3; i++) {
    if (parallel) {
      angles_monopoole[i] =
          make_monopole_angles_parallel(angles_su3[i], inverse_laplacian);
    } else {
      angles_monopoole[i] =
          make_monopole_angles(angles_su3[i], inverse_laplacian);
    }
    angles_su3[i].erase(angles_su3[i].begin(), angles_su3[i].end());
  }

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "decomposition time: " << search_time << std::endl;

  cos_sum = 0;
  for (int i = 0; i < angles_monopoole[0].size(); i++) {
    cos_sum += cos(angles_monopoole[0][i]);
  }
  cos_sum = cos_sum / angles_monopoole[0].size();

  cout << "cos sum = " << cos_sum << endl;

  // start_time = clock();

  // cout << "decomposition started" << endl;

  // vector<double> monopole_angles =
  //     make_monopole_angles(conf_angles_U1, inverse_laplacian);

  // cout << "decomposition ended" << endl;

  // end_time = clock();
  // search_time = end_time - start_time;
  // cout << "make_monopole_angles time: " << search_time * 1. / CLOCKS_PER_SEC
  //      << endl;

  // write_double_angles(path_conf_monopole, monopole_angles);

  // get_monopoless_optimized(conf_su3.array, monopole_angles);

  // write_double_su2(path_conf_monopoless, conf_su3.array);
}
