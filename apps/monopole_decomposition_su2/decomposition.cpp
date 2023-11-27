#include "../../lib/cpu/include/decomposition.h"
#include "../../lib/cpu/include/Landau_U1.h"
#include "../../lib/cpu/include/basic_observables.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/matrix.h"
#include "../../lib/cpu/include/monopoles.h"

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

void monopoles_test(vector<double> angles) {
  vector<double> J = calculate_current(angles);
  vector<loop *> LL = calculate_clusters(J);
  map<int, int> lengths_unwrapped;
  int length;
  vector<int> lengths_mu;

  for (int i = 0; i < LL.size(); i++) {
    length = cluster_length(LL[i]);
    lengths_mu = length_mu(LL[i]);
    if (lengths_mu[0] == 0 && lengths_mu[1] == 0 && lengths_mu[2] == 0 &&
        lengths_mu[3] == 0) {
      lengths_unwrapped[length]++;
    }
  }

  for (auto it = lengths_unwrapped.cbegin(); it != lengths_unwrapped.cend();
       ++it) {
    cout << it->first << "," << it->second << endl;
  }
}

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

  // read configuration
  bool convert = 0;
  get_data(conf_su2, path_conf, conf_format, bytes_skip, convert);

  cout.precision(17);

  double tolerance_maximal = 1e-5;
  double tolerance_average = 1e-7;
  int OR_steps = 4;

  std::vector<double> inverse_laplacian_real;
  std::vector<double> inverse_laplacian_imag;

  vector<complex_t> conf_complex = convert_to_complex(conf_su2.array);
  vector<complex_t> gauge_complex = generate_gauge_complex_uniform();
  // vector<complex_t> gauge_complex = generate_gauge_complex_unity();

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

  cout << "functional after applying gauge "
       << Landau_functional_complex(conf_complex) << endl;

  vector<double> conf_angles_U1 = convert_complex_to_angles(conf_complex);

  gauge_complex.clear();
  gauge_complex.shrink_to_fit();
  conf_complex.clear();
  conf_complex.shrink_to_fit();

  vector<double> inverse_laplacian =
      read_inverse_laplacian(path_inverse_laplacian);

  start_time = omp_get_wtime();

  cout << "decomposition started" << endl;

  vector<vector<int>> dirac_plakets =
      calculate_monopole_plaket_singular(conf_angles_U1);

  vector<double> monopole_angles =
      make_monopole_angles_parallel(dirac_plakets, inverse_laplacian);

  cout << "decomposition ended" << endl;

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "make_monopole_angles time: " << search_time << endl;

  cout << "monopoles for monopole part:" << endl;
  monopoles_test(monopole_angles);
  cout << endl;

  write_double_angles(path_conf_monopole, monopole_angles);

  get_monopoless_optimized(conf_su2.array, monopole_angles);

  vector<double> monopoless_angles = convert_to_angles(conf_su2.array);
  cout << "monopoles for monopoless part:" << endl;
  monopoles_test(monopoless_angles);

  write_double_su2(path_conf_monopoless, conf_su2.array);
}
