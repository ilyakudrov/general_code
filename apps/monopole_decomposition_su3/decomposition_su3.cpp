#include "../../lib/cpu/include/abelian_projection_su3.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/decomposition.h"
#include "../../lib/cpu/include/matrix.h"
#include "../../lib/cpu/include/monopoles.h"

#include <iostream>
#include <omp.h>
#include <sstream>

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

  bool compensate_dirac = true;
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
    } else if (string(argv[i]) == "-compensate_dirac") {
      istringstream(string(argv[++i])) >> compensate_dirac;
    } else
      cout << "unknown parameter " << argv[i] << endl;
  }

  cout << "path_conf " << path_conf << endl;
  cout << "conf_format " << conf_format << endl;
  cout << "path_conf_monopole " << path_conf_monopole << endl;
  cout << "path_conf_monopoless " << path_conf_monopoless << endl;
  cout << "path_inverse_laplacian " << path_inverse_laplacian << endl;
  cout << "bytes_skip " << bytes_skip << endl;
  cout << "parallel " << parallel << endl;
  cout << "compensate_dirac " << compensate_dirac << endl;

  cout << "x_size " << x_size << endl;
  cout << "y_size " << y_size << endl;
  cout << "z_size " << z_size << endl;
  cout << "t_size " << t_size << endl;
  cout << "ml5_conf_num " << ml5_conf_num << endl;

  data<su3> conf_su3;

  // read configuration
  bool convert = 0;
  get_data(conf_su3, path_conf, conf_format, bytes_skip, convert);

  cout.precision(17);

  std::vector<std::vector<double>> angles_su3 = make_angles_SU3(conf_su3.array);

  vector<double> inverse_laplacian =
      read_inverse_laplacian(path_inverse_laplacian);

  vector<vector<vector<double>>> monopole_plakets(3);
  vector<vector<vector<int>>> dirac_plakets(3);

  if (compensate_dirac) {
    make_plakets_both(angles_su3, monopole_plakets, dirac_plakets);
    for (int i = 0; i < monopole_plakets.size(); i++) {
      monopole_plakets[i].erase(monopole_plakets[i].begin(),
                                monopole_plakets[i].end());
      angles_su3[i].erase(angles_su3[i].begin(), angles_su3[i].end());
    }
  } else {
    for (int i = 0; i < monopole_plakets.size(); i++) {
      dirac_plakets[i] = calculate_monopole_plaket_singular(angles_su3[i]);
      angles_su3[i].erase(angles_su3[i].begin(), angles_su3[i].end());
    }
  }

  vector<vector<double>> angles_monopole(3);

  start_time = omp_get_wtime();

  for (int i = 0; i < 3; i++) {
    if (parallel) {
      angles_monopole[i] =
          make_monopole_angles_parallel(dirac_plakets[i], inverse_laplacian);
    } else {
      angles_monopole[i] =
          make_monopole_angles(dirac_plakets[i], inverse_laplacian);
    }
    dirac_plakets[i].erase(dirac_plakets[i].begin(), dirac_plakets[i].end());
  }

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "decomposition time: " << search_time << std::endl;

  make_unitary(angles_monopole);

  write_double_angles_su3(path_conf_monopole, angles_monopole);

  get_monopoless_optimized_su3(conf_su3.array, angles_monopole);

  write_double_su3(path_conf_monopoless, conf_su3.array);
}
