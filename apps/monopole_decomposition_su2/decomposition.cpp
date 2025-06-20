#include "../../lib/cpu/include/decomposition.h"
#include "../../lib/cpu/include/Landau_U1.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/matrix.h"
#include "../../lib/cpu/include/monopoles.h"

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

  int x_size1;
  int y_size1;
  int z_size1;
  int t_size1;

  string path_conf;
  string conf_format;
  string file_precision;
  int bytes_skip = 0;
  string path_conf_monopole;
  string path_conf_monopoless;
  string path_inverse_laplacian;

  // read parameters
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "--conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "--file_precision") {
      file_precision = argv[++i];
    } else if (string(argv[i]) == "--path_conf") {
      path_conf = argv[++i];
    } else if (string(argv[i]) == "--path_conf_monopole") {
      path_conf_monopole = argv[++i];
    } else if (string(argv[i]) == "--path_conf_monopoless") {
      path_conf_monopoless = argv[++i];
    } else if (string(argv[i]) == "--path_inverse_laplacian") {
      path_inverse_laplacian = argv[++i];
    } else if (string(argv[i]) == "--bytes_skip") {
      bytes_skip = stoi(string(argv[++i]));
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
  cout << "path_conf_monopole " << path_conf_monopole << endl;
  cout << "path_conf_monopoless " << path_conf_monopoless << endl;
  cout << "path_inverse_laplacian " << path_inverse_laplacian << endl;
  cout << "bytes_skip " << bytes_skip << endl;

  cout << "x_size " << x_size1 << endl;
  cout << "y_size " << y_size1 << endl;
  cout << "z_size " << z_size1 << endl;
  cout << "t_size " << t_size1 << endl;

  // read configuration
  bool convert = 0;
  Data::LatticeData<DataPatternLexicographical, su2> conf_su2(
      {x_size1, y_size1, z_size1, t_size1});
  Data::read_data_convert(conf_su2, path_conf, conf_format, bytes_skip,
                          file_precision, convert);

  DataPatternLexicographical data_pattern_conf(conf_su2.lat_dim);
  DataPatternLexicographical data_pattern_laplacian(
      {conf_su2.lat_dim[0] / 2 + 1, conf_su2.lat_dim[1] / 2 + 1,
       conf_su2.lat_dim[2] / 2 + 1, conf_su2.lat_dim[3] / 2 + 1});

  cout.precision(17);

  double tolerance_maximal = 1e-9;
  double tolerance_average = 1e-13;
  int OR_steps = 4;

  vector<std::complex<double>> conf_complex = convert_to_complex(conf_su2);
  vector<std::complex<double>> gauge_complex =
      generate_gauge_complex_uniform(data_pattern_conf);

  cout << "initial Landau U1 functional "
       << Landau_functional_complex(conf_complex) << endl;

  start_time = omp_get_wtime();

  make_maximization_final(gauge_complex, conf_complex, OR_steps,
                          tolerance_maximal, tolerance_average);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "make_maximization_final time: " << search_time << endl;

  normalize_complex(gauge_complex);
  apply_gauge_Landau_complex(gauge_complex, conf_complex, data_pattern_conf);
  apply_gauge_Landau(gauge_complex, conf_su2);

  cout << "functional after applying gauge "
       << Landau_functional_complex(conf_complex) << endl;

  vector<double> conf_angles_U1 = convert_complex_to_angles(conf_complex);

  gauge_complex.clear();
  gauge_complex.shrink_to_fit();
  conf_complex.clear();
  conf_complex.shrink_to_fit();

  vector<double> inverse_laplacian =
      read_inverse_laplacian(path_inverse_laplacian, data_pattern_laplacian);

  start_time = omp_get_wtime();

  cout << "decomposition started" << endl;

  vector<vector<int>> dirac_plakets =
      calculate_monopole_plaket_singular(conf_angles_U1, data_pattern_conf);

  vector<double> monopole_angles =
      make_monopole_angles(dirac_plakets, inverse_laplacian, data_pattern_conf,
                           data_pattern_laplacian);

  cout << "decomposition ended" << endl;

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "make_monopole_angles time: " << search_time << endl;

  write_double_angles(path_conf_monopole, monopole_angles, data_pattern_conf);
  get_monopoless_optimized(conf_su2.array, monopole_angles);
  FilePatternLexicographical<4, su2> file_pattern_lexicographical;
  conf_su2.write_data(path_conf_monopoless, file_pattern_lexicographical);
}
