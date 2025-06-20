#include "../../lib/cpu/include/abelian_projection_su3.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/decomposition.h"
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
  cout << "file_precision " << file_precision << endl;
  cout << "path_conf_monopole " << path_conf_monopole << endl;
  cout << "path_conf_monopoless " << path_conf_monopoless << endl;
  cout << "path_inverse_laplacian " << path_inverse_laplacian << endl;
  cout << "bytes_skip " << bytes_skip << endl;
  cout << "x_size " << x_size1 << endl;
  cout << "y_size " << y_size1 << endl;
  cout << "z_size " << z_size1 << endl;
  cout << "t_size " << t_size1 << endl;

  Data::LatticeData<DataPatternLexicographical, su3> conf_su3(
      {x_size1, y_size1, z_size1, t_size1});
  bool convert_su3 = 0;
  Data::read_data_convert(conf_su3, path_conf, conf_format, bytes_skip,
                          file_precision, convert_su3);
  Data::LatticeData<DataPatternLexicographical, su3_angles> conf(
      {x_size1, y_size1, z_size1, t_size1});
  bool convert_angles = 1;
  Data::read_data_convert(conf, path_conf, conf_format, bytes_skip,
                          file_precision, convert_angles);
  cout.precision(17);
  DataPatternLexicographical data_pattern_conf(conf.lat_dim);
  DataPatternLexicographical data_pattern_laplacian(
      {conf.lat_dim[0] / 2 + 1, conf.lat_dim[1] / 2 + 1,
       conf.lat_dim[2] / 2 + 1, conf.lat_dim[3] / 2 + 1});

  vector<double> inverse_laplacian =
      read_inverse_laplacian(path_inverse_laplacian, data_pattern_laplacian);
  vector<vector<vector<double>>> monopole_plakets(3);
  vector<vector<vector<int>>> dirac_plakets(3);
  vector<vector<double>> angles_monopole(3);
  make_plakets_both(conf, monopole_plakets, dirac_plakets);

  start_time = omp_get_wtime();
  for (int i = 0; i < 3; i++) {
    angles_monopole[i] =
        make_monopole_angles(dirac_plakets[i], inverse_laplacian,
                             data_pattern_conf, data_pattern_laplacian);
    dirac_plakets[i].erase(dirac_plakets[i].begin(), dirac_plakets[i].end());
  }
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "decomposition time: " << search_time << std::endl;

  make_unitary(angles_monopole);
  write_double_angles_su3(path_conf_monopole, angles_monopole,
                          data_pattern_conf);
  get_monopoless_optimized_su3(conf_su3.array, angles_monopole);
  FilePatternLexicographical<4, su3> file_pattern_lexicographical;
  conf_su3.write_data(path_conf_monopoless, file_pattern_lexicographical);
}
