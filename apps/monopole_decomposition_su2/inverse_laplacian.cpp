#include "../../lib/cpu/include/decomposition.h"

// #include <complex>
// #include <fstream>
#include <iostream>
// #include <map>
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
  double omp_time;

  int x_size1;
  int y_size1;
  int z_size1;
  int t_size1;

  string path_inverse_laplacian;

  // read parameters
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "--path_inverse_laplacian") {
      path_inverse_laplacian = argv[++i];
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

  cout << "path_inverse_laplacian " << path_inverse_laplacian << endl;
  cout << "x_size " << x_size1 << endl;
  cout << "y_size " << y_size1 << endl;
  cout << "z_size " << z_size1 << endl;
  cout << "t_size " << t_size1 << endl;

  // x_size = x_size1;
  // y_size = y_size1;
  // z_size = z_size1;
  // t_size = t_size1;

  DataPatternLexicographical data_pattern({x_size1, y_size1, z_size1, t_size1});
  DataPatternLexicographical data_pattern_laplacian(
      {data_pattern.lat_dim[0] / 2 + 1, data_pattern.lat_dim[1] / 2 + 1,
       data_pattern.lat_dim[2] / 2 + 1, data_pattern.lat_dim[3] / 2 + 1});
  cout.precision(17);
  omp_time = omp_get_wtime();
  std::vector<double> inverse_momenta = make_inverse_momenta(data_pattern);
  vector<double> inverse_laplacian =
      calculate_inverse_laplacian_spatially_symmetric(
          data_pattern, data_pattern_laplacian, inverse_momenta);
  std::cout << "inverse_laplacian time: " << omp_get_wtime() - omp_time
            << std::endl;
  write_inverse_laplacian(inverse_laplacian, path_inverse_laplacian,
                          data_pattern_laplacian);
}