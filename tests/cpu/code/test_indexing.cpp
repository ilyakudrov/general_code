#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/indexing.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/smearing.h"
#include "../../../lib/cpu/include/test_observables.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <omp.h>

int x_size;
int y_size;
int z_size;
int t_size;
int size1;
int size2;

#define MATRIX_TYPE su3

using namespace std;

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double search_time;

  x_size = 64;
  y_size = 64;
  z_size = 64;
  t_size = 20;
  size1 = x_size * y_size;
  size2 = x_size * y_size * z_size;

  std::cout.precision(17);

  Data::data<MATRIX_TYPE> conf1;
  Data::data<Eigen::Matrix3cd> conf4;

  string conf_path1 = "../../confs/su3/QCD/140MeV/nt20/conf.0501";
  string conf_format1 = "ildg";
  int bytes_skip = 0;
  bool convert = 0;

  int R_min = 1;
  int R_max = 32;
  int T_min = 1;
  int T_max = 2;

  map<tuple<int, int>, double> wilson_loops;

  get_data(conf1, conf_path1, conf_format1, bytes_skip, convert);

  start_time = omp_get_wtime();
  cout << "plaket indexed " << plaket_indexed(conf1.array) << endl;
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "plaket indexed su3 time: " << search_time << std::endl;

  start_time = omp_get_wtime();
  wilson_loops =
      wilson_parallel_indexed(conf1.array, R_min, R_max, T_min, T_max);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "wilson loop parallel indexed su3 time: " << search_time
            << std::endl;
  std::cout << "wilson_loops:" << std::endl;
  for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
    cout << get<0>(it->first) << "," << get<1>(it->first) << "," << it->second
         << endl;
  }

  start_time = omp_get_wtime();
  wilson_loops = wilson_parallel_indexed_single_rxt(conf1.array, R_min, R_max,
                                                    T_min, T_max);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "wilson loop parallel indexed single_rxt su3 time: "
            << search_time << std::endl;
  std::cout << "wilson_loops:" << std::endl;
  for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
    cout << get<0>(it->first) << "," << get<1>(it->first) << "," << it->second
         << endl;
  }

  std::vector<std::vector<MATRIX_TYPE>> conf_separated =
      separate_wilson(conf1.array);
  conf1.array.clear();
  conf1.array.shrink_to_fit();
  start_time = omp_get_wtime();
  cout << "plaket parallel " << plaket_parallel(conf_separated) << endl;
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "plaket parallel su3 time: " << search_time << std::endl;

  start_time = omp_get_wtime();
  wilson_loops = wilson_parallel(conf_separated, R_min, R_max, T_min, T_max);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "wilson loop parallel su3 time: " << search_time << std::endl;
  std::cout << "wilson_loops:" << std::endl;
  for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
    cout << get<0>(it->first) << "," << get<1>(it->first) << "," << it->second
         << endl;
  }
  for (int i = 0; i < 4; i++) {
    conf_separated[i].clear();
    conf_separated[i].shrink_to_fit();
  }

  get_data(conf4, conf_path1, conf_format1, bytes_skip, convert);

  start_time = omp_get_wtime();
  wilson_loops = wilson_parallel_indexed_single_rxt(conf4.array, R_min, R_max,
                                                    T_min, T_max);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "wilson loop parallel indexed single_rxt Eigen::Matrix3cd time: "
            << search_time << std::endl;
  std::cout << "wilson_loops:" << std::endl;
  for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
    cout << get<0>(it->first) << "," << get<1>(it->first) << "," << it->second
         << endl;
  }
  conf4.array.clear();
  conf4.array.shrink_to_fit();
}