// #include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/indexing.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/matrix_test.h"
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
  Data::data<su3_test> conf2;
  Data::data<su3_eigen> conf3;
  Data::data<Eigen::Matrix3cd> conf4;

  string conf_path1 = "../../confs/su3/QCD/140MeV/nt20/conf.0501";
  string conf_format1 = "ildg";
  int bytes_skip = 0;
  bool convert = 0;

  int R_min = 1;
  int R_max = 32;
  int T_min = 1;
  int T_max = 10;

  // get_data(conf1, conf_path1, conf_format1, bytes_skip, convert);
  // get_data(conf2, conf_path1, conf_format1, bytes_skip, convert);
  // get_data(conf3, conf_path1, conf_format1, bytes_skip, convert);
  get_data(conf4, conf_path1, conf_format1, bytes_skip, convert);

  map<tuple<int, int>, double> wilson_loops;

  // start_time = omp_get_wtime();
  // cout << "plaket indexed " << plaket_indexed(conf1.array) <<
  // endl; end_time = omp_get_wtime(); search_time = end_time -
  // start_time; std::cout << "plaket indexed time: " << search_time
  // << std::endl;

  // start_time = omp_get_wtime();
  // wilson_loops =
  //     wilson_parallel_indexed(conf1.array, R_min, R_max,
  //     T_min, T_max);
  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // std::cout << "wilson loop parallel indexed time: " <<
  // search_time
  //           << std::endl;
  // std::cout << "wilson_loops:" << std::endl;
  // for (auto it = wilson_loops.begin(); it !=
  // wilson_loops.end(); it++) {
  //   cout << get<0>(it->first) << "," << get<1>(it->first) <<
  //   "," << it->second
  //        << endl;
  // }

  // start_time = omp_get_wtime();
  // wilson_loops =
  // wilson_parallel_indexed_single_rxt(conf1.array, R_min,
  // R_max,
  //                                                   T_min,
  //                                                   T_max);
  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // std::cout << "wilson loop parallel indexed single_rxt
  // time: " << search_time
  //           << std::endl;
  // std::cout << "wilson_loops:" << std::endl;
  // for (auto it = wilson_loops.begin(); it !=
  // wilson_loops.end(); it++) {
  //   cout << get<0>(it->first) << "," << get<1>(it->first)
  //   << "," << it->second
  //        << endl;
  // }

  start_time = omp_get_wtime();
  wilson_loops = wilson_parallel_indexed_single_rxt(conf4.array, R_min, R_max,
                                                    T_min, T_max);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "wilson loop parallel indexed single_rxt su3_test time: "
            << search_time << std::endl;
  std::cout << "wilson_loops:" << std::endl;
  for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
    cout << get<0>(it->first) << "," << get<1>(it->first) << "," << it->second
         << endl;
  }

  // start_time = omp_get_wtime();
  // wilson_loops = wilson_parallel_indexed_single_rxt(conf4.array, R_min,
  // R_max,
  //                                                   T_min, T_max);
  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // std::cout << "wilson loop parallel indexed single_rxt su3_test time: "
  //           << search_time << std::endl;
  // std::cout << "wilson_loops:" << std::endl;
  // for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
  //   cout << get<0>(it->first) << "," << get<1>(it->first) << "," <<
  //   it->second
  //        << endl;
  // }

  // std::vector<std::vector<MATRIX_TYPE>> conf_separated =
  //     separate_wilson(conf1.array);
  // start_time = omp_get_wtime();
  // cout << "plaket parallel " << plaket_parallel(conf_separated) << endl;
  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // std::cout << "plaket parallel time: " << search_time << std::endl;

  // start_time = omp_get_wtime();
  // wilson_loops = wilson_parallel(conf_separated, R_min, R_max, T_min, T_max);
  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // std::cout << "wilson loop parallel time: " << search_time << std::endl;
  // std::cout << "wilson_loops:" << std::endl;
  // for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
  //   cout << get<0>(it->first) << "," << get<1>(it->first) << "," <<
  //   it->second
  //        << endl;
  // }
}