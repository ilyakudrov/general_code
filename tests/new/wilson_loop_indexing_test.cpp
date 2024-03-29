#include "../../lib/cpu/include/basic_observables.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/flux_tube.h"
#include "../../lib/cpu/include/link.h"
#include "../../lib/cpu/include/matrix.h"
#include "../../lib/cpu/include/smearing.h"
#include "include/indexing.h"
#include "include/wilson_loop_indexing.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <stdio.h>
#include <tuple>
#include <unordered_map>
#include <vector>

int x_size;
int y_size;
int z_size;
int t_size;

#define MATRIX_TYPE abelian

using namespace std;

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double search_time;

  x_size = 48;
  y_size = 48;
  z_size = 48;
  t_size = 48;

  std::cout.precision(17);

  link1 link(x_size, y_size, z_size, t_size);
  // data<su2> conf;
  data<MATRIX_TYPE> conf;
  // data<abelian> conf;
  // string conf_path = "../confs/qc2dstag/40^4/mu0.00/CONF0201";
  string conf_path =
      "../confs/su2_suzuki/48^4/beta2.7/monopole/CON_MON_MAG_003.LAT";
  // string conf_path = "../confs/su2_suzuki/24^4/beta2.4/CON_fxd_MAG_001.LAT";
  // string conf_path = "../confs/SU3_conf/nt14/conf.0501";
  // conf.read_double(conf_path, 8);
  conf.read_float(conf_path, 4);
  // conf.read_double_qc2dstag(conf_path);
  // conf.read_ildg(conf_path);

  // int wilson_line_length = 10;
  // int wilson_lines_direction = 3;

  // start_time = clock();

  // std::vector<MATRIX_TYPE> wilson_line =
  //     wilson_lines(conf.array, wilson_lines_direction, wilson_line_length);

  // end_time = clock();
  // search_time = end_time - start_time;
  // std::cout << "wilson_lines time: " << search_time * 1. / CLOCKS_PER_SEC
  //           << std::endl;

  // std::cout << "wilson lines sum_tr "
  //           << std::accumulate(wilson_line.begin(), wilson_line.end(), 0.,
  //                              [](double a, MATRIX_TYPE B) { return a +
  //                              B.tr(); })
  //           << std::endl;

  // std::vector<std::vector<MATRIX_TYPE>> separated =
  // separate_wilson(conf.array);

  // int test_start = 94;
  // int test_end = 98;

  // start_time = clock();

  // wilson_line =
  //     wilson_lines_test1(separated[wilson_lines_direction],
  //     wilson_line_length,
  //                        link.lattice_size[wilson_lines_direction]);

  // end_time = clock();
  // search_time = end_time - start_time;
  // std::cout << "wilson_lines_test1 time: " << search_time * 1. /
  // CLOCKS_PER_SEC
  //           << std::endl;

  // std::cout << "wilson lines test1 sum_tr "
  //           << std::accumulate(wilson_line.begin(), wilson_line.end(), 0.,
  //                              [](double a, MATRIX_TYPE B) { return a +
  //                              B.tr(); })
  //           << std::endl;

  // start_time = clock();

  // wilson_line =
  //     wilson_lines_test2(separated[wilson_lines_direction],
  //     wilson_line_length,
  //                        link.lattice_size[wilson_lines_direction]);

  // end_time = clock();
  // search_time = end_time - start_time;
  // std::cout << "wilson_lines_test2 time: " << search_time * 1. /
  // CLOCKS_PER_SEC
  //           << std::endl;

  // std::cout << "wilson lines test2 sum_tr "
  //           << std::accumulate(wilson_line.begin(), wilson_line.end(), 0.,
  //                              [](double a, MATRIX_TYPE B) { return a +
  //                              B.tr(); })
  //           << std::endl;

  // std::vector<MATRIX_TYPE> wilson_line1 = wilson_lines_test2(
  //     separated[wilson_lines_direction], wilson_line_length + 1,
  //     link.lattice_size[wilson_lines_direction]);

  // std::cout << "wilson lines length + 1 sum_tr "
  //           << std::accumulate(wilson_line1.begin(), wilson_line1.end(), 0.,
  //                              [](double a, MATRIX_TYPE B) { return a +
  //                              B.tr(); })
  //           << std::endl;

  // start_time = clock();

  // wilson_line = wilson_lines_increase_test(
  //     separated[wilson_lines_direction], wilson_line, wilson_line_length + 1,
  //     link.lattice_size[wilson_lines_direction]);

  // end_time = clock();
  // search_time = end_time - start_time;
  // std::cout << "wilson_lines_increase_test time: "
  //           << search_time * 1. / CLOCKS_PER_SEC << std::endl;

  // std::cout << "wilson_lines_increase_test sum_tr "
  //           << std::accumulate(wilson_line.begin(), wilson_line.end(), 0.,
  //                              [](double a, MATRIX_TYPE B) { return a +
  //                              B.tr(); })
  //           << std::endl;

  // std::vector<std::vector<MATRIX_TYPE>> separated_unchanged =
  //     separate_wilson_unchanged(conf.array);

  // start_time = clock();

  // int size3 = 0;

  // if (wilson_lines_direction < 3)
  //   size3 = link.multiplier[wilson_lines_direction + 1] / 4;
  // else
  //   size3 = link.multiplier[wilson_lines_direction] / 4 *
  //   link.lattice_size[3];

  // wilson_line = wilson_lines_test3(
  //     separated_unchanged[wilson_lines_direction], wilson_line_length,
  //     link.multiplier[wilson_lines_direction] / 4, size3);

  // end_time = clock();
  // search_time = end_time - start_time;
  // std::cout << "wilson_lines_test3 time: " << search_time * 1. /
  // CLOCKS_PER_SEC
  //           << std::endl;

  // std::cout << "wilson_lines_test3 sum_tr "
  //           << std::accumulate(wilson_line.begin(), wilson_line.end(), 0.,
  //                              [](double a, MATRIX_TYPE B) { return a +
  //                              B.tr(); })
  //           << std::endl;

  int T_min = 1, T_max = 4;
  int R_min = 1, R_max = 4;

  std::vector<double> vec_wilson;
  start_time = omp_get_wtime();

  vec_wilson = wilson(conf.array, R_min, R_max, T_min, T_max);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "on-axis wilson time: " << search_time << std::endl;

  std::cout << "wilson_loops:" << std::endl;
  for (int T = T_min; T <= T_max; T++) {
    for (int R = R_min; R <= R_max; R++) {
      std::cout << "T = " << T << " R = " << R << " "
                << vec_wilson[(R - R_min) + (T - T_min) * (R_max - R_min + 1)]
                << std::endl;
    }
  }

  std::vector<std::vector<MATRIX_TYPE>> separated_unchanged =
      separate_wilson_unchanged(conf.array);

  int length_R = 20;
  int length_T = 10;

  std::vector<std::vector<MATRIX_TYPE>> wilson_lines(4);
  for (int i = 0; i < 3; i++) {
    wilson_lines[i] =
        wilson_lines_test3(separated_unchanged[i], length_R,
                           link.multiplier[i] / 4, link.multiplier[i + 1] / 4);
  }
  wilson_lines[3] = wilson_lines_test3(
      separated_unchanged[3], length_T, link.multiplier[3] / 4,
      link.multiplier[3] / 4 * link.lattice_size[3]);

  std::cout << wilson_loop_test_time(wilson_lines, length_R, length_T, 4)
            << std::endl;

  start_time = omp_get_wtime();

  std::map<std::tuple<int, int>, double> wilson_loops_new =
      wilson_parallel(separated_unchanged, R_min, R_max, T_min, T_max);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "wilson_parallel time: " << search_time << std::endl;

  for (auto it = wilson_loops_new.begin(); it != wilson_loops_new.end(); it++) {
    cout << "T = " << get<0>(it->first) << ", R = " << get<1>(it->first)
         << ", wilson loop new = " << it->second << endl;
  }
}
