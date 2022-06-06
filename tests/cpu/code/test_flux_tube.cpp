#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/flux_tube.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <omp.h>
#include <stdio.h>
#include <unordered_map>
#include <vector>

int x_size;
int y_size;
int z_size;
int t_size;

#define MATRIX_TYPE su2

int main(int argc, char *argv[]) {
  x_size = 40;
  y_size = 40;
  z_size = 40;
  t_size = 40;
  double start_time;
  double end_time;
  double search_time;
  link1 link(x_size, y_size, z_size, t_size);
  data<MATRIX_TYPE> conf_qc2dstag;
  data<MATRIX_TYPE> conf_smeared;
  std::string path_qc2dstag = "../../confs/qc2dstag/mu0.05/s0/CONF0201";
  std::string path_smeared =
      "../../confs/qc2dstag/HYP_APE/mu0.05/s0/smeared_0201";

  conf_qc2dstag.read_double_qc2dstag(path_qc2dstag);
  // conf_smeared.read_double(path_smeared, 0);
  conf_smeared.read_double_qc2dstag(path_qc2dstag);

  // std::cout << "schwinger wilson flux tube test" << std::endl;

  int R = 10;
  int T = 10;

  // start_time = omp_get_wtime();
  // std::vector<MATRIX_TYPE> plaket_electric =
  //     calculate_plaket_schwinger_time(conf_qc2dstag.array);
  // end_time = clock();
  // search_time = end_time - start_time;
  // std::cout << "calculate_plaket_schwinger_time time: "
  //           << search_time * 1. / CLOCKS_PER_SEC << std::endl;

  // start_time = omp_get_wtime();
  // std::map<int, double> schwinger_electric =
  // wilson_plaket_schwinger_electric(
  //     conf_qc2dstag.array, plaket_electric, -5, R + 5, T, R);
  // end_time = clock();
  // search_time = end_time - start_time;
  // std::cout << "wilson_plaket_schwinger_electric time: "
  //           << search_time * 1. / CLOCKS_PER_SEC << std::endl;

  // for (auto it = schwinger_electric.begin(); it != schwinger_electric.end();
  //      ++it) {
  //   std::cout << "d = " << it->first << " schwinger_electric = " <<
  //   it->second
  //             << std::endl;
  // }

  std::cout << "wilson plaket flux tube test" << std::endl;

  R = 4;
  T = 4;
  int d_min = -5;
  int d_max = R + 5;

  start_time = omp_get_wtime();
  std::vector<double> wilson_vec =
      calculate_wilson_loop_tr(conf_smeared.array, R, T);
  std::vector<double> plaket_time_vec =
      calculate_plaket_time_tr(conf_qc2dstag.array);
  std::vector<double> plaket_space_vec =
      calculate_plaket_space_tr(conf_qc2dstag.array);
  double wilson_loop = 0;
  double plaket = 0;
  for (int i = 0; i < wilson_vec.size(); i++) {
    wilson_loop += wilson_vec[i];
  }
  for (int i = 0; i < plaket_time_vec.size(); i++) {
    plaket += plaket_time_vec[i];
  }

  std::cout << "wilson_loop " << wilson_loop / wilson_vec.size() << std::endl;
  std::cout << "plaket " << plaket / plaket_time_vec.size() << std::endl;
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "prepare time: " << search_time << std::endl;

  start_time = omp_get_wtime();
  std::map<int, double> electric = wilson_plaket_correlator_electric(
      wilson_vec, plaket_time_vec, R, T, 0, d_min, d_max);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "longitudinal  electric time: " << search_time << std::endl;

  start_time = omp_get_wtime();
  std::map<int, double> magnetic = wilson_plaket_correlator_magnetic(
      wilson_vec, plaket_space_vec, R, T, 0, d_min, d_max);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "longitudinal magnetic time: " << search_time << std::endl;

  for (auto it = electric.begin(); it != electric.end(); ++it) {
    std::cout << "longitudinal electric R = " << R << " T = " << T << " d "
              << it->first << " " << it->second << std::endl;
  }
  std::cout << std::endl;
  for (auto it = magnetic.begin(); it != magnetic.end(); ++it) {
    std::cout << "longitudinal magnetic R = " << R << " T = " << T << " d "
              << it->first << " " << it->second << std::endl;
  }

  start_time = omp_get_wtime();
  std::map<int, double> electric_trans = wilson_plaket_correlator_electric_x(
      wilson_vec, plaket_time_vec, R, T, 5, R / 2);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "transversal electric time: " << search_time << std::endl;

  start_time = omp_get_wtime();
  std::map<int, double> magnetic_trans = wilson_plaket_correlator_magnetic_x(
      wilson_vec, plaket_space_vec, R, T, 5, R / 2);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "transversal magnetic time: " << search_time << std::endl;

  for (auto it = electric_trans.begin(); it != electric_trans.end(); ++it) {
    std::cout << "transversal electric R = " << R << " T = " << T << " d "
              << it->first << " " << it->second << std::endl;
  }
  std::cout << std::endl;
  for (auto it = magnetic_trans.begin(); it != magnetic_trans.end(); ++it) {
    std::cout << "transversal magnetic R = " << R << " T = " << T << " d "
              << it->first << " " << it->second << std::endl;
  }

  std::vector<std::vector<MATRIX_TYPE>> separated_unchanged =
      separate_wilson(conf_qc2dstag.array);

  int T_min = 4, T_max = 4;
  int R_min = 4, R_max = 4;

  std::map<std::tuple<int, int, int>, double> flux_tube_new;

  std::vector<double> plaket_time_tr = plaket_aver_tr_time(separated_unchanged);

  std::vector<double> plaket_space_tr =
      plaket_aver_tr_space(separated_unchanged);

  start_time = omp_get_wtime();

  flux_tube_new =
      wilson_plaket_correlator(plaket_time_tr, separated_unchanged, T_min,
                               T_max, R_min, R_max, 5, 0, "longitudinal");

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "flux tube longitudinal electric new time: " << search_time
            << std::endl;

  for (auto it = flux_tube_new.begin(); it != flux_tube_new.end(); it++) {
    std::cout << "T = " << std::get<0>(it->first)
              << ", R = " << std::get<1>(it->first)
              << ", d = " << std::get<2>(it->first)
              << ", longitudinal electric flux_tube = " << it->second
              << std::endl;
  }

  start_time = omp_get_wtime();

  flux_tube_new =
      wilson_plaket_correlator(plaket_space_tr, separated_unchanged, T_min,
                               T_max, R_min, R_max, 5, 0, "longitudinal");

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "flux tube longitudinal magnetic new time: " << search_time
            << std::endl;

  for (auto it = flux_tube_new.begin(); it != flux_tube_new.end(); it++) {
    std::cout << "T = " << std::get<0>(it->first)
              << ", R = " << std::get<1>(it->first)
              << ", d = " << std::get<2>(it->first)
              << ", longitudinal magnetic flux_tube = " << it->second
              << std::endl;
  }

  start_time = omp_get_wtime();

  flux_tube_new =
      wilson_plaket_correlator(plaket_time_tr, separated_unchanged, T_min,
                               T_max, R_min, R_max, 5, R / 2, "transversal");

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "flux tube transversal electric new time: " << search_time
            << std::endl;

  for (auto it = flux_tube_new.begin(); it != flux_tube_new.end(); it++) {
    std::cout << "T = " << std::get<0>(it->first)
              << ", R = " << std::get<1>(it->first)
              << ", d = " << std::get<2>(it->first)
              << ", transversal electric flux_tube = " << it->second
              << std::endl;
  }

  start_time = omp_get_wtime();

  flux_tube_new =
      wilson_plaket_correlator(plaket_space_tr, separated_unchanged, T_min,
                               T_max, R_min, R_max, 5, R / 2, "transversal");

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "flux tube transversal magnetic new time: " << search_time
            << std::endl;

  for (auto it = flux_tube_new.begin(); it != flux_tube_new.end(); it++) {
    std::cout << "T = " << std::get<0>(it->first)
              << ", R = " << std::get<1>(it->first)
              << ", d = " << std::get<2>(it->first)
              << ", transversal magnetic flux_tube = " << it->second
              << std::endl;
  }

  // std::vector<int> steps = {1, x_size, x_size * y_size,
  //                           x_size * y_size * z_size,
  //                           x_size * y_size * z_size * t_size};

  // std::vector<double> wilson_tr;
  // std::vector<MATRIX_TYPE> space_lines;
  // std::vector<MATRIX_TYPE> time_lines;

  // std::vector<std::vector<MATRIX_TYPE>> separated_unchanged1 =
  //     separate_wilson(conf_qc2dstag.array);

  // std::vector<double> plaket_tr = plaket_aver_tr_time(separated_unchanged1);

  // // start_time = omp_get_wtime();

  // space_lines = wilson_lines(separated_unchanged1[0], 10, steps[0],
  // steps[1]); time_lines = wilson_lines(separated_unchanged1[3], 10, steps[3],
  // steps[4]);

  // // end_time = omp_get_wtime();
  // // search_time = end_time - start_time;
  // // std::cout << "test time: " << search_time << std::endl;

  // wilson_tr = wilson_plane_tr(space_lines, time_lines, steps[0], steps[1],
  //                             steps[3], steps[4], 10, 10);

  // int main_coordinate_min = -5;
  // int main_coordinate_max = 15;

  // std::vector<double> correlator(main_coordinate_max - main_coordinate_min +
  // 1,
  //                                0.0);

  // start_time = omp_get_wtime();

  // wilson_plaket_correlator_plane_longitudinal(
  //     correlator, wilson_tr, plaket_tr, steps[0], steps[1], 0,
  //     main_coordinate_min, main_coordinate_max);

  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // std::cout << "test time: " << search_time << std::endl;
}