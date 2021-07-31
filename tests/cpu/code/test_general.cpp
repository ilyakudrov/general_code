#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/flux_tube.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/result.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <stdio.h>
#include <unordered_map>
#include <vector>

int x_size;
int y_size;
int z_size;
int t_size;

int main(int argc, char *argv[]) {
  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  y_size = 32;
  z_size = 32;
  t_size = 32;
  x_size = 32;

  std::cout.precision(17);

  link1 link(x_size, y_size, z_size, t_size);
  // data<su2> conf_qc2dstag;
  data<abelian> conf_qc2dstag;
  //   std::string path_qc2dstag =
  //   "../../confs/qc2dstag/40^4/mu0.05/s0/CONF0201"; std::string path_qc2dstag
  //   = "/home/ilya/soft/lattice/mag/confs/fixated/su2/"
  //                          "qc2dstag/32^4/mu0.00/conf_fixed_0201";
  //   std::string path_qc2dstag =
  //   "/home/ilya/soft/lattice/general_code/tests/confs/"
  //                          "qc2dstag/32^4/mu0.00/CONF0201";
  // std::string path_qc2dstag =
  //     "/home/ilya/soft/lattice/general_code/tests/confs/"
  //     "decomposed/monopole/qc2dstag/40^4/mu0.05/s0/conf_monopole_0202";
  // std::string path_qc2dstag =
  // "../../confs/decomposed/monopole/qc2dstag/40^4/"
  //                             "mu0.05/s0/conf_monopole_0201";
  // std::string path_qc2dstag =
  // "../../confs/decomposed/monopole/qc2dstag/40^4/"
  //                             "mu0.05/s0/conf_monopole_0201";
  // std::string path_qc2dstag = "/home/ilya/soft/lattice/decomposition/test/"
  //                             "confs/32^4/CON_32^3x32_031.LAT";
  std::string path_qc2dstag =
      "../../confs/su2_dynam/monopole/CON_MON_MAG_031.LAT";
  // std::string path_qc2dstag = "/home/ilya/soft/lattice/decomposition/test/"
  //                             "confs/monopoless/32^4/CON_OFFD_031";
  // std::string path_qc2dstag = "../../confs/decomposed/monopole/"
  //                             "qc2dstag/40^4/mu0.05/s0/conf_monopole_0201";
  // std::string path_qc2dstag =
  //     "../../confs/smeared/qc2dstag/40^4/mu0.05/s0/conf_APE_alpha=0.7_0202";

  // conf_qc2dstag.read_double(path_qc2dstag);
  // conf_qc2dstag.read_double_fortran(path_qc2dstag);
  conf_qc2dstag.read_float_fortran(path_qc2dstag);

  std::vector<double> conf_abelian = read_abelian_fortran_float(path_qc2dstag);
  std::cout << "wilson abelian test " << wilson_abelian(conf_abelian, 1, 1)
            << std::endl;
  std::cout << "wilson abelian test " << wilson_abelian(conf_abelian, 1, 18)
            << std::endl;
  std::cout << "wilson abelian test " << wilson_abelian(conf_abelian, 18, 1)
            << std::endl;
  std::cout << "wilson abelian test " << wilson_abelian(conf_abelian, 18, 18)
            << std::endl;

  int T_min = 14, T_max = 18;
  int R_min = 14, R_max = 18;

  std::vector<FLOAT> vec_wilson;
  start_time = clock();

  vec_wilson = wilson(conf_qc2dstag.array, R_min, R_max, T_min, T_max);

  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "wilson time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;
  std::cout << "wilson_loops:" << std::endl;
  for (int T = T_min; T <= T_max; T++) {
    for (int R = R_min; R <= R_max; R++) {
      std::cout << "T = " << T << " R = " << R << " "
                << vec_wilson[(R - R_min) + (T - T_min) * (R_max - R_min + 1)]
                << std::endl;
    }
  }

  std::vector<std::vector<int>> directions;
  directions = generate_directions(4);

  start_time = clock();

  std::vector<wilson_result> wilson_offaxis_result =
      wilson_offaxis(conf_qc2dstag.array, directions, 14, 18, 14, 18);

  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "offaxis wilson loop time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;

  for (int i = 0; i < wilson_offaxis_result.size(); i++) {
    std::cout << "wilson loop: " << wilson_offaxis_result[i].wilson_loop
              << " time size: " << wilson_offaxis_result[i].time_size
              << " space size: " << wilson_offaxis_result[i].space_size
              << std::endl;
  }

  wilson_offaxis_reduce(wilson_offaxis_result);

  std::cout << std::endl << "after reduction" << std::endl;

  for (int i = 0; i < wilson_offaxis_result.size(); i++) {
    std::cout << "wilson loop: " << wilson_offaxis_result[i].wilson_loop
              << " time size: " << wilson_offaxis_result[i].time_size
              << " space size: " << wilson_offaxis_result[i].space_size
              << std::endl;
  }

  start_time = clock();
  // cout << "MAG_functional_su2 " << MAG_functional_su2(conf_qc2dstag.array)
  //      << endl;
  end_time = clock();
  search_time = end_time - start_time;
  std::cout << " time: " << search_time * 1. / CLOCKS_PER_SEC << std::endl;

  start_time = clock();
  std::cout << "qc2dstag plaket " << plaket(conf_qc2dstag.array) / 2
            << std::endl;
  std::cout << "qc2dstag plaket_time " << plaket_time(conf_qc2dstag.array) / 2
            << std::endl;
  std::cout << "qc2dstag plaket_space " << plaket_space(conf_qc2dstag.array) / 2
            << std::endl;
  std::cout << "qc2dstag polyakov " << polyakov(conf_qc2dstag.array) / 2
            << std::endl;
  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "plaket and staff time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;

  double start;
  double end;

  int R = 8;
  int T = 4;
  int d_min = -5;
  int d_max = 0;

  start_time = clock();
  std::vector<FLOAT> wilson_vec =
      calculate_wilson_loop_tr(conf_qc2dstag.array, R, T);
  std::vector<FLOAT> plaket_time_vec =
      calculate_plaket_time_tr(conf_qc2dstag.array);
  std::vector<FLOAT> plaket_space_vec =
      calculate_plaket_space_tr(conf_qc2dstag.array);
  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "prepare time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;

  start_time = clock();
  std::vector<FLOAT> electric = wilson_plaket_correlator_electric(
      wilson_vec, plaket_time_vec, R, T, 0, d_min, d_max);
  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "electric time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;

  start_time = clock();
  std::vector<FLOAT> magnetic = wilson_plaket_correlator_magnetic(
      wilson_vec, plaket_space_vec, R, T, 0, d_min, d_max);
  end_time = clock();
  search_time = end_time - start_time;

  std::cout << "magnetic time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;
  for (int i = 0; i < electric.size(); i++) {
    std::cout << "electric R = " << R << " T= " << T << " " << electric[i]
              << std::endl;
    std::cout << "magnetic R = " << R << " T= " << T << " " << magnetic[i]
              << std::endl;
  }

  int x_trans_min = -12;
  int x_trans_max = 12;

  start_time = clock();
  std::vector<FLOAT> electric_x = wilson_plaket_correlator_electric_x(
      wilson_vec, plaket_time_vec, R, T, x_trans_min, x_trans_max, R / 2);
  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "electric_x time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;

  start_time = clock();
  std::vector<FLOAT> magnetic_x = wilson_plaket_correlator_magnetic_x(
      wilson_vec, plaket_space_vec, R, T, x_trans_min, x_trans_max, R / 2);
  end_time = clock();
  search_time = end_time - start_time;

  std::cout << "magnetic_x time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;
  for (int i = 0; i < electric_x.size(); i++) {
    std::cout << "electric_x R = " << R << " T= " << T << " " << electric_x[i]
              << std::endl;
    std::cout << "magnetic_x R = " << R << " T= " << T << " " << magnetic_x[i]
              << std::endl;
  }
}
