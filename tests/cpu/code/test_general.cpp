#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/flux_tube.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/result.h"
#include "../../../lib/cpu/include/smearing.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <stdio.h>
#include <tuple>
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

  x_size = 24;
  y_size = 24;
  z_size = 24;
  t_size = 24;

  std::cout.precision(17);

  link1 link(x_size, y_size, z_size, t_size);
  data<su2> conf_qc2dstag;
  // data<su3_full> conf_qc2dstag;
  // data<abelian> conf_qc2dstag;
  // std::string path_qc2dstag = "../../confs/qc2dstag/40^4/mu0.05/s0/CONF0201";
  // std::string path_qc2dstag = "../../confs/SU3_conf/nt6/conf.0501";
  // std::string path_qc2dstag = "../../confs/SU3_conf/nt6/conf.0502";
  // std::string path_qc2dstag =
  //     "../../confs/SU3_conf/nt6/conf.single_gaugefixed_0501.ildg_2";
  std::string path_qc2dstag =
      "../../confs/su2_suzuki/monopoless/CON_OFF_MAG_001.LAT";
  // std::string path_qc2dstag =
  // "../../confs/su2_dynam/32^4/CON_32^3x32_001.LAT";

  // conf_qc2dstag.read_double_qc2dstag(path_qc2dstag);
  // conf_qc2dstag.read_double(path_qc2dstag);
  conf_qc2dstag.read_double_fortran(path_qc2dstag, 8);

  for (int i = 0; i < 10; i++) {
    std::cout << conf_qc2dstag.array[i].module() << std::endl;
  }
  // conf_qc2dstag.read_float_fortran(path_qc2dstag);

  // std::string path_out_test =
  //     "../../confs/SU3_conf/nt6/conf.gaugefixed_0501_test";
  // conf_qc2dstag.write_double(path_out_test);

  // for (int i = 0; i < 10; i++) {
  //   std::cout << conf_qc2dstag.array[i] << std::endl;
  // }

  // polyakov correlator
  /*std::map<int, FLOAT> polyakov_correlator =
      polyakov_loop_correlator(conf_qc2dstag.array, 4, 16);

  for (auto it = polyakov_correlator.begin(); it != polyakov_correlator.end();
       ++it) {
    std::cout << "distance: " << it->first
              << " polyakov_correlator: " << it->second << std::endl;
  }
  std::cout << std::endl;

  polyakov_correlator.clear();

  double alpha_APE = 0.7;

  std::map<std::tuple<int, int>, std::vector<su2>> APE_2d =
      smearing_APE_2d(conf_qc2dstag.array, alpha_APE);

  for (int i = 0; i < 10; i++) {
    smearing_APE_2d_continue(APE_2d, alpha_APE);
  }

  // spatial wilson_lines
  std::map<std::tuple<int, int>, FLOAT> wilson_spat =
      wilson_spatial(conf_qc2dstag.array, APE_2d, 6, 8, 4, 8);

  for (auto it = wilson_spat.begin(); it != wilson_spat.end(); ++it) {
    std::cout << "distance: (" << std::get<0>(it->first) << ", "
              << std::get<1>(it->first) << ")"
              << " wilson_spatial: " << it->second << std::endl;
  }
  std::cout << std::endl;*/

  // plakets and polyakov loop
  start_time = clock();
  std::cout << "qc2dstag plaket " << 1 - plaket(conf_qc2dstag.array) / 2
            << std::endl;
  std::cout << "qc2dstag plaket_time "
            << 1 - plaket_time(conf_qc2dstag.array) / 2 << std::endl;
  std::cout << "qc2dstag plaket_space "
            << 1 - plaket_space(conf_qc2dstag.array) / 2 << std::endl;
  std::cout << "qc2dstag polyakov " << polyakov(conf_qc2dstag.array) / 3
            << std::endl;
  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "plaket and staff time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;

  // on-axis wilson loops
  int T_min = 1, T_max = 3;
  int R_min = 1, R_max = 3;

  std::vector<FLOAT> vec_wilson;
  start_time = clock();

  vec_wilson = wilson(conf_qc2dstag.array, R_min, R_max, T_min, T_max);

  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "on-axis wilson time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;
  std::cout << "wilson_loops:" << std::endl;
  for (int T = T_min; T <= T_max; T++) {
    for (int R = R_min; R <= R_max; R++) {
      std::cout << "T = " << T << " R = " << R << " "
                << vec_wilson[(R - R_min) + (T - T_min) * (R_max - R_min + 1)] /
                       2
                << std::endl;
    }
  }

  // off-axis wilson loops
  std::vector<std::vector<int>> directions;
  directions = generate_directions(4);

  start_time = clock();

  std::vector<wilson_result> wilson_offaxis_result =
      wilson_offaxis(conf_qc2dstag.array, directions, 2, 3, 2, 3);

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
}
