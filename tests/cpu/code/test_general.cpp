#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/flux_tube.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/mag.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/smearing.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <numeric>
#include <stdio.h>
#include <tuple>
#include <unordered_map>
#include <vector>

int x_size;
int y_size;
int z_size;
int t_size;

#define MATRIX_TYPE su3

using namespace std;

int main(int argc, char *argv[]) {
  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  x_size = 36;
  y_size = 36;
  z_size = 36;
  t_size = 36;

  std::cout.precision(17);

  link1 link(x_size, y_size, z_size, t_size);
  data<MATRIX_TYPE> conf;
  // string conf_path = "../../confs/decomposed/test/monopoless";
  // string conf_path =
  // "../../confs/decomposed/monopoless/qc2dstag/40^4/mu0.05/"
  //                    "s0/conf_monopoless_0201";
  // string conf_path = "../../confs/su2/mag/su2_suzuki/24^4/beta2.4/conf_0001";
  // string conf_path =
  //     "../../confs/su2/monopole/su2_suzuki/24^4/beta2.4/CON_MON_MAG_001.LAT";
  // string conf_path = "../../../apps/monopole_decomposition_su2/test/result/"
  //                    "monopoless_24^4";
  // string conf_path =
  //     "../../confs/su2/monopole/su2_suzuki/24^4/beta2.4/CON_MON_MAG_001.LAT";
  // string conf_path =
  //     "../../confs/su2/monopoless/su2_suzuki/24^4/beta2.4/conf_monopoless_0001";
  // string conf_path =
  //     "../../confs/su2/monopole/su2_suzuki/24^4/beta2.4/conf_monopole_0001";
  // string conf_path =
  //     "../../../apps/monopole_decomposition_su2/test/result/monopole_40^4";
  // string conf_path =
  //     "../../confs/su2/monopole/qc2dstag/40^4/mu0.00/conf_monopole_0201";
  // string conf_path =
  //     "../../confs/su2/monopole/su2_suzuki/48^4/beta2.8/conf_monopole_0035";
  // string conf_path =
  //     "../../confs/su2/monopoless/su2_suzuki/48^4/beta2.8/conf_monopoless_0035";
  // string conf_path = "../../confs/su3/Landau_U1/gluodynamics/24^4/beta6.0/"
  //                    "conf_Landau_gaugefixed_0001";
  // string conf_path = "../../confs/su3/gluodynamics/32^4/beta6.2/CONF0001";
  // string conf_path =
  //     "../../confs/su3/mag/gluodynamics/32^4/beta6.2/CONFDP_gaugefixed_0001";
  // string conf_path =
  // "../../confs/su3/Landau_U1xU1/gluodynamics/32^4/beta6.2/"
  //                    "conf_Landau_gaugefixed_0001";
  string conf_path = "../../confs/su3/gluodynamics/36^4/beta6.3/CONF0001";
  // conf.read_double(conf_path, 0);
  conf.read_double_qc2dstag(conf_path);
  // conf.read_ildg(conf_path);
  // conf.read_float(conf_path, 8);
  // conf.read_double_convert_abelian(conf_path, 8);
  // conf.read_double_qc2dstag_convert_abelian(conf_path);

  // cout << "initial functional " << MAG_functional_su2(conf.array) << endl;

  // plakets and polyakov loop
  start_time = clock();
  std::cout << "qc2dstag plaket " << plaket(conf.array) << std::endl;
  std::cout << "qc2dstag plaket_time " << plaket_time(conf.array) << std::endl;
  std::cout << "qc2dstag plaket_space " << plaket_space(conf.array)
            << std::endl;
  std::cout << "qc2dstag polyakov " << polyakov(conf.array) << std::endl;
  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "plaket and staff time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;

  // cout << "MAG functional " << MAG_functional_su2(conf.array) << endl;

  std::vector<std::vector<MATRIX_TYPE>> conf_separated =
      separate_wilson(conf.array);
  // std::vector<std::vector<MATRIX_TYPE>> conf_separated =
  // separate_3(conf.array);

  cout << "plaket parallel " << plaket_parallel(conf_separated) << endl;
  cout << "plaket time parallel " << plaket_time_parallel(conf_separated)
       << endl;
  cout << "plaket space parallel " << plaket_space_parallel(conf_separated)
       << endl;

  // on-axis wilson loops
  int T_min = 1, T_max = 8;
  int R_min = 1, R_max = 8;

  std::vector<double> vec_wilson;
  start_time = clock();

  vec_wilson = wilson(conf.array, R_min, R_max, T_min, T_max);

  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "on-axis wilson time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;
  std::cout << "wilson_loops:" << std::endl;
  for (int T = T_min; T <= T_max; T++) {
    for (int R = R_min; R <= R_max; R++) {
      std::cout << "T = " << T << " R = " << R << " "
                << vec_wilson[(R - R_min) + (T - T_min) * (R_max - R_min + 1)]
                << std::endl;
    }
  }

  // off-axis wilson loops
  std::vector<std::vector<int>> directions;
  directions = generate_directions(4);

  start_time = clock();

  std::vector<wilson_result> wilson_offaxis_result =
      wilson_offaxis(conf.array, directions, 0.9, 4, 1, 4);

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
