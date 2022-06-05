#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/flux_tube.h"
#include "../../../lib/cpu/include/link.h"
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

#define MATRIX_TYPE su2

using namespace std;

int main(int argc, char *argv[]) {
  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  x_size = 40;
  y_size = 40;
  z_size = 40;
  t_size = 40;

  std::cout.precision(17);

  link1 link(x_size, y_size, z_size, t_size);
  data<su2> conf;
  // data<su3_full> conf;
  // data<su2> conf;
  // data<abelian> conf;
  // string conf_path = "../../confs/su2_suzuki/monopoless/"
  //                    "HYP0_alpha=1_1_0.5_APE100_alpha=0.5/conf_0001";
  // string conf_path = "../../confs/su2_suzuki/monopoless/CON_OFF_MAG_001.LAT";
  // string conf_path = "../../confs/su3/conf.0501";
  string conf_path = "../../confs/qc2dstag/mu0.05/s0/CONF0201";
  // string conf_path =
  //     "../../confs/su2_suzuki/48^4/beta2.7/monopole/CON_MON_MAG_003.LAT";
  // string conf_path =
  //     "../../confs/su2_suzuki/48^4/beta2.7/monopoless/CON_OFF_MAG_003.LAT";
  // string conf_path = "../../confs/qc2dstag/40^4/mu0.00/CONF0201";
  // conf.read_double(conf_path, 4);
  conf.read_double_qc2dstag(conf_path);
  // conf.read_ildg(conf_path);
  // conf.read_float(conf_path, 4);
  // conf.read_double_convert_abelian(conf_path, 8);
  // conf.read_double_qc2dstag_convert_abelian(conf_path);

  // for (int mu = 0; mu < 4; mu++) {
  // link.move(mu, 1);
  //   cout << conf.array[link.place + mu] << endl;
  // }

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

  std::vector<std::vector<MATRIX_TYPE>> conf_separated =
      separate_wilson(conf.array);

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
