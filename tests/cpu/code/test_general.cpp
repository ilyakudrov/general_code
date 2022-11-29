#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/decomposition.h"
#include "../../../lib/cpu/include/flux_tube.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/mag.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/monopoles.h"
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

  x_size = 64;
  y_size = 64;
  z_size = 64;
  t_size = 4;

  std::cout.precision(17);

  data<MATRIX_TYPE> conf;
  // string conf_path =
  //     "../../confs/SU3_conf/gluodynamics/16^4/su3_mag_u1.01001.lat";
  // string conf_path =
  //     "../../confs/decomposed/monopoless/gluodynamics/16^4/MLS_conf.01001.lat";
  // string conf_path = "../../confs/decomposed/monopole/gluodynamics/16^4/"
  //                    "CON_MON_MAG_01001.LAT";
  // string conf_path =
  //     "/home/ilya/soft/lattice/general_code/apps/monopole_decomposition_su3/"
  //     "test/result/conf_monopole_16_1001_non-unitary";
  // string conf_path =
  //     "/home/ilya/soft/lattice/general_code/apps/monopole_decomposition_su3/"
  //     "test/result/conf_monopole_16_1001";
  // string conf_path = "../../confs/decomposed/monopole/gluodynamics/16^4/old/"
  //                    "CON_MON_MAG_01001.LAT";
  string conf_path = "../../confs/MA_gauge/su3/QCD/140MeV/nt4/steps_2000/"
                     "copies=1/tol=1e-13/conf.SP_gaugefixed_0501";

  // conf.read_double(conf_path, 0);
  // conf.read_double_vitaly(conf_path, 4);
  // conf.read_double_qc2dstag(conf_path);
  conf.read_ildg(conf_path);
  // conf.read_float(conf_path, 4);
  // conf.read_double_convert_abelian(conf_path, 8);
  // conf.read_double_qc2dstag_convert_abelian(conf_path);

  cout << "mag functional " << mag_functional_su3(conf.array) << endl;

  // double sum;
  // double sum_aver = 0;
  // double sum_max = 0;
  // for (int i = 0; i < conf.array.size(); i++) {
  //   sum = 0;
  //   for (int j = 0; j < 3; j++) {
  //     sum += atan2(conf.array[i].matrix[j].imag,
  //     conf.array[i].matrix[j].real);
  //   }
  //   sum_aver += sum;
  //   if (sum > sum_max) {
  //     sum_max = sum;
  //   }
  //   // cout << "sum " << sum << endl;
  // }
  // sum_aver = sum_aver / conf.array.size();
  // cout << "sum_aver = " << sum_aver << " sum_max =  " << sum_max << endl;

  // vector<vector<double>> angles = convert_to_angles(conf.array);

  // string conf_path1 =
  //     "../../confs/SU3_conf/gluodynamics/16^4/su3_mag_u1.01001.lat";

  // data<su3> conf1;
  // conf1.read_double_vitaly(conf_path1, 4);
  // std::cout << "plaket su3 " << plaket(conf1.array) << std::endl;

  // get_monopoless_optimized_su3(conf1.array, angles);
  // std::cout << "plaket monopoless1 " << plaket(conf1.array) << std::endl;

  // string conf_path2 =
  //     "../../confs/decomposed/monopoless/gluodynamics/16^4/MLS_conf.01001.lat";
  // data<su3> conf2;
  // conf2.read_double_vitaly(conf_path2, 4);
  // std::cout << "plaket monopoless2 " << plaket(conf2.array) << std::endl;

  // plakets and polyakov loop
  start_time = clock();
  std::cout << "plaket " << plaket(conf.array) << std::endl;
  std::cout << "plaket_time " << plaket_time(conf.array) << std::endl;
  std::cout << "plaket_space " << plaket_space(conf.array) << std::endl;
  std::cout << "polyakov " << polyakov_loop(conf.array) << std::endl;
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
  int T_min = 1, T_max = 4;
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
