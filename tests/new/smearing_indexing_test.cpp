#include "../../lib/cpu/include/basic_observables.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/flux_tube.h"
#include "../../lib/cpu/include/link.h"
#include "../../lib/cpu/include/matrix.h"
#include "../../lib/cpu/include/smearing.h"
#include "include/indexing.h"
#include "include/smearing_indexing.h"
#include "include/wilson_loop_indexing.h"

#include <cmath>
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

#define MATRIX_TYPE su2

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
  // data<su3_full> conf;
  data<MATRIX_TYPE> conf;
  data<MATRIX_TYPE> conf1;
  // data<abelian> conf;
  string conf_path = "../confs/su2/su2_suzuki/48^4/beta2.8/CON_fxd_MAG_035.LAT";
  // string conf_path = "../confs/su2/CONF0201";
  // string conf_path = "../confs/su2_suzuki/24^4/beta2.4/CON_fxd_MAG_001.LAT";
  // string conf_path = "../confs/SU3_conf/nt14/conf.0501";
  conf.read_double(conf_path, 4);
  // conf.read_double_qc2dstag(conf_path);
  // conf.read_ildg(conf_path);

  double alpha_APE = 0.5;

  std::cout << "plaket before smearing " << plaket(conf.array) << std::endl;
  std::cout << "plaket time before smearing " << plaket_time(conf.array)
            << std::endl;
  std::cout << "plaket space before smearing " << plaket_space(conf.array)
            << std::endl;
  std::cout << "polyakov loop before smearing " << polyakov(conf.array)
            << std::endl;

  int T_min = 1, T_max = 5;
  int R_min = 1, R_max = 5;

  std::vector<double> vec_wilson =
      wilson(conf.array, R_min, R_max, T_min, T_max);

  std::cout << "wilson_loops before smearing:" << std::endl;
  for (int T = T_min; T <= T_max; T++) {
    for (int R = R_min; R <= R_max; R++) {
      std::cout << "T = " << T << " R = " << R << " "
                << vec_wilson[(R - R_min) + (T - T_min) * (R_max - R_min + 1)]
                << std::endl;
    }
  }

  start_time = omp_get_wtime();

  conf1.array = smearing1_APE(conf.array, alpha_APE);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "smearing old time: " << search_time << std::endl;

  // for (int i = 0; i < 16; i++) {
  //   std::cout << "conf.array after smearing, i = " << i << " " <<
  //   conf1.array[i]
  //             << std::endl;
  // }

  vec_wilson = wilson(conf1.array, R_min, R_max, T_min, T_max);

  std::cout << "wilson_loops after smearing:" << std::endl;
  for (int T = T_min; T <= T_max; T++) {
    for (int R = R_min; R <= R_max; R++) {
      std::cout << "T = " << T << " R = " << R << " "
                << vec_wilson[(R - R_min) + (T - T_min) * (R_max - R_min + 1)]
                << std::endl;
    }
  }

  std::cout << "plaket after smearing " << plaket(conf1.array) << std::endl;
  std::cout << "plaket time after smearing " << plaket_time(conf1.array)
            << std::endl;
  std::cout << "plaket space after smearing " << plaket_space(conf1.array)
            << std::endl;
  std::cout << "polyakov loop after smearing " << polyakov(conf1.array)
            << std::endl;

  std::vector<std::vector<MATRIX_TYPE>> conf_separated =
      separate_smearing_unchanged(conf.array);

  smearing_APE_test1(conf_separated, alpha_APE);

  cout << "plaket test " << plaket(conf1.array) << endl;

  int length_R = 1;
  int length_T = 1;

  std::vector<std::vector<MATRIX_TYPE>> wilson_lines(4);
  for (int i = 0; i < 3; i++) {
    wilson_lines[i] =
        wilson_lines_test3(conf_separated[i], length_R, link.multiplier[i] / 4,
                           link.multiplier[i + 1] / 4);
  }
  wilson_lines[3] =
      wilson_lines_test3(conf_separated[3], length_T, link.multiplier[3] / 4,
                         link.multiplier[3] / 4 * link.lattice_size[3]);

  std::cout << wilson_loop_test_time(wilson_lines, length_R, length_T, 4)
            << std::endl;

  double alpha1 = 0.75, alpha2 = 0.6, alpha3 = 0.3;

  std::vector<std::vector<su2>> smearing_first;
  std::vector<std::vector<su2>> smearing_second;

  start_time = omp_get_wtime();

  smearing_first = smearing_first_full(conf.array, alpha3);
  smearing_second = smearing_second_full(conf.array, smearing_first, alpha2);
  conf1.array = smearing_HYP(conf.array, smearing_second, alpha1);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "smearing HYP old time: " << search_time << std::endl;

  std::cout << "plaket after HYP smearing " << plaket(conf1.array) << std::endl;
  std::cout << "plaket time after HYP smearing " << plaket_time(conf1.array)
            << std::endl;
  std::cout << "plaket space after HYP smearing " << plaket_space(conf1.array)
            << std::endl;
  std::cout << "polyakov loop after HYP smearing " << polyakov(conf1.array)
            << std::endl;

  conf_separated = separate_smearing_unchanged(conf.array);

  start_time = omp_get_wtime();

  smearing_HYP_test1(conf_separated, alpha1, alpha2, alpha3);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "smearing_HYP_test1 time: " << search_time << std::endl;

  for (int i = 0; i < 3; i++) {
    wilson_lines[i] =
        wilson_lines_test3(conf_separated[i], length_R, link.multiplier[i] / 4,
                           link.multiplier[i + 1] / 4);
  }
  wilson_lines[3] =
      wilson_lines_test3(conf_separated[3], length_T, link.multiplier[3] / 4,
                         link.multiplier[3] / 4 * link.lattice_size[3]);

  std::cout << wilson_loop_test_time(wilson_lines, length_R, length_T, 4)
            << std::endl;
}
