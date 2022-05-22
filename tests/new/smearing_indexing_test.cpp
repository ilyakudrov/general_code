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
  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  x_size = 24;
  y_size = 24;
  z_size = 24;
  t_size = 24;

  std::cout.precision(17);

  link1 link(x_size, y_size, z_size, t_size);
  // data<su2> conf;
  // data<su3_full> conf;
  data<MATRIX_TYPE> conf;
  data<MATRIX_TYPE> conf1;
  // data<abelian> conf;
  // string conf_path = "../confs/qc2dstag/40^4/mu0.00/CONF0201";
  string conf_path = "../confs/su2_suzuki/24^4/beta2.4/CON_fxd_MAG_001.LAT";
  // string conf_path = "../confs/su3/conf.0501";
  conf.read_double(conf_path, 8);
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

  int T_min = 1, T_max = 1;
  int R_min = 1, R_max = 1;

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

  start_time = clock();

  conf1.array = smearing1_APE(conf.array, alpha_APE);

  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "wilson_lines time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;

  // for (int i = 0; i < 16; i++) {
  //   std::cout << "conf.array after smearing, i = " << i << " " <<
  //   conf1.array[i]
  //             << std::endl;
  // }

  vec_wilson = wilson(conf1.array, R_min, R_max, T_min, T_max);

  std::cout << "wilson_loops after smearing:" << std::endl;
  for (int T = T_min; T <= T_max; T++) {
    for (int R = R_min; R <= R_max; R++) {
      std::cout
          << "T = " << T << " R = " << R << " "
          << 1 - vec_wilson[(R - R_min) + (T - T_min) * (R_max - R_min + 1)]
          << std::endl;
    }
  }

  std::cout << "plaket after smearing " << 1 - plaket(conf1.array) << std::endl;
  std::cout << "plaket time after smearing " << 1 - plaket_time(conf1.array)
            << std::endl;
  std::cout << "plaket space after smearing " << 1 - plaket_space(conf1.array)
            << std::endl;
  std::cout << "polyakov loop after smearing " << polyakov(conf1.array)
            << std::endl;

  std::vector<std::vector<MATRIX_TYPE>> conf_separated =
      separate_smearing_unchanged(conf.array);

  smearing_APE_test1(conf_separated, alpha_APE);

  // for (int mu = 0; mu < 4; mu++) {
  //   for (int i = 0; i < 4; i++) {
  //     std::cout << "conf_separated, mu = " << i << " " <<
  //     conf_separated[i][mu]
  //               << std::endl;
  //   }
  // }

  /*double test_difference = 0;
  double test_difference_max = 0;
  double test_difference_sum = 0;
  for (int i = 0; i < conf_separated[0].size(); i++) {
    // for (int i = 0; i < 100; i++) {
    for (int mu = 0; mu < 4; mu++) {
      test_difference =
          sqrt((conf_separated[mu][i] - conf1.array[i * 4 + mu]).module());
      // (conf_separated[mu][i] - conf1.array[i * 4 + mu]).a0;
      test_difference_sum += test_difference;
      if (test_difference_max < test_difference)
        test_difference_max = test_difference;
      // std::cout << test_difference << std::endl;
      // if (test_difference > 1e-12)
      //   std::cout << "difference if too big " << test_difference <<
      //   std::endl;
    }
  }
  std::cout << "difference max = " << test_difference_max << std::endl;
  std::cout << "difference sum = " << test_difference_sum << std::endl;
  std::cout << "difference aver = "
            << test_difference_sum / (conf_separated[0].size() * 4)
            << std::endl;*/

  // std::vector<std::vector<MATRIX_TYPE>> separated_unchanged =
  //     separate_wilson_unchanged(conf1.array);

  cout << "plaket test " << 1 - plaket(conf1.array) << endl;

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

  std::cout << 1 - wilson_loop_test_time(wilson_lines, length_R, length_T, 4)
            << std::endl;
}
