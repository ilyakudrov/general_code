#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/smearing.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <stdio.h>

int x_size = 24;
int y_size = 24;
int z_size = 24;
int t_size = 24;

int main(int argc, char *argv[]) {
  double alpha1 = 0.75;
  double alpha2 = 0.6;
  double alpha3 = 0.3;
  double alpha_APE = 0.5;
  double stout_rho = 0.15;

  std::vector<std::vector<su2>> smearing_first;
  std::vector<std::vector<su2>> smearing_second;

  link1 link(x_size, y_size, z_size, t_size);
  link1 link_double(x_size, y_size, z_size, t_size);
  data<su2> conf;
  data<su2> smeared;
  data<abelian> conf_abelian;
  data<abelian> conf_abelian_smeared;
  std::string path_su2 =
      "../../confs/su2_suzuki/24^4/beta2.4/CON_fxd_MAG_001.LAT";
  std::string path_abelian = "../../confs/su2/time_32/mu0.00/conf_0001.fl";
  conf.read_double(path_su2, 8);
  // conf_abelian.read_float(path_abelian, 0);
  double aver[2];

  std::cout.precision(17);
  su2 A;
  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  // std::map<std::tuple<int, int>, std::vector<su2>> APE_2d_test =
  //     smearing_APE_2d(conf.array, alpha_APE);

  // smearing_APE_2d_continue(APE_2d_test, alpha_APE);

  for (int i = 0; i < 4; i++) {
    std::cout << conf.array[i] << std::endl;
  }

  std::cout << "plaket: " << 1 - plaket(conf.array) << std::endl;
  std::cout << "plaket_time: " << 1 - plaket_time(conf.array) << std::endl;
  std::cout << "plaket_space: " << 1 - plaket_space(conf.array) << std::endl;
  std::cout << "polyakov_loop: " << polyakov(conf.array) << std::endl;

  start_time = clock();

  double start;
  double end;

  smearing_first = smearing_first_full(conf.array, alpha3);
  smearing_second = smearing_second_full(conf.array, smearing_first, alpha2);
  smeared.array = smearing_HYP(conf.array, smearing_second, alpha1);

  std::cout << "plaket: " << 1 - plaket(smeared.array) << std::endl;
  std::cout << "plaket_time: " << 1 - plaket_time(smeared.array) << std::endl;
  std::cout << "plaket_space: " << 1 - plaket_space(smeared.array) << std::endl;
  std::cout << "polyakov_loop: " << polyakov(smeared.array) << std::endl;

  smeared.array = smearing_HYP_refresh(conf, alpha1, alpha2, alpha3);

  int T_min = 1, T_max = 3;
  int R_min = 1, R_max = 3;

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
                << vec_wilson[(R - R_min) + (T - T_min) * (R_max - R_min + 1)] /
                       2
                << std::endl;
    }
  }

  smeared.array = smearing1_APE(conf.array, alpha_APE);

  for (int i = 0; i < 50; i++) {
    conf.array = smearing1_APE(conf.array, alpha_APE);
  }

  for (int i = 0; i < 4; i++) {
    std::cout << conf.array[i] << std::endl;
  }

  std::cout << "plaket: " << 1 - plaket(conf.array) << std::endl;
  std::cout << "plaket_time: " << 1 - plaket_time(conf.array) << std::endl;
  std::cout << "plaket_space: " << 1 - plaket_space(conf.array) << std::endl;
  std::cout << "polyakov_loop: " << polyakov(conf.array) << std::endl;

  start_time = clock();

  vec_wilson = wilson(smeared.array, R_min, R_max, T_min, T_max);

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

  start_time = clock();

  smeared.array = smearing_stout(conf, stout_rho);

  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "working time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;
}
