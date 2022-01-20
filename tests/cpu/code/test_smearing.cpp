#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/result.h"
#include "../../../lib/cpu/include/smearing.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <stdio.h>

int x_size = 32;
int y_size = 32;
int z_size = 32;
int t_size = 32;

int main(int argc, char *argv[]) {
  double alpha1 = 1;
  double alpha2 = 1;
  double alpha3 = 0.5;
  double alpha_APE = 0.5;
  double stout_rho = 0.15;

  std::vector<std::vector<su2>> smearing_first(9, std::vector<su2>(1));
  std::vector<std::vector<su2>> smearing_second(6, std::vector<su2>(1));

  link1 link(x_size, y_size, z_size, t_size);
  link1 link_double(x_size, y_size, z_size, t_size);
  data<su2> conf;
  data<su2> smeared;
  data<abelian> conf_abelian;
  data<abelian> conf_abelian_smeared;
  std::string path_su2 = "../../confs/su2_dynam/32^4/CON_32^3x32_001.LAT";
  std::string path_abelian = "../../confs/su2/time_32/mu0.00/conf_0001.fl";
  conf.read_float_fortran(path_su2);
  conf_abelian.read_float_fortran(path_abelian);
  double aver[2];

  std::cout.precision(10);
  su2 A;
  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  std::map<std::tuple<int, int>, std::vector<su2>> APE_2d_test =
      smearing_APE_2d(conf.array, alpha_APE);

  smearing_APE_2d_continue(APE_2d_test, alpha_APE);

  start_time = clock();

  double start;
  double end;

  smearing_first = smearing_first_full(conf.array, alpha3);
  smearing_second = smearing_second_full(conf.array, smearing_first, alpha2);
  smeared.array = smearing_HYP(conf.array, smearing_second, alpha1);

  for (int i = 0; i < 4; i++) {
    std::cout << "smearing_HYP test " << smeared.array[i] << std::endl;
  }
  std::cout << "right:" << std::endl;
  std::cout << "-0.110943 0.629543 -0.583767 0.500582" << std::endl;
  std::cout << "0.318657 -0.218438 -0.415944 0.823245" << std::endl;
  std::cout << "0.174894 0.667862 0.674797 0.260808" << std::endl;
  std::cout << "-0.704066 -0.13241 0.453476 -0.530206" << std::endl;

  smeared.array = smearing_HYP_refresh(conf, alpha1, alpha2, alpha3);

  for (int i = 0; i < 4; i++) {
    A = smeared.array[i];
    std::cout << "smearing_HYP_refresh test " << A.a0 << " " << A.a1 << " "
              << A.a2 << " " << A.a3 << std::endl;
  }

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

  for (int i = 0; i < 4; i++) {
    A = smeared.array[i];
    std::cout << "smearing_APE test " << A.a0 << " " << A.a1 << " " << A.a2
              << " " << A.a3 << std::endl;
  }
  std::cout << "right:" << std::endl;
  std::cout << "-0.0833883 0.641689 -0.604086 0.465147" << std::endl;
  std::cout << "0.336007 -0.177589 -0.426239 0.820903" << std::endl;
  std::cout << "0.125318 0.626699 0.730398 0.240963" << std::endl;
  std::cout << "-0.632818 0.13761 0.745239 0.158818" << std::endl;

  smeared.array = smearing_stout(conf, stout_rho);

  for (int i = 0; i < 4; i++) {
    A = smeared.array[i];
    std::cout << "smearing_stout test " << A.a0 << " " << A.a1 << " " << A.a2
              << " " << A.a3 << std::endl;
  }
  std::cout << "right:" << std::endl;
  std::cout << "-0.0185375 0.684286 -0.660024 0.30948" << std::endl;
  std::cout << "0.363087 -0.0991993 -0.575234 0.726246" << std::endl;
  std::cout << "0.172742 0.633665 0.73173 0.182208" << std::endl;
  std::cout << "-0.816237 -0.0404862 0.571983 -0.0703786" << std::endl;

  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "working time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;
}
