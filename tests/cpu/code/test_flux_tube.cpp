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
  x_size = 40;
  y_size = 40;
  z_size = 40;
  t_size = 40;
  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;
  link1 link(x_size, y_size, z_size, t_size);
  data<su2> conf_qc2dstag;
  data<su2> conf_smeared;
  std::string path_qc2dstag = "../../confs/qc2dstag/mu0.05/s0/CONF0201";
  std::string path_smeared =
      "../../confs/qc2dstag/HYP_APE/mu0.05/s0/smeared_0201";

  conf_qc2dstag.read_double_qc2dstag(path_qc2dstag);
  conf_smeared.read_double(path_smeared);

  std::cout << "schwinger wilson flux tube test" << std::endl;

  int R = 10;
  int T = 10;
  start_time = clock();
  std::vector<std::vector<std::vector<su2>>> schwinger_lines_short;
  for (int d = 1; d <= R + 5; d++) {
    schwinger_lines_short.push_back(
        calculate_schwinger_lines_short(conf_qc2dstag.array, d));
  }
  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "schwinger_lines_short time: "
            << search_time * 1. / CLOCKS_PER_SEC << std::endl;

  start_time = clock();
  std::vector<su2> plaket_electric =
      calculate_plaket_schwinger_time(conf_qc2dstag.array);
  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "calculate_plaket_schwinger_time time: "
            << search_time * 1. / CLOCKS_PER_SEC << std::endl;

  start_time = clock();
  std::vector<FLOAT> schwinger_electric = wilson_plaket_schwinger_electric(
      conf_qc2dstag.array, schwinger_lines_short, plaket_electric, -5, R + 5, T,
      R);
  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "wilson_plaket_schwinger_electric time: "
            << search_time * 1. / CLOCKS_PER_SEC << std::endl;

  for (int i = 0; i < schwinger_electric.size(); i++) {
    std::cout << "d = " << i - 5
              << " schwinger_electric = " << schwinger_electric[i] << std::endl;
  }

  std::cout << "wilson plaket flux tube test" << std::endl;

  R = 14;
  T = 8;
  int d_min = -10;
  int d_max = R;

  start_time = clock();
  std::vector<FLOAT> wilson_vec =
      calculate_wilson_loop_tr(conf_smeared.array, R, T);
  std::vector<FLOAT> plaket_time_vec =
      calculate_plaket_time_tr(conf_qc2dstag.array);
  std::vector<FLOAT> plaket_space_vec =
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
    std::cout << "electric R = " << R << " T= " << T << " d " << i + d_min
              << " " << electric[i] << std::endl;
  }
  std::cout << std::endl;
  for (int i = 0; i < magnetic.size(); i++) {
    std::cout << "magnetic R = " << R << " T= " << T << " d " << i + d_min
              << " " << magnetic[i] << std::endl;
  }
}
