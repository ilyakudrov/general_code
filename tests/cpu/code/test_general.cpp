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
#include <omp.h>
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
  double start_time;
  double end_time;
  double search_time;

  x_size = 24;
  y_size = 24;
  z_size = 24;
  t_size = 24;

  std::cout.precision(17);

  data<MATRIX_TYPE> conf;

  string conf_path = "../../confs/MA_gauge/su3/gluodynamics/24^4/beta6.0/"
                     "steps_500/copies=3/conf_gaugefixed_0001";
  // string conf_path =
  //     "../../confs/MA_gauge/su2/su2_suzuki/48^4/beta2.6/conf_0001";
  // string conf_path = "../../confs/Landau_U1xU1/gluodynamics/24^4/beta6.0/"
  //                    "steps_500/copies=3/conf_Landau_gaugefixed_0001";
  // string conf_path =
  // "/home/ilya/soft/lattice/general_code/apps/smearing/test/"
  //                    "result/smeared_0001";
  string conf_format = "double_qc2dstag";
  int bytes_skip = 0;
  bool convert = 0;

  get_data(conf, conf_path, conf_format, bytes_skip, convert);

  // plakets and polyakov loop
  start_time = omp_get_wtime();
  std::cout << "plaket " << plaket(conf.array) << std::endl;
  std::cout << "plaket_time " << plaket_time(conf.array) << std::endl;
  std::cout << "plaket_space " << plaket_space(conf.array) << std::endl;
  std::cout << "polyakov " << polyakov_loop(conf.array) << std::endl;
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "plaket and staff time: " << search_time << std::endl;

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
  start_time = omp_get_wtime();

  map<tuple<int, int>, double> wilson_loops =
      wilson_adjoint_parallel(conf_separated, R_min, R_max, T_min, T_max);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "on-axis wilson time: " << search_time << std::endl;
  std::cout << "wilson_loops adjoint:" << std::endl;
  for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
    cout << get<0>(it->first) << "," << get<1>(it->first) << "," << it->second
         << endl;
  }

  // off-axis wilson loops
  std::vector<std::vector<int>> directions;
  directions = generate_directions(4);

  start_time = omp_get_wtime();

  std::vector<wilson_result> wilson_offaxis_result =
      wilson_offaxis(conf.array, directions, 0.9, 4, 1, 4);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "offaxis wilson loop time: " << search_time << std::endl;

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
