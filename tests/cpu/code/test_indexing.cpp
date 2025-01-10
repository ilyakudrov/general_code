#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/indexing.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/smearing.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <omp.h>

int x_size;
int y_size;
int z_size;
int t_size;
int size1;
int size2;

#define MATRIX_TYPE su3

using namespace std;

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double search_time;

  x_size = 64;
  y_size = 64;
  z_size = 64;
  t_size = 4;
  size1 = x_size * y_size;
  size2 = x_size * y_size * z_size;

  std::cout.precision(17);

  Data::data<MATRIX_TYPE> conf1;
  Data::data<MATRIX_TYPE> conf2;

  string conf_path1 = "../../confs/su3/QCD/140MeV/nt4/conf.0501";
  string conf_format1 = "ildg";
  // string conf_path1 =
  // "../../confs/MAG/su3/gluodynamics/40^4/beta6.4/steps_0/"
  //                     "copies=20/s1/conf_gaugefixed_0002_1";
  // string conf_format1 = "double";
  // string conf_path1 = "../../confs/su2/qc2dstag/40^4/mu0.00/CONF0201";
  // string conf_format1 = "double_qc2dstag";
  int bytes_skip = 0;
  bool convert = 0;

  int R_min = 1;
  int R_max = 4;
  int T_min = 1;
  int T_max = 4;

  map<tuple<int, int>, double> wilson_loops;

  get_data(conf1, conf_path1, conf_format1, bytes_skip, convert);
  conf2.array = conf1.array;

  start_time = omp_get_wtime();
  smearing_HYP_indexed(conf1.array, 1, 1, 0.5);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "smearing_HYP_parallel time: " << search_time << std::endl;

  start_time = omp_get_wtime();
  for (int i = 0; i < 10; i++) {
    smearing_APE_indexed(conf1.array, 0.5);
    smearing_APE_indexed(conf2.array, 0.5);
  }
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "smearing_APE_indexed time: " << search_time << std::endl;
  std::cout << conf1.array[4] << std::endl;
  std::cout << conf1.array[5] << std::endl;
  std::cout << conf1.array[6] << std::endl;
  std::cout << conf1.array[7] << std::endl;
  double trace_sum = 0;
  for (int i = 0; i < conf1.array.size(); i++) {
    trace_sum += conf1.array[i].tr();
  }
  std::cout << "trace " << trace_sum << std::endl;

  start_time = omp_get_wtime();
  wilson_loops =
      wilson_gevp_indexed(conf1.array, conf2.array, R_min, R_max, T_min, T_max);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "wilson_gevp_indexed time: " << search_time << std::endl;
  for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
    std::cout << get<0>(it->first) << ", " << get<1>(it->first) << ", "
              << it->second << endl;
  }

  get_data(conf1, conf_path1, conf_format1, bytes_skip, convert);
  conf2.array = conf1.array;
  std::vector<std::vector<MATRIX_TYPE>> conf_separated1 =
      separate_wilson(conf1.array);
  conf1.array.clear();
  conf1.array.shrink_to_fit();
  std::vector<std::vector<MATRIX_TYPE>> conf_separated2 =
      separate_wilson(conf2.array);
  conf2.array.clear();
  conf2.array.shrink_to_fit();

  start_time = omp_get_wtime();
  smearing_HYP_parallel(conf_separated1, 1, 1, 0.5);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "smearing_HYP_parallel time: " << search_time << std::endl;

  start_time = omp_get_wtime();
  for (int i = 0; i < 10; i++) {
    smearing_APE_parallel(conf_separated1, 0.5);
    smearing_APE_parallel(conf_separated2, 0.5);
  }
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "smearing_APE_parallel time: " << search_time << std::endl;
  std::cout << conf_separated1[0][1] << std::endl;
  std::cout << conf_separated1[1][1] << std::endl;
  std::cout << conf_separated1[2][1] << std::endl;
  std::cout << conf_separated1[3][1] << std::endl;
  trace_sum = 0;
  for (int mu = 0; mu < 4; mu++) {
    for (int i = 0; i < conf_separated1[mu].size(); i++) {
      trace_sum += conf_separated1[mu][i].tr();
    }
  }
  std::cout << "trace " << trace_sum << std::endl;

  start_time = omp_get_wtime();
  wilson_loops = wilson_gevp_parallel(conf_separated1, conf_separated2, R_min,
                                      R_max, T_min, T_max);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "wilson_gevp_parallel time: " << search_time << std::endl;
  for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
    std::cout << get<0>(it->first) << ", " << get<1>(it->first) << ", "
              << it->second << endl;
  }
}