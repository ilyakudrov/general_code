#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
// #include "../../../lib/cpu/include/decomposition.h"
// #include "../../../lib/cpu/include/flux_tube.h"
#include "../../../lib/cpu/include/link.h"
// #include "../../../lib/cpu/include/mag.h"
#include "../../../lib/cpu/include/matrix.h"
// #include "../../../lib/cpu/include/monopoles.h"
#include "../../../lib/cpu/include/smearing.h"

// #include <algorithm>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <omp.h>
#include <vector>

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

  x_size = 40;
  y_size = 40;
  z_size = 40;
  t_size = 40;
  size1 = x_size * y_size;
  size2 = x_size * y_size * z_size;

  std::cout.precision(17);

  Data::data<MATRIX_TYPE> conf1;
  // data<MATRIX_TYPE> conf2;

  //   string conf_path1 = "../../confs/su3/gluodynamics/24^4/beta6.0/CONF0175";
  string conf_path1 = "../../confs/MAG/su3/gluodynamics/40^4/beta6.4/steps_0/"
                      "copies=20/s1/conf_gaugefixed_0002_1";
  // string conf_path1 =
  //     "/home/ilya/soft/lattice/general_code/tests/confs/monopole/su2/"
  //     "qc2dstag/40^4/mu0.00/conf_monopole_0201";
  string conf_path2 = "../../confs/su3/QCD/140MeV/nt20/conf.0501";
  string conf_format1 = "double";
  string conf_format2 = "ildg";
  int bytes_skip = 0;
  bool convert = 0;

  get_data(conf1, conf_path1, conf_format1, bytes_skip, convert);
  // vector<float> conf_full = read_full_ml5(conf_path1, 1);
  // conf1.read_float_ml5(conf_full, 0);
  // get_data(conf2, conf_path2, conf_format2, bytes_skip, convert);

  // plakets and polyakov loop
  start_time = omp_get_wtime();
  std::cout << "plaket " << plaket(conf1.array) << std::endl;
  std::cout << "plaket_time " << plaket_time(conf1.array) << std::endl;
  std::cout << "plaket_space " << plaket_space(conf1.array) << std::endl;
  std::cout << "polyakov " << polyakov_loop(conf1.array) << std::endl;
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "plaket and staff time: " << search_time << std::endl;
  // start_time = omp_get_wtime();
  // std::cout << "plaket " << plaket(conf2.array) << std::endl;
  // std::cout << "plaket_time " << plaket_time(conf2.array) << std::endl;
  // std::cout << "plaket_space " << plaket_space(conf2.array) << std::endl;
  // std::cout << "polyakov " << polyakov_loop(conf2.array) << std::endl;
  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // std::cout << "plaket and staff time: " << search_time << std::endl;

  // cout << "MAG functional " << MAG_functional_su2(conf1.array) << endl;

  std::vector<std::vector<MATRIX_TYPE>> conf_separated =
      separate_wilson(conf1.array);
  start_time = omp_get_wtime();
  std::cout << "polyakov_parallel " << polyakov_loop_parallel(conf_separated)
            << endl;
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "polyakov_parallel time: " << search_time << std::endl;
  // std::vector<std::vector<MATRIX_TYPE>> conf_separated =
  // separate_3(conf1.array);

  // on-axis wilson loops
  int T_min = 1, T_max = 20;
  int R_min = 1, R_max = 20;

  // std::vector<double> wilson_test =
  //     wilson(conf2.array, R_min, R_max, T_min, T_max);

  // std::cout << std::endl;
  // for (int i = 0; i < wilson_test.size(); i++) {
  //   std::cout << wilson_test[i] << std::endl;
  // }

  // off-axis wilson loops

  // start_time = omp_get_wtime();

  // map<tuple<int, double>, double> wilson_loops =
  //     wilson_offaxis_result(conf1.array, 0.9, 4, 1, 4);

  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // std::cout << "offaxis wilson loop time: " << search_time << std::endl;

  // for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
  //   std::cout << get<0>(it->first) << "," << get<1>(it->first) << ","
  //             << it->second << endl;
  // }
}
