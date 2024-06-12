#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/decomposition.h"
#include "../../../lib/cpu/include/flux_tube.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/mag.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/monopoles.h"
#include "../../../lib/cpu/include/smearing.h"

#include <algorithm>
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

#define MATRIX_TYPE su3_abelian

using namespace std;

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double search_time;

  x_size = 16;
  y_size = 16;
  z_size = 16;
  t_size = 16;

  std::cout.precision(17);

  data<MATRIX_TYPE> conf1;
  data<MATRIX_TYPE> conf2;

  string conf_path1 = "../../confs/MAG/su3/gluodynamics/16^4/beta6.0/steps_0/"
                      "copies=20/conf_gaugefixed_01001.lime_1";
  // string conf_path1 =
  //     "/home/ilya/soft/lattice/general_code/tests/confs/monopole/su2/"
  //     "qc2dstag/40^4/mu0.00/conf_monopole_0201";
  string conf_path2 =
      "/home/ilya/soft/lattice/general_code/tests/confs/"
      "monopoless/su2/qc2dstag/40^4/mu0.00/conf_monopoless_0201";
  string conf_format1 = "ildg";
  string conf_format2 = "double_qc2dstag";
  int bytes_skip = 0;
  bool convert = 1;

  get_data(conf1, conf_path1, conf_format1, bytes_skip, convert);
  // vector<float> conf_full = read_full_ml5(conf_path1, 1);
  // conf1.read_float_ml5(conf_full, 0);
  // get_data(conf2, conf_path2, conf_format2, bytes_skip, convert);

  // link1 link(x_size, y_size, z_size, t_size);
  // double a = 0.2;
  // su3 A, B;
  // // staple in dir 1
  // A = conf1.array[link.place + 1];
  // link.move(1, 1);
  // A = A * conf1.array[link.place];
  // link.move(1, -1);
  // link.move(0, 1);
  // B = A ^ conf1.array[link.place + 1];
  // link.move(0, -1);
  // link.move(1, -1);
  // A = conf1.array[link.place + 1].conj();
  // A = A * conf1.array[link.place];
  // link.move(1, 1);
  // B = B + A * conf1.array[link.place + 1];
  // link.move(0, -1);
  // link.move(1, -1);
  // // staple in dir 2
  // A = conf1.array[link.place + 2];
  // link.move(2, 1);
  // A = A * conf1.array[link.place];
  // link.move(2, -1);
  // link.move(0, 1);
  // B = B + A ^ conf1.array[link.place + 2];
  // link.move(0, -1);
  // link.move(2, -1);
  // A = conf1.array[link.place + 2].conj();
  // A = A * conf1.array[link.place];
  // link.move(2, 1);
  // B = B + A * conf1.array[link.place + 2];
  // link.move(0, -1);
  // link.move(2, -1);
  // B = (1 - a) * conf1.array[link.place] + (a / 4) * B;

  // cout << B << endl;
  // cout << B.determinant() << endl << endl;
  // cout << B.proj1() << endl;
  // cout << B.proj1().determinant() << endl << endl;
  // cout << B.proj() << endl;
  // cout << B.proj().determinant() << endl << endl;
  // cout << B.proj().proj() << endl;
  // cout << B.proj().proj().determinant() << endl << endl;
  // cout << B.proj().proj().proj() << endl;
  // cout << B.proj().proj().proj().determinant() << endl << endl;
  // cout << B.proj().proj().proj().proj() << endl;
  // cout << B.proj().proj().proj().proj().determinant() << endl << endl;
  // su3 C1 = B.proj1();
  // su3 C2 = B.proj().proj().proj().proj();
  // cout << (C1 ^ C1) << endl;
  // cout << (C2 ^ C2) << endl;
  // cout << (C1 - C2) << endl;

  // plakets and polyakov loop
  start_time = omp_get_wtime();
  std::cout << "plaket " << plaket(conf1.array) << std::endl;
  std::cout << "plaket_time " << plaket_time(conf1.array) << std::endl;
  std::cout << "plaket_space " << plaket_space(conf1.array) << std::endl;
  std::cout << "polyakov " << polyakov_loop(conf1.array) << std::endl;
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "plaket and staff time: " << search_time << std::endl;
  start_time = omp_get_wtime();
  std::cout << "plaket " << plaket(conf2.array) << std::endl;
  std::cout << "plaket_time " << plaket_time(conf2.array) << std::endl;
  std::cout << "plaket_space " << plaket_space(conf2.array) << std::endl;
  std::cout << "polyakov " << polyakov_loop(conf2.array) << std::endl;
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "plaket and staff time: " << search_time << std::endl;

  // cout << "MAG functional " << MAG_functional_su2(conf1.array) << endl;

  std::vector<std::vector<MATRIX_TYPE>> conf_separated =
      separate_wilson(conf2.array);
  // std::vector<std::vector<MATRIX_TYPE>> conf_separated =
  // separate_3(conf1.array);

  cout << "plaket parallel " << plaket_parallel(conf_separated) << endl;
  cout << "plaket time parallel " << plaket_time_parallel(conf_separated)
       << endl;
  cout << "plaket space parallel " << plaket_space_parallel(conf_separated)
       << endl;

  // on-axis wilson loops
  int T_min = 1, T_max = 16;
  int R_min = 1, R_max = 16;

  start_time = omp_get_wtime();

  map<tuple<int, int>, double> wilson_loops =
      wilson_parallel(conf_separated, R_min, R_max, T_min, T_max);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "on-axis wilson time: " << search_time << std::endl;
  std::cout << "wilson_loops adjoint:" << std::endl;
  for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
    cout << get<0>(it->first) << "," << get<1>(it->first) << "," << it->second
         << endl;
  }

  std::vector<double> wilson_test =
      wilson(conf2.array, R_min, R_max, T_min, T_max);

  std::cout << std::endl;
  for (int i = 0; i < wilson_test.size(); i++) {
    std::cout << wilson_test[i] << std::endl;
  }

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
