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

using namespace std;

int x_size;
int y_size;
int z_size;
int t_size;

int main(int argc, char *argv[]) {
  x_size = 32;
  y_size = 32;
  z_size = 32;
  t_size = 32;
  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;
  link1 link(x_size, y_size, z_size, t_size);
  data<su2> conf;
  data<abelian> conf_abelian;
  data<su2> conf_offd;
  data<su2> conf_qc2dstag;
  string path_su2 = "../../confs/su2/time_32/mu0.00/conf_0001.fl";
  string path_abelian = "../../confs/su2/abelian/CON_MON_MAG_031.LAT";
  string path_offd = "../../confs/offd/mu0.00/CON_OFF_MAG_033.LAT";
  string path_qc2dstag = "../../confs/qc2dstag/mu0.05/CONF0001";
  string path_ml5 = "/home/ilya/code/test/qc2dstag_reading/MAG_conf/beta5.ml5";
  conf.read_float(path_su2);
  conf_abelian.read_float_fortran(path_abelian);
  conf_offd.read_double_fortran(path_offd);

  x_size = 8;
  y_size = 8;
  z_size = 8;
  t_size = 2;

  data<su2> conf_ml5;
  vector<float> array_ml5 = read_full_ml5(path_ml5, 1000);
  for (int i = 0; i < 1; i++) {
    conf_ml5.read_float_ml5(array_ml5, i);
    cout << "ml5 plaket " << plaket(conf_ml5.array) / 2 << endl;
    cout << "ml5 plaket_time " << plaket_time(conf_ml5.array) / 2 << endl;
    cout << "ml5 plaket_space " << plaket_space(conf_ml5.array) / 2 << endl;
    cout << "ml5 polyakov " << polyakov(conf_ml5.array) / 2 << endl;
    cout << "MAG_functional ml5 " << MAG_functional_su2(conf_ml5.array) << endl;

    int T_min = 1, T_max = 2;
    int R_min = 1, R_max = 8;

    vector<FLOAT> vec_wilson;
    start_time = clock();

    vec_wilson = wilson(conf_ml5.array, R_min, R_max, T_min, T_max);

    end_time = clock();
    search_time = end_time - start_time;
    cout << "wilson time: " << search_time * 1. / CLOCKS_PER_SEC << endl;
    cout << "wilson_loops:" << endl;
    for (int T = T_min; T <= T_max; T++) {
      for (int R = R_min; R <= R_max; R++) {
        cout << "T = " << T << " R = " << R << " "
             << vec_wilson[(R - R_min) + (T - T_min) * (R_max - R_min + 1)]
             << endl;
      }
    }
  }

  x_size = 40;
  y_size = 40;
  z_size = 40;
  t_size = 40;

  conf_qc2dstag.read_double_qc2dstag(path_qc2dstag);

  vector<vector<int>> directions;
  directions = generate_directions(4);

  start_time = clock();

  vector<wilson_result> wilson_offaxis_result =
      wilson_offaxis(conf_qc2dstag.array, directions, 4.9, 6, 10, 10);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "offaxis wilson loop time: " << search_time * 1. / CLOCKS_PER_SEC
       << endl;

  for (int i = 0; i < wilson_offaxis_result.size(); i++) {
    cout << "wilson loop: " << wilson_offaxis_result[i].wilson_loop
         << " time size: " << wilson_offaxis_result[i].time_size
         << " space size: " << wilson_offaxis_result[i].space_size << endl;
  }

  wilson_offaxis_reduce(wilson_offaxis_result);

  cout << endl << "after reuction" << endl;

  for (int i = 0; i < wilson_offaxis_result.size(); i++) {
    cout << "wilson loop: " << wilson_offaxis_result[i].wilson_loop
         << " time size: " << wilson_offaxis_result[i].time_size
         << " space size: " << wilson_offaxis_result[i].space_size << endl;
  }

  start_time = clock();
  cout << "MAG_functional_su2 " << MAG_functional_su2(conf_qc2dstag.array)
       << endl;
  end_time = clock();
  search_time = end_time - start_time;
  cout << " time: " << search_time * 1. / CLOCKS_PER_SEC << endl;

  start_time = clock();
  cout << "qc2dstag plaket " << plaket(conf_qc2dstag.array) / 2 << endl;
  cout << "qc2dstag plaket_time " << plaket_time(conf_qc2dstag.array) / 2
       << endl;
  cout << "qc2dstag plaket_space " << plaket_space(conf_qc2dstag.array) / 2
       << endl;
  cout << "qc2dstag polyakov " << polyakov(conf_qc2dstag.array) / 2 << endl;
  end_time = clock();
  search_time = end_time - start_time;
  cout << "plaket and staff time: " << search_time * 1. / CLOCKS_PER_SEC
       << endl;

  double start;
  double end;

  int T_min = 10, T_max = 10;
  int R_min = 5, R_max = 6;

  vector<FLOAT> vec_wilson;
  start_time = clock();

  vec_wilson = wilson(conf_qc2dstag.array, R_min, R_max, T_min, T_max);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "wilson time: " << search_time * 1. / CLOCKS_PER_SEC << endl;
  cout << "wilson_loops:" << endl;
  for (int T = T_min; T <= T_max; T++) {
    for (int R = R_min; R <= R_max; R++) {
      cout << "T = " << T << " R = " << R << " "
           << vec_wilson[(R - R_min) + (T - T_min) * (R_max - R_min + 1)]
           << endl;
    }
  }

  int R = 20;
  int T = 8;
  int d_min = -5;
  int d_max = 0;

  start_time = clock();
  vector<FLOAT> wilson_vec =
      calculate_wilson_loop_tr(conf_qc2dstag.array, R, T);
  vector<FLOAT> plaket_time_vec = calculate_plaket_time_tr(conf_qc2dstag.array);
  vector<FLOAT> plaket_space_vec =
      calculate_plaket_space_tr(conf_qc2dstag.array);
  end_time = clock();
  search_time = end_time - start_time;
  cout << "prepare time: " << search_time * 1. / CLOCKS_PER_SEC << endl;

  start_time = clock();
  result electric = wilson_plaket_correlator_electric(
      wilson_vec, plaket_time_vec, R, T, 0, d_min, d_max);
  end_time = clock();
  search_time = end_time - start_time;
  cout << "electric time: " << search_time * 1. / CLOCKS_PER_SEC << endl;

  start_time = clock();
  result magnetic = wilson_plaket_correlator_magnetic(
      wilson_vec, plaket_space_vec, R, T, 0, d_min, d_max);
  end_time = clock();
  search_time = end_time - start_time;

  cout << "magnetic time: " << search_time * 1. / CLOCKS_PER_SEC << endl;
  for (int i = 0; i < electric.array.size(); i++) {
    cout << "electric R = " << R << " T= " << T << " " << electric.array[i]
         << endl;
    cout << "magnetic R = " << R << " T= " << T << " " << magnetic.array[i]
         << endl;
  }

  int x_trans_min = -12;
  int x_trans_max = 12;

  start_time = clock();
  result electric_x = wilson_plaket_correlator_electric_x(
      wilson_vec, plaket_time_vec, R, T, x_trans_min, x_trans_max, R / 2);
  end_time = clock();
  search_time = end_time - start_time;
  cout << "electric_x time: " << search_time * 1. / CLOCKS_PER_SEC << endl;

  start_time = clock();
  result magnetic_x = wilson_plaket_correlator_magnetic_x(
      wilson_vec, plaket_space_vec, R, T, x_trans_min, x_trans_max, R / 2);
  end_time = clock();
  search_time = end_time - start_time;

  cout << "magnetic_x time: " << search_time * 1. / CLOCKS_PER_SEC << endl;
  for (int i = 0; i < electric_x.array.size(); i++) {
    cout << "electric_x R = " << R << " T= " << T << " " << electric_x.array[i]
         << endl;
    cout << "magnetic_x R = " << R << " T= " << T << " " << magnetic_x.array[i]
         << endl;
  }

  x_size = 32;
  y_size = 32;
  z_size = 32;
  t_size = 32;
  link1 link_abelian(x_size, y_size, z_size, t_size);
  cout.precision(10);

  start_time = clock();

  cout << "abelian" << endl;
  cout << "test plaket " << plaket(conf_abelian.array) / 2 << endl;
  cout << "test plaket_time " << plaket_time(conf_abelian.array) / 2 << endl;
  cout << "test plaket_space " << plaket_space(conf_abelian.array) / 2 << endl;
  cout << "test polyakov_loop " << polyakov(conf_abelian.array) / 2 << endl;
  vector<FLOAT> vec_wilson_abelian;
  vec_wilson_abelian = wilson(conf_abelian.array, 10, 10, 6, 6);
  cout << "test wilson_loop_R = 10_T = 6 " << vec_wilson_abelian[0] << endl;

  end_time = clock();
  search_time = end_time - start_time;
  cout << "abelian working time: " << search_time * 1. / CLOCKS_PER_SEC << endl;

  start_time = clock();

  cout << "offd" << endl;
  cout << "test plaket " << plaket(conf_offd.array) / 2 << endl;
  cout << "test plaket_time " << plaket_time(conf_offd.array) / 2 << endl;
  cout << "test plaket_space " << plaket_space(conf_offd.array) / 2 << endl;
  cout << "test polyakov_loop " << polyakov(conf_offd.array) / 2 << endl;

  vector<FLOAT> vec_wilson_offd;
  vec_wilson_offd = wilson(conf_offd.array, 10, 10, 6, 6);
  cout << "test wilson_loop_R=10_T=6 " << vec_wilson_offd[0] << endl;

  end_time = clock();
  search_time = end_time - start_time;
  cout << "offd working time: " << search_time * 1. / CLOCKS_PER_SEC << endl;

  start_time = clock();

  cout << "su2" << endl;
  cout << "test plaket " << plaket(conf.array) / 2 << " right: 0.6769540066 "
       << endl;
  cout << "test plaket_time " << plaket_time(conf.array) / 2
       << " right: 0.6770628794" << endl;
  cout << "test plaket_space " << plaket_space(conf.array) / 2
       << " right: 0.6768451339" << endl;
  cout << "test polyakov_loop " << polyakov(conf.array) / 2
       << " right: -0.004586235468" << endl;

  vector<FLOAT> vec_wilson_su2;
  vec_wilson_su2 = wilson(conf.array, 10, 10, 6, 6);
  cout << "test wilson_loop_R=10_T=6 " << vec_wilson_su2[0]
       << " right: 0.001178588784" << endl;

  end_time = clock();
  search_time = end_time - start_time;
  cout << "working time: " << search_time * 1. / CLOCKS_PER_SEC << endl;
}
