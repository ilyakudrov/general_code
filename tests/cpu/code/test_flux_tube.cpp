#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/flux_tube.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <omp.h>
#include <stdio.h>
#include <unordered_map>
#include <vector>

using namespace std;

int x_size;
int y_size;
int z_size;
int t_size;

#define MATRIX_TYPE su2

int main(int argc, char *argv[]) {
  x_size = 40;
  y_size = 40;
  z_size = 40;
  t_size = 40;

  double start_time;
  double end_time;
  double search_time;

  data<MATRIX_TYPE> conf_qc2dstag;
  data<MATRIX_TYPE> conf_smeared;
  string path_qc2dstag = "../../confs/su2/qc2dstag/40^4/mu0.00/CONF0201";
  string path_smeared = "../../confs/smeared/qc2dstag/40^4/mu0.00/"
                        "HYP0_alpha=1_1_0.5_APE_alpha=0.5/smeared_0201";

  string conf_format = "double_qc2dstag";
  int bytes_skip = 0;
  bool convert = 0;
  get_data(conf_qc2dstag, path_qc2dstag, conf_format, bytes_skip, convert);

  string conf_format_smeared = "double_qc2dstag";
  int bytes_skip_smeared = 0;
  bool convert_smeared = 0;
  get_data(conf_smeared, path_smeared, conf_format_smeared, bytes_skip_smeared,
           convert_smeared);

  cout << "schwinger wilson flux tube test" << endl;

  int R = 4;
  int T = 4;

  start_time = omp_get_wtime();
  vector<MATRIX_TYPE> plaket_schwinger_electric =
      calculate_plaket_time(conf_qc2dstag.array);
  end_time = clock();
  search_time = end_time - start_time;
  cout << "calculate_plaket_schwinger_time time: "
       << search_time * 1. / CLOCKS_PER_SEC << endl;

  vector<vector<MATRIX_TYPE>> schwinger_lines_short(R + 5,
                                                    vector<MATRIX_TYPE>());
  for (int d = 0; d < R + 5; d++) {
    schwinger_lines_short[d] =
        calculate_schwinger_lines_short(conf_smeared.array, d + 1);
  }

  start_time = omp_get_wtime();
  map<int, double> schwinger_electric =
      wilson_plaket_schwinger_electric_longitudinal(
          conf_smeared.array, plaket_schwinger_electric, schwinger_lines_short,
          -5, R + 5, T, R);
  end_time = clock();
  search_time = end_time - start_time;
  cout << "wilson_plaket_schwinger_electric_longitudinal time: "
       << search_time * 1. / CLOCKS_PER_SEC << endl;

  for (auto it = schwinger_electric.begin(); it != schwinger_electric.end();
       ++it) {
    cout << "d = " << it->first << " schwinger_electric = " << it->second
         << endl;
  }

  cout << "wilson plaket flux tube test" << endl;

  R = 8;
  T = 10;
  int d_min = -5;
  int d_max = R + 5;

  start_time = omp_get_wtime();
  vector<double> wilson_vec =
      calculate_wilson_loop_tr(conf_smeared.array, R, T);
  vector<double> plaket_time_vec =
      calculate_plaket_time_tr(conf_qc2dstag.array);
  vector<double> plaket_space_vec =
      calculate_plaket_space_tr(conf_qc2dstag.array);
  double wilson_loop = 0;
  double plaket = 0;
  for (int i = 0; i < wilson_vec.size(); i++) {
    wilson_loop += wilson_vec[i];
  }
  for (int i = 0; i < plaket_time_vec.size(); i++) {
    plaket += plaket_time_vec[i];
  }

  cout << "wilson_loop " << wilson_loop / wilson_vec.size() << endl;
  cout << "plaket " << plaket / plaket_time_vec.size() << endl;
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "prepare time: " << search_time << endl;

  start_time = omp_get_wtime();
  map<int, double> electric = wilson_plaket_correlator_electric(
      wilson_vec, plaket_time_vec, R, T, 0, d_min, d_max);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "longitudinal  electric time: " << search_time << endl;

  start_time = omp_get_wtime();
  map<int, double> magnetic = wilson_plaket_correlator_magnetic(
      wilson_vec, plaket_space_vec, R, T, 0, d_min, d_max);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "longitudinal magnetic time: " << search_time << endl;

  for (auto it = electric.begin(); it != electric.end(); ++it) {
    cout << "longitudinal electric R = " << R << " T = " << T << " d "
         << it->first << " " << it->second << endl;
  }
  cout << endl;
  for (auto it = magnetic.begin(); it != magnetic.end(); ++it) {
    cout << "longitudinal magnetic R = " << R << " T = " << T << " d "
         << it->first << " " << it->second << endl;
  }

  start_time = omp_get_wtime();
  map<int, double> electric_trans = wilson_plaket_correlator_electric_x(
      wilson_vec, plaket_time_vec, R, T, 5, R / 2);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "transversal electric time: " << search_time << endl;

  start_time = omp_get_wtime();
  map<int, double> magnetic_trans = wilson_plaket_correlator_magnetic_x(
      wilson_vec, plaket_space_vec, R, T, 5, R / 2);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "transversal magnetic time: " << search_time << endl;

  for (auto it = electric_trans.begin(); it != electric_trans.end(); ++it) {
    cout << "transversal electric R = " << R << " T = " << T << " d "
         << it->first << " " << it->second << endl;
  }
  cout << endl;
  for (auto it = magnetic_trans.begin(); it != magnetic_trans.end(); ++it) {
    cout << "transversal magnetic R = " << R << " T = " << T << " d "
         << it->first << " " << it->second << endl;
  }

  vector<vector<MATRIX_TYPE>> separated_unchanged =
      separate_wilson(conf_qc2dstag.array);
  vector<vector<MATRIX_TYPE>> separated_smeared =
      separate_wilson(conf_smeared.array);

  int T_min = 10, T_max = 10;
  int R_min = 8, R_max = 8;

  map<tuple<int, int, int>, double> flux_tube_new;

  vector<double> plaket_time_tr = plaket_aver_tr_time(separated_unchanged);

  vector<double> plaket_space_tr = plaket_aver_tr_space(separated_unchanged);

  start_time = omp_get_wtime();

  flux_tube_new =
      wilson_plaket_correlator(plaket_time_tr, separated_smeared, T_min, T_max,
                               R_min, R_max, 5, 0, "longitudinal");

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "flux tube longitudinal electric new time: " << search_time << endl;

  for (auto it = flux_tube_new.begin(); it != flux_tube_new.end(); it++) {
    cout << "T = " << get<0>(it->first) << ", R = " << get<1>(it->first)
         << ", d = " << get<2>(it->first)
         << ", longitudinal electric flux_tube = " << it->second << endl;
  }

  start_time = omp_get_wtime();

  flux_tube_new =
      wilson_plaket_correlator(plaket_space_tr, separated_smeared, T_min, T_max,
                               R_min, R_max, 5, 0, "longitudinal");

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "flux tube longitudinal magnetic new time: " << search_time << endl;

  for (auto it = flux_tube_new.begin(); it != flux_tube_new.end(); it++) {
    cout << "T = " << get<0>(it->first) << ", R = " << get<1>(it->first)
         << ", d = " << get<2>(it->first)
         << ", longitudinal magnetic flux_tube = " << it->second << endl;
  }

  start_time = omp_get_wtime();

  flux_tube_new =
      wilson_plaket_correlator(plaket_time_tr, separated_smeared, T_min, T_max,
                               R_min, R_max, 5, R / 2, "transversal");

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "flux tube transversal electric new time: " << search_time << endl;

  for (auto it = flux_tube_new.begin(); it != flux_tube_new.end(); it++) {
    cout << "T = " << get<0>(it->first) << ", R = " << get<1>(it->first)
         << ", d = " << get<2>(it->first)
         << ", transversal electric flux_tube = " << it->second << endl;
  }

  start_time = omp_get_wtime();

  flux_tube_new =
      wilson_plaket_correlator(plaket_space_tr, separated_smeared, T_min, T_max,
                               R_min, R_max, 5, R / 2, "transversal");

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "flux tube transversal magnetic new time: " << search_time << endl;

  for (auto it = flux_tube_new.begin(); it != flux_tube_new.end(); it++) {
    cout << "T = " << get<0>(it->first) << ", R = " << get<1>(it->first)
         << ", d = " << get<2>(it->first)
         << ", transversal magnetic flux_tube = " << it->second << endl;
  }
}