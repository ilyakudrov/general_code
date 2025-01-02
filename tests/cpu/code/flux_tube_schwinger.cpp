// #include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
// #include "../../../lib/cpu/include/flux_tube.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
// #include <iostream>
// #include <map>
#include <omp.h>
// #include <vector>

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

  Data::data<MATRIX_TYPE> conf_qc2dstag;
  Data::data<MATRIX_TYPE> conf_smeared;
  string path_qc2dstag = "../../confs/su2/qc2dstag/40^4/mu0.00/CONF0201";
  string path_smeared = "../../confs/smeared/qc2dstag/40^4/mu0.00/"
                        "HYP0_alpha=1_1_0.5_APE_alpha=0.5/smeared_0201";

  string conf_format = "double_qc2dstag";
  int bytes_skip = 0;
  bool convert = 0;
  get_data(conf_qc2dstag, path_qc2dstag, conf_format, bytes_skip, convert);

  string conf_format_smeared = "double";
  int bytes_skip_smeared = 0;
  bool convert_smeared = 0;
  get_data(conf_smeared, path_smeared, conf_format_smeared, bytes_skip_smeared,
           convert_smeared);

  int R_min = 2;
  int R_max = 4;
  int T_min = 2;
  int T_max = 4;

  // cout << "schwinger wilson flux tube test" << endl;

  // start_time = omp_get_wtime();

  // std::map<std::tuple<int, int, int>, double> flux_tube_schwinger =
  //     flux_schwinger_electric_longitudinal(conf_qc2dstag.array,
  //                                          conf_smeared.array, T_min, T_max,
  //                                          R_min, R_max, 5);

  // end_time = clock();
  // search_time = end_time - start_time;
  // cout << "flux_schwinger_electric_longitudinal: "
  //      << search_time * 1. / CLOCKS_PER_SEC << endl;

  // for (auto it = flux_tube_schwinger.begin(); it !=
  // flux_tube_schwinger.end();
  //      it++) {
  //   cout << "T = " << get<0>(it->first) << ", R = " << get<1>(it->first)
  //        << ", d = " << get<2>(it->first)
  //        << ", longitudinal electric flux_tube schwinger = " << it->second
  //        << endl;
  // }

  // vector<vector<MATRIX_TYPE>> separated_smeared =
  //     separate_wilson(conf_smeared.array);

  // map<tuple<int, int>, double> wilson_loops =
  //     wilson_parallel(separated_smeared, R_min, R_max, T_min, T_max);

  // for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
  //   cout << "T = " << get<0>(it->first) << ", R = " << get<1>(it->first)
  //        << ", wilson_loop = " << it->second << endl;
  // }

  // map<tuple<int, int, int>, double> flux_tube_wilson;

  // start_time = omp_get_wtime();

  // vector<vector<MATRIX_TYPE>> separated_plaket =
  //     separate_wilson(conf_qc2dstag.array);

  // vector<double> plaket_time_tr = plaket_aver_tr_time(separated_plaket);
  // double plaket_time_aver = 0;
  // for (int i = 0; i < plaket_time_tr.size(); i++) {
  //   plaket_time_aver += plaket_time_tr[i];
  // }
  // plaket_time_aver = plaket_time_aver / plaket_time_tr.size();

  // cout << "plaket time " << plaket_time_aver << endl;

  // flux_tube_wilson =
  //     wilson_plaket_correlator(plaket_time_tr, separated_smeared, T_min,
  //     T_max,
  //                              R_min, R_max, 5, 0, "longitudinal");

  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // cout << "flux tube longitudinal electric time: " << search_time << endl;

  // for (auto it = flux_tube_wilson.begin(); it != flux_tube_wilson.end();
  // it++) {
  //   cout << "T = " << get<0>(it->first) << ", R = " << get<1>(it->first)
  //        << ", d = " << get<2>(it->first)
  //        << ", longitudinal electric flux_tube = " << it->second << endl;
  // }
}