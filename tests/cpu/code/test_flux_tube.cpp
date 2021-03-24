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
  string path_qc2dstag = "../../confs/qc2dstag/mu0.05/s0/CONF0201";
  string path_smeared = "../../confs/qc2dstag/HYP_APE/mu0.05/s0/smeared_0201";

  conf_qc2dstag.read_double_qc2dstag(path_qc2dstag);
  conf_smeared.read_double(path_smeared);

  int R = 14;
  int T = 8;
  int d_min = -10;
  int d_max = R;

  start_time = clock();
  vector<FLOAT> wilson_vec = calculate_wilson_loop_tr(conf_smeared.array, R, T);
  vector<FLOAT> plaket_time_vec = calculate_plaket_time_tr(conf_qc2dstag.array);
  vector<FLOAT> plaket_space_vec =
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
  end_time = clock();
  search_time = end_time - start_time;
  cout << "prepare time: " << search_time * 1. / CLOCKS_PER_SEC << endl;

  start_time = clock();
  vector<FLOAT> electric = wilson_plaket_correlator_electric(
      wilson_vec, plaket_time_vec, R, T, 0, d_min, d_max);
  end_time = clock();
  search_time = end_time - start_time;
  cout << "electric time: " << search_time * 1. / CLOCKS_PER_SEC << endl;

  start_time = clock();
  vector<FLOAT> magnetic = wilson_plaket_correlator_magnetic(
      wilson_vec, plaket_space_vec, R, T, 0, d_min, d_max);
  end_time = clock();
  search_time = end_time - start_time;

  cout << "magnetic time: " << search_time * 1. / CLOCKS_PER_SEC << endl;
  for (int i = 0; i < electric.size(); i++) {
    cout << "electric R = " << R << " T= " << T << " d " << i + d_min << " "
         << electric[i] << endl;
  }
  cout << endl;
  for (int i = 0; i < magnetic.size(); i++) {
    cout << "magnetic R = " << R << " T= " << T << " d " << i + d_min << " "
         << magnetic[i] << endl;
  }
}
