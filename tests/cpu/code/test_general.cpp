#include "data.h"
#include "link.h"
#include "matrix.h"
#include "observables.h"
#include "result.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <stdio.h>
#include <vector>

using namespace std;

int x_size = 32;
int y_size = 32;
int z_size = 32;
int t_size = 32;

int main(int argc, char *argv[]) {
  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;
  link1<su2> link(x_size, y_size, z_size, t_size);
  data<su2> conf;
  data<abelian> conf_abelian;
  data<su2> conf_offd;
  char const *path1 = "../../confs/su2/time_32/mu0.00/conf_0001.fl";
  char const *path_abelian = "../../confs/su2/abelian/CON_MON_MAG_031.LAT";
  char const *path_offd = "/home/ilya/lattice/conf/offd/CON_OFF_MAG_033.LAT";
  conf.read_float(path1);
  conf_abelian.read_float_fortran(path_abelian);
  conf_offd.read_double_fortran(path_offd);

  link1<abelian> link_abelian(x_size, y_size, z_size, t_size);
  start_time = clock();

  end_time = clock();
  search_time = end_time - start_time;

  int r = 1;
  int t = 1;
  result res_plaket_time;
  result res_plaket_space;
  result res_wilson;
  result res_correlator_electric;
  result res_correlator_magnetic;
  res_plaket_time.array = calculate_plaket_time_tr(conf_abelian.array);
  res_plaket_space.array = calculate_plaket_space_tr(conf_abelian.array);
  res_wilson.array = calculate_wilson_loop_tr(conf_abelian.array, r, t);
  FLOAT aver[2];
  res_plaket_time.average(aver);

  int d_min = -10;
  int d_max = 10;
  int x_trans = 0;
  res_correlator_electric = wilson_plaket_correlator_electric(
      res_wilson.array, res_plaket_time.array, r, t, x_trans, d_min, d_max);
  res_correlator_magnetic = wilson_plaket_correlator_magnetic(
      res_wilson.array, res_plaket_space.array, r, t, x_trans, d_min, d_max);

  link1<su2> link_test(32, 32, 32, 32);
  link_test.move(2, 31);
  cout << link_test.coordinate[1] << endl;
  link_test.move(2, 10);
  cout << link_test.coordinate[1] << endl;

  cout.precision(10);

  start_time = clock();

  cout << "abelian" << endl;
  cout << "test plaket " << plaket(conf_abelian.array) / 2 << endl;
  cout << "test plaket_time " << plaket_time(conf_abelian.array) / 2 << endl;
  cout << "test plaket_space " << plaket_space(conf_abelian.array) / 2 << endl;
  cout << "test polyakov_loop " << polyakov(conf_abelian.array) / 2 << endl;
  cout << "test wilson_loop_R=10_T=6 " << wilson(conf_abelian.array, 10, 6)
       << endl;

  end_time = clock();
  search_time = end_time - start_time;
  cout << "abelian working time: " << search_time * 1. / CLOCKS_PER_SEC << endl;

  start_time = clock();

  cout << "offd" << endl;
  cout << "test plaket " << plaket(conf_offd.array) / 2 << endl;
  cout << "test plaket_time " << plaket_time(conf_offd.array) / 2 << endl;
  cout << "test plaket_space " << plaket_space(conf_offd.array) / 2 << endl;
  cout << "test polyakov_loop " << polyakov(conf_offd.array) / 2 << endl;
  cout << "test wilson_loop_R=10_T=6 " << wilson(conf_offd.array, 10, 6)
       << endl;

  end_time = clock();
  search_time = end_time - start_time;
  cout << "offd working time: " << search_time * 1. / CLOCKS_PER_SEC << endl;

  start_time = clock();

  cout << "su2" << endl;
  cout << "test plaket " << plaket(conf.array) / 2 << " right: 0.6769540066"
       << endl;
  cout << "test plaket_time " << plaket_time(conf.array) / 2
       << " right: 0.6770628794" << endl;
  cout << "test plaket_space " << plaket_space(conf.array) / 2
       << " right: 0.6768451339" << endl;
  cout << "test polyakov_loop " << polyakov(conf.array) / 2
       << " right: -0.004586235468" << endl;
  cout << "test wilson_loop_R=10_T=6 " << wilson(conf.array, 10, 6)
       << " right: 0.001178588784" << endl;

  end_time = clock();
  search_time = end_time - start_time;
  cout << "working time: " << search_time * 1. / CLOCKS_PER_SEC << endl;
}
