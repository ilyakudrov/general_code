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
  link1<su2> link(x_size, y_size, z_size, t_size);
  data<su2> conf;
  data<abelian> conf_abelian;
  data<su2> conf_offd;
  data<su2> conf_qc2dstag;
  char const *path1 = "../../confs/su2/time_32/mu0.00/conf_0001.fl";
  char const *path_abelian = "../../confs/su2/abelian/CON_MON_MAG_031.LAT";
  char const *path_offd = "/home/ilya/lattice/conf/offd/CON_OFF_MAG_033.LAT";
  char const *path_qc2dstag =
      "/home/ilya/lattice/general_code/tests/confs/qc2dstag/mu0.05/CONF0001";
  conf.read_float(path1);
  conf_abelian.read_float_fortran(path_abelian);
  conf_offd.read_double_fortran(path_offd);
  x_size = 40;
  y_size = 40;
  z_size = 40;
  t_size = 40;
  conf_qc2dstag.read_double_qc2dstag(path_qc2dstag);

  cout << "qc2dstag plaket " << plaket(conf_qc2dstag.array) / 2 << endl;
  cout << "qc2dstag plaket_time " << plaket_time(conf_qc2dstag.array) / 2
       << endl;
  cout << "qc2dstag plaket_space " << plaket_space(conf_qc2dstag.array) / 2
       << endl;
  cout << "qc2dstag wilson " << wilson(conf_qc2dstag.array, 1, 1) / 2 << endl;
  cout << "qc2dstag wilson " << wilson(conf_qc2dstag.array, 10, 10) / 2 << endl;
  cout << "qc2dstag polyakov " << polyakov(conf_qc2dstag.array) / 2 << endl;
  cout << "module " << conf_qc2dstag.array[0].module() << endl;

  link1<su2> link_qc2dstag(x_size, y_size, z_size, t_size);
  link_qc2dstag.move(1, 2);
  cout << "qc2dstag plaket_mu "
       << link_qc2dstag.plaket_mu(conf_qc2dstag.array, 2) << endl;

  x_size = 32;
  y_size = 32;
  z_size = 32;
  t_size = 32;
  link1<abelian> link_abelian(x_size, y_size, z_size, t_size);
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
