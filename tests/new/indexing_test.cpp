#include "../../lib/cpu/include/basic_observables.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/flux_tube.h"
#include "../../lib/cpu/include/link.h"
#include "../../lib/cpu/include/matrix.h"
#include "../../lib/cpu/include/smearing.h"
#include "include/indexing.h"

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

using namespace std;

int main(int argc, char *argv[]) {
  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  x_size = 64;
  y_size = 64;
  z_size = 64;
  t_size = 14;

  std::cout.precision(17);

  link1 link(x_size, y_size, z_size, t_size);
  // data<su2> conf;
  // data<su3_full> conf;
  data<su3> conf;
  // data<abelian> conf;
  // string conf_path = "../confs/qc2dstag/40^4/mu0.05/s0/CONF0201";
  string conf_path = "../confs/SU3_conf/nt14/conf.0501";
  // conf.read_double(conf_path, 8);
  // conf.read_double_qc2dstag(conf_path);
  conf.read_ildg(conf_path);

  // plakets and polyakov loop
  start_time = clock();

  std::cout << "qc2dstag plaket_time " << plaket(conf.array) << std::endl;

  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "plaket time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;

  std::vector<std::vector<su3>> separated = separate_3(conf.array);
  // std::vector<std::vector<su3_full>> separated = separate_3(conf.array, 3);
  // std::vector<std::vector<su2>> separated = separate_3(conf.array, 3);

  // start_time = clock();

  // std::cout << "qc2dstag plaket_time " << plaket_time_test(separated)
  //           << std::endl;

  // end_time = clock();
  // search_time = end_time - start_time;
  // std::cout << "plaket new time: " << search_time * 1. / CLOCKS_PER_SEC
  //           << std::endl;

  // start_time = clock();

  // std::cout << "qc2dstag plaket_time " << plaket_time_test1(separated)
  //           << std::endl;

  // end_time = clock();
  // search_time = end_time - start_time;
  // std::cout << "plaket1 new time: " << search_time * 1. / CLOCKS_PER_SEC
  //           << std::endl;

  // start_time = clock();

  // std::cout << "qc2dstag plaket_time " << plaket_time_test3(separated)
  //           << std::endl;

  // end_time = clock();
  // search_time = end_time - start_time;
  // std::cout << "plaket3 new time: " << search_time * 1. / CLOCKS_PER_SEC
  //           << std::endl;

  // start_time = clock();

  // std::cout << "qc2dstag plaket_time " << plaket_time_test4(separated)
  //           << std::endl;

  // end_time = clock();
  // search_time = end_time - start_time;
  // std::cout << "plaket4 new time: " << search_time * 1. / CLOCKS_PER_SEC
  //           << std::endl;

  // start_time = clock();

  // std::cout << "qc2dstag plaket_time " << plaket_time_test5(conf.array)
  //           << std::endl;

  // end_time = clock();
  // search_time = end_time - start_time;
  // std::cout << "plaket5 new time: " << search_time * 1. / CLOCKS_PER_SEC
  //           << std::endl;

  // start_time = clock();

  // std::cout << "qc2dstag plaket_time " << plaket_time_test6(conf.array)
  //           << std::endl;

  // end_time = clock();
  // search_time = end_time - start_time;
  // std::cout << "plaket6 new time: " << search_time * 1. / CLOCKS_PER_SEC
  //           << std::endl;

  // start_time = clock();

  // std::cout << "qc2dstag plaket_time " << plaket_time_test7(conf.array)
  //           << std::endl;

  // end_time = clock();
  // search_time = end_time - start_time;
  // std::cout << "plaket7 new time: " << search_time * 1. / CLOCKS_PER_SEC
  //           << std::endl;

  // start_time = clock();

  // std::cout << "qc2dstag plaket_time " << plaket_time_test8(separated)
  //           << std::endl;

  // end_time = clock();
  // search_time = end_time - start_time;
  // std::cout << "plaket8 new time: " << search_time * 1. / CLOCKS_PER_SEC
  //           << std::endl;

  // start_time = clock();

  // std::cout << "qc2dstag plaket " << plaket(conf.array) << std::endl;

  // end_time = clock();
  // search_time = end_time - start_time;
  // std::cout << "plaket_full time: " << search_time * 1. / CLOCKS_PER_SEC
  //           << std::endl;

  // start_time = clock();

  // std::cout << "qc2dstag plaket_space " << plaket_space(conf.array)
  //           << std::endl;

  // end_time = clock();
  // search_time = end_time - start_time;
  // std::cout << "plaket_space time: " << search_time * 1. / CLOCKS_PER_SEC
  //           << std::endl;

  // start_time = clock();

  // std::cout << "qc2dstag plaket_time_test9 " << plaket_time_test9(separated,
  // 32)
  //           << std::endl;

  // end_time = clock();
  // search_time = end_time - start_time;
  // std::cout << "plaket_time_test9 time: " << search_time * 1. /
  // CLOCKS_PER_SEC
  //           << std::endl;

  // start_time = omp_get_wtime();

  // std::cout << "qc2dstag plaket_time_test10 "
  //           << plaket_time_test10(separated, 32) << std::endl;
  std::cout << plaket_time_test10(separated, 4) << std::endl;

  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // std::cout << "plaket_time_test10 time: " << search_time * 1. /
  // CLOCKS_PER_SEC
  //           << std::endl;
}
