#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/smearing.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <omp.h>
#include <vector>

int x_size;
int y_size;
int z_size;
int t_size;

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

  std::cout.precision(17);

  data<MATRIX_TYPE> conf;
  string conf_path =
      "../../confs/Coulomb_su3/QCD/140MeV/nt4/conf_Coulomb_gaugefixed_0501";
  // string conf_path = "../../confs/smeared/QCD/140MeV/nt4/smeared_0501";
  // conf.read_double(conf_path, 0);
  // conf.read_double_qc2dstag(conf_path);
  conf.read_ildg(conf_path);
  // conf.read_float(conf_path, 8);
  // conf.read_double_convert_abelian(conf_path, 8);
  // conf.read_double_qc2dstag_convert_abelian(conf_path);

  // plakets and polyakov loop
  start_time = clock();
  std::cout << "plaket " << plaket(conf.array) << std::endl;
  std::cout << "plaket_time " << plaket_time(conf.array) << std::endl;
  std::cout << "plaket_space " << plaket_space(conf.array) << std::endl;
  std::cout << "polyakov " << polyakov_loop(conf.array) << std::endl;
  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "plaket and staff time: " << search_time * 1. / CLOCKS_PER_SEC
            << std::endl;

  std::vector<std::vector<MATRIX_TYPE>> conf_separated =
      separate_wilson(conf.array);

  int D_max = 8;

  start_time = omp_get_wtime();

  //   std::vector<double> polyakov_correlator_vec =
  //       polyakov_loop_correlator_singlet(conf.array, D_max);
  std::vector<double> polyakov_correlator_vec =
      polyakov_loop_correlator(conf.array, D_max);

  std::map<double, double> polyakov_correlator =
      polyakov_average_directions(polyakov_correlator_vec, D_max);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "polyakov_loop_correlator_singlet time: " << search_time
            << std::endl;

  for (auto it = polyakov_correlator.begin(); it != polyakov_correlator.end();
       it++) {
    std::cout << "D = " << it->first << ", correlator = " << it->second
              << std::endl;
  }
}
