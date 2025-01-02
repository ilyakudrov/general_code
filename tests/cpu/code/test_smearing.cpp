#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/smearing.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <omp.h>

#define MATRIX_TYPE su3

int x_size = 64;
int y_size = 64;
int z_size = 64;
int t_size = 20;

using namespace std;

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double search_time;

  double alpha1 = 1;
  double alpha2 = 1;
  double alpha3 = 0.5;
  double alpha_APE = 0.5;
  double stout_rho = 0.15;

  Data::data<MATRIX_TYPE> conf;
  Data::data<MATRIX_TYPE> smeared;
  // string conf_path =
  //     "../../confs/su3/gluodynamics/16^4/beta6.0/b6p00_L16x16x16x16.01001.lime";
  string conf_path = "../../confs/su3/QCD/140MeV/nt20/conf.0501";

  string conf_format = "ildg";
  int bytes_skip = 0;
  bool convert = 0;

  get_data(conf, conf_path, conf_format, bytes_skip, convert);

  std::cout.precision(17);

  std::cout << "plaket before smearing: " << plaket(conf.array) << std::endl;
  std::cout << "plaket_time before smearing: " << plaket_time(conf.array)
            << std::endl;
  std::cout << "plaket_space before smearing: " << plaket_space(conf.array)
            << std::endl;
  std::cout << "polyakov_loop before smearing: " << polyakov_loop(conf.array)
            << std::endl;

  std::vector<std::vector<MATRIX_TYPE>> conf_separated =
      separate_smearing(conf.array);

  // start_time = omp_get_wtime();

  // for (int i = 0; i < 10; i++) {
  //   smearing_APE_parallel(conf_separated, alpha_APE);
  // }

  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // std::cout << "smearing APE time: " << search_time << std::endl;

  // std::cout << "plaket after smearing: " << plaket_parallel(conf_separated)
  //           << std::endl;
  // std::cout << "plaket_time after smearing: "
  //           << plaket_time_parallel(conf_separated) << std::endl;
  // std::cout << "plaket_space after smearing: "
  //           << plaket_space_parallel(conf_separated) << std::endl;

  // conf_separated = separate_smearing(conf.array);

  start_time = omp_get_wtime();

  smearing_HYP_new(conf_separated, alpha1, alpha2, alpha3);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "smearing HYP time: " << search_time << std::endl;

  std::cout << "plaket after smearing: " << plaket_parallel(conf_separated)
            << std::endl;
  std::cout << "plaket_time after smearing: "
            << plaket_time_parallel(conf_separated) << std::endl;
  std::cout << "plaket_space after smearing: "
            << plaket_space_parallel(conf_separated) << std::endl;
}