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
#include <stdio.h>

#define MATRIX_TYPE su3

int x_size = 64;
int y_size = 64;
int z_size = 64;
int t_size = 16;

using namespace std;

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double search_time;

  double alpha1 = 0.1;
  double alpha2 = 0.1;
  double alpha3 = 0.5;
  double alpha_APE = 0.5;
  double stout_rho = 0.15;

  data<MATRIX_TYPE> conf;
  data<MATRIX_TYPE> smeared;
  // string conf_path = "../../confs/su2/CONF0201";
  string conf_path =
      "../../confs/Coulomb_su3/QCD/140MeV/nt16/conf_Coulomb_gaugefixed_0501";
  // string conf_path = "../../confs/SU3_conf/nt14/conf.0501";
  // string conf_path = "../../confs/qc2dstag/40^4/mu0.00/CONF0201";
  // conf.read_double(path_su2, 8);
  // conf_abelian.read_float(path_abelian, 0);
  // conf.read_double_qc2dstag(conf_path);
  conf.read_ildg(conf_path);

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

  start_time = omp_get_wtime();

  for (int i = 0; i < 10; i++) {

    start_time = omp_get_wtime();

    smearing_APE_new(conf_separated, alpha_APE);
  }

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "smearing APE time: " << search_time << std::endl;

  std::cout << "smearing APE time: " << search_time << std::endl;

  std::cout << "plaket after smearing: " << plaket_parallel(conf_separated)
            << std::endl;
  std::cout << "plaket_time after smearing: "
            << plaket_time_parallel(conf_separated) << std::endl;
  std::cout << "plaket_space after smearing: "
            << plaket_space_parallel(conf_separated) << std::endl;

  conf_separated = separate_smearing(conf.array);

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