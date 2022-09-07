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

int x_size = 24;
int y_size = 24;
int z_size = 24;
int t_size = 24;

using namespace std;

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double search_time;

  double alpha1 = 0.75;
  double alpha2 = 0.6;
  double alpha3 = 0.3;
  double alpha_APE = 0.5;
  double stout_rho = 0.15;

  data<MATRIX_TYPE> conf;
  data<MATRIX_TYPE> smeared;
  // string conf_path = "../../confs/su2/CONF0201";
  string conf_path = "../../confs/MA_gauge/su3/gluodynamics/24^4/beta6.0/"
                     "CONFDP_gaugefixed_0001";
  // string conf_path = "../../confs/SU3_conf/nt14/conf.0501";
  // string conf_path = "../../confs/qc2dstag/40^4/mu0.00/CONF0201";
  // conf.read_double(path_su2, 8);
  // conf_abelian.read_float(path_abelian, 0);
  conf.read_double_qc2dstag(conf_path);
  // conf.read_ildg(conf_path);

  std::cout.precision(17);

  /*link1 link(x_size, y_size, z_size, t_size);

  MATRIX_TYPE A = conf.array[link.place];
  link.move(0, 1);
  A = A * conf.array[link.place + 1];
  link.move(1, 1);
  link.move(0, -1);
  A = A ^ conf.array[link.place];
  link.move(1, -1);
  A = conf.array[link.place + 1] + 0.5 * A;

  cout << A << endl;

  su3 B = A.proj();

  cout << B << endl;

  cout << B.determinant() << endl;
  su3 C = B ^ B;
  cout << C << endl;*/

  std::cout << "plaket before smearing: " << plaket(conf.array) << std::endl;
  std::cout << "plaket_time before smearing: " << plaket_time(conf.array)
            << std::endl;
  std::cout << "plaket_space before smearing: " << plaket_space(conf.array)
            << std::endl;
  std::cout << "polyakov_loop before smearing: " << polyakov(conf.array)
            << std::endl;

  std::vector<std::vector<MATRIX_TYPE>> conf_separated =
      separate_smearing(conf.array);

  for (int i = 0; i < 3; i++) {

    start_time = omp_get_wtime();

    smearing_APE_new(conf_separated, alpha_APE);

    end_time = omp_get_wtime();
    search_time = end_time - start_time;
    std::cout << "smearing APE time: " << search_time << std::endl;

    std::cout << "plaket after smearing: " << plaket_parallel(conf_separated)
              << std::endl;
    std::cout << "plaket_time after smearing: "
              << plaket_time_parallel(conf_separated) << std::endl;
    std::cout << "plaket_space after smearing: "
              << plaket_space_parallel(conf_separated) << std::endl;
  }

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