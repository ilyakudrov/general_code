#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/loop.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/monopoles.h"
#include "../../../lib/cpu/include/result.h"

#include <cstring>
#include <ctime>
#include <iostream>
#include <vector>

int x_size;
int y_size;
int z_size;
int t_size;

int main(int argc, char *argv[]) {
  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  x_size = 40;
  y_size = 40;
  z_size = 40;
  t_size = 40;

  std::cout.precision(17);

  int data_size = 4 * x_size * y_size * z_size * t_size;

  // string path_abelian = "../../confs/abelian/mu0.00/CON_32^3x32_0001.LAT";
  // string path_abelian = "../../confs/mon_wl/mu0.00/CON_MON_MAG_031.LAT";
  // string path_abelian = "../../confs/qc2dstag/mu0.05/s0/CONF0201";
  // string path_abelian =
  //     "/home/ilya/soft/lattice/decomposition/CON_MON_MAG_001.LAT";
  std::string path_abelian = "../../confs/decomposed/monopole/qc2dstag/40^4/"
                             "mu0.05/s0/conf_monopole_0202";

  // std::vector<FLOAT> angles =
  // read_float_fortran_convet_abelian(path_abelian); std::vector<FLOAT> angles
  // = read_angles_float_fortran(path_abelian); std::vector<FLOAT> angles =
  // read_double_qc2dstag_convet_abelian(path_abelian); std::vector<FLOAT>
  // angles = read_angles_double_fortran(path_abelian);
  std::vector<FLOAT> angles = read_angles_double_fortran(path_abelian);

  std::vector<FLOAT> J = calculate_current(angles);

  std::vector<loop *> LL = calculate_clusters(J);

  std::cout << "number of clusters " << LL.size() << std::endl;

  int length;

  for (int i = 0; i < LL.size(); i++) {
    length = 0;
    cluster_length(LL[i], length);
    std::vector<int> lengths_mu = {0, 0, 0, 0};
    length_mu(LL[i], lengths_mu);
  }
}