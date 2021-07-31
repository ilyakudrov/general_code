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
#include <map>
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

  // std::string path_abelian =
  // "../../confs/abelian/mu0.00/CON_32^3x32_0001.LAT"; string path_abelian =
  // "../../confs/mon_wl/mu0.00/CON_MON_MAG_031.LAT"; string path_abelian =
  // "../../confs/qc2dstag/mu0.05/s0/CONF0201"; string path_abelian =
  //     "/home/ilya/soft/lattice/decomposition/CON_MON_MAG_001.LAT";
  // std::string path_abelian =
  // "/home/ilya/soft/lattice/decomposition/test/confs/"
  //                            "monopole/40^4/conf_monopole_0201";
  std::string path_abelian = "../../confs/decomposed/monopole/qc2dstag/40^4/"
                             "mu0.05/s0/conf_monopole_0201";
  // std::string path_abelian =
  // "../../confs/decomposed/monopoless/qc2dstag/40^4/"
  //                            "mu0.05/s0/conf_monopoless_0201";

  // std::vector<FLOAT> angles =
  // read_float_fortran_convet_abelian(path_abelian); std::vector<FLOAT> angles
  // = read_angles_float_fortran(path_abelian); std::vector<FLOAT> angles =
  // read_double_qc2dstag_convet_abelian(path_abelian); std::vector<FLOAT>
  std::vector<FLOAT> angles = read_angles_double_fortran(path_abelian);
  // std::vector<FLOAT> angles =
  // read_float_fortran_convet_abelian(path_abelian); std::vector<FLOAT> angles
  // = read_double_fortran_convet_abelian(path_abelian);

  std::vector<FLOAT> J = calculate_current(angles);
  // for (int i = 0; i < J.size(); i++) {
  //   if (J[i] > 1.3 || J[i] < -1.3)
  //     std::cout << J[i] << " " << std::endl;
  // }
  std::vector<FLOAT> J_test = calculate_current(angles);

  std::vector<loop *> LL = calculate_clusters(J);

  std::cout << "number of clusters " << LL.size() << std::endl;

  // std::cout << LL[0]->coordinate[0] << " " << LL[0]->coordinate[1] << " "
  //           << LL[0]->coordinate[2] << " " << LL[0]->coordinate[3] <<
  //           std::endl;

  // cluster_sites(LL[0]);

  // bool include;
  // int coordinate[4] = {22, 0, 1, 0};
  // for (int i = 0; i < LL.size(); i++) {
  //   include = false;
  //   check_for_coordinate(LL[i], coordinate, include);
  //   if (include)
  //     std::cout << "coordinate included in cluster " << i << std::endl;
  // }

  int length;
  int bool_test = true;

  // link1 link(x_size, y_size, z_size, t_size);

  // link.go_update(21, 0, 0, 0);
  // std::cout << J_test[link.place] << std::endl;
  // link.go_update(21, 0, 0, 0);
  // std::cout << J_test[link.place + 2] << std::endl;
  // link.go_update(22, 0, 1, 31);
  // std::cout << J_test[link.place + 3] << std::endl;
  // link.go_update(22, 0, 0, 0);
  // std::cout << J_test[link.place + 2] << std::endl;

  std::map<int, int> lengths;
  std::vector<int> lengths_mu;
  int length_mu_test;

  for (int i = 0; i < LL.size(); i++) {
    length = cluster_length(LL[i]);
    lengths[length]++;
    lengths_mu = length_mu(LL[i]);

    for (int j = 0; j < 4; j++) {
      length_mu_test = lengths_mu[j];
      while (length_mu_test < 0 || length_mu_test > 31) {
        if (length_mu_test < 0)
          length_mu_test += 32;
        if (length_mu_test > 31)
          length_mu_test -= 32;
      }
      if (length_mu_test != 0)
        std::cout << "not multiple of 32: " << i << " " << length << " "
                  << length_mu_test << " " << lengths_mu[0] << " "
                  << lengths_mu[1] << " " << lengths_mu[2] << " "
                  << lengths_mu[3] << std::endl;
    }
  }

  for (auto it = lengths.cbegin(); it != lengths.cend(); ++it) {
    std::cout << it->first << " " << it->second << "\n";
  }
}