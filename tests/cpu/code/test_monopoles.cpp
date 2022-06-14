#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/loop.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/monopoles.h"

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

  x_size = 32;
  y_size = 32;
  z_size = 32;
  t_size = 32;

  std::cout.precision(17);

  int data_size = 4 * x_size * y_size * z_size * t_size;

  // std::string path_abelian =
  // "../../confs/abelian/mu0.00/CON_32^3x32_0001.LAT"; string path_abelian =
  // "../../confs/mon_wl/mu0.00/CON_MON_MAG_031.LAT";
  // std::string path_abelian =
  //     "../../confs/MA_gauge/qc2dstag/40^4/mu0.05/s0/conf_abelian_0201";
  std::string path_abelian =
      "../../confs/su2_suzuki/32^4/beta2.6/CON_fxd_MAG_001.LAT";

  // string path_abelian =
  //     "/home/ilya/soft/lattice/decomposition/CON_MON_MAG_001.LAT";
  // std::string path_abelian =
  // "/home/ilya/soft/lattice/decomposition/test/confs/"
  //                            "monopole/40^4/conf_monopole_0201";
  // std::string path_abelian = "../../confs/decomposed/monopole/qc2dstag/40^4/"
  //                            "mu0.05/s0/conf_monopole_0203";
  // std::string path_abelian =
  // "../../confs/decomposed/monopoless/qc2dstag/40^4/"
  //                            "mu0.05/s0/conf_monopoless_0201";

  data<abelian> conf;
  conf.read_double_convert_abelian(path_abelian, 8);
  std::vector<double> angles = convert_abelian_to_abelian(conf.array);
  // std::vector<double> angles
  // = read_angles_float_fortran(path_abelian); std::vector<double> angles =
  // read_double_qc2dstag_convet_abelian(path_abelian);
  // std::vector<double> angles = read_angles_double_fortran(path_abelian);
  // std::vector<double> angles =
  // read_float_fortran_convet_abelian(path_abelian); std::vector<double>
  // angles = read_double_fortran_convet_abelian(path_abelian);

  std::vector<double> J = calculate_current(angles);
  std::vector<loop *> LL = calculate_clusters(J);

  std::cout << "number of clusters " << LL.size() << std::endl;

  int length;

  std::map<int, int> lengths;
  std::map<int, int> windings;
  std::vector<int> lengths_mu;
  int length_mu_test;

  for (int i = 0; i < LL.size(); i++) {
    length = cluster_length(LL[i]);
    lengths[length]++;
    lengths_mu = length_mu(LL[i]);

    // for (int j = 0; j < 4; j++) {
    //   length_mu_test = lengths_mu[j];
    //   while (length_mu_test < 0 || length_mu_test > 31) {
    //     if (length_mu_test < 0)
    //       length_mu_test += 32;
    //     if (length_mu_test > 31)
    //       length_mu_test -= 32;
    //   }
    //   if (length_mu_test != 0)
    //     std::cout << "not multiple of 32: " << i << " " << length << " "
    //               << length_mu_test << " " << lengths_mu[0] << " "
    //               << lengths_mu[1] << " " << lengths_mu[2] << " "
    //               << lengths_mu[3] << std::endl;
    // }

    for (int j = 0; j < 4; j++) {
      if (lengths_mu[j] != 0) {
        std::cout << "winding occured " << std::endl;
        windings[abs(lengths_mu[j])]++;
      }
    }

    if (length == 34) {
      std::cout << "length: " << length
                << " ;variation: " << cluster_variation(LL[i]) / length
                << std::endl;

      std::cout << "length: " << length
                << " ;number of sites: " << site_number(LL[i]) << std::endl;
    }
  }

  for (auto it = lengths.cbegin(); it != lengths.cend(); ++it) {
    std::cout << it->first << " " << it->second << "\n";
  }

  std::cout << std::endl;
  std::cout << "windings: " << std::endl;
  std::cout << std::endl;

  for (auto it = windings.begin(); it != windings.end(); ++it) {
    std::cout << it->first << " " << it->second << "\n";
  }
}