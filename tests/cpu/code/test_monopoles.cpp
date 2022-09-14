#include "../../../lib/cpu/include/Landau_U1.h"
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

using namespace std;

int x_size;
int y_size;
int z_size;
int t_size;

int main(int argc, char *argv[]) {
  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  x_size = 24;
  y_size = 24;
  z_size = 24;
  t_size = 24;

  std::cout.precision(17);

  std::string path_abelian =
      "../../confs/su2/mag/su2_suzuki/24^4/beta2.4/conf_0001";

  data<abelian> conf;
  // data<su2> conf;
  // conf.read_double_convert_abelian(path_abelian, 0);
  // conf.read_double(path_abelian, 0);
  // conf.read_float(path_abelian, 8);
  conf.read_double_convert_abelian(path_abelian, 0);
  std::vector<double> angles = convert_abelian_to_abelian(conf.array);
  // std::vector<double> angles
  // = read_angles_float_fortran(path_abelian); std::vector<double> angles =
  // read_double_qc2dstag_convet_abelian(path_abelian);
  // std::vector<double> angles = read_angles_double_fortran(path_abelian);
  // std::vector<double> angles =
  // read_float_fortran_convet_abelian(path_abelian); std::vector<double>
  // angles = read_double_fortran_convet_abelian(path_abelian);

  // std::vector<double> J = calculate_current(angles);
  std::vector<int> J = calculate_current_singular(angles);

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