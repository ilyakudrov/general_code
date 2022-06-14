#include "../../../lib/cpu/include/abelian_projection_su3.h"
#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/flux_tube.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/monopoles.h"
#include "../../../lib/cpu/include/smearing.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
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
  t_size = 6;

  std::cout.precision(17);

  link1 link(x_size, y_size, z_size, t_size);
  // data<su2> conf;
  data<su3> conf;
  // data<abelian> conf;
  // std::string path_conf = "../../confs/SU3_conf/16^4/bqcd_mag.01001.lat";
  // std::string path_conf = "../../confs/SU3_conf/16^4/su3_mag_u1.01001.lat";
  std::string path_conf =
      "../../confs/SU3_conf/nt6/steps_330/conf.SP_gaugefixed_0501.ildg";

  // conf.read_double(path_conf, 4);
  // conf.read_double(path_conf);
  conf.read_ildg(path_conf);

  cout << plaket(conf.array) << endl;

  // string laplacian_path = "../../confs/inverse_laplacian/ALPHA16x16_d.LAT";

  // vector<double> inverse_laplacian = read_inverse_laplacian(laplacian_path);

  // link1 link_laplace(x_size / 2 + 1, y_size / 2 + 1, z_size / 2 + 1,
  //                    t_size / 2 + 1);

  // for (int x = 0; x < 2; x++) {
  //   for (int y = 0; y < 2; y++) {
  //     for (int z = 0; z < 2; z++) {
  //       for (int t = 0; t < 5; t++) {
  //         link_laplace.go_update(x, y, z, t);
  //         cout << x << " " << y << " " << z << " " << t << " "
  //              << inverse_laplacian[link_laplace.place / 4] << endl;
  //       }
  //     }
  //   }
  // }

  vector<vector<double>> angles = make_angles_SU3(conf.array);

  for (int i = 0; i < 10; i++) {
    cout << angles[0][i] << " " << angles[1][i] << " " << angles[2][i] << endl;
  }

  for (int color = 0; color < 3; color++) {

    cout << "color: " << color << endl;

    std::vector<double> J = calculate_current(angles[color]);
    std::vector<loop *> LL = calculate_clusters(J);

    std::cout << "number of clusters " << LL.size() << std::endl;

    int length;

    std::map<int, int> lengths;
    std::map<int, int> wrappings;
    std::map<int, int> wrapped;
    std::vector<int> lengths_mu;
    int length_mu_test;

    int total_length = 0;

    for (int i = 0; i < LL.size(); i++) {
      length = cluster_length(LL[i]);
      total_length += length;
      lengths[length]++;
      lengths_mu = length_mu(LL[i]);

      for (int j = 0; j < 4; j++) {
        if (lengths_mu[j] != 0) {
          wrappings[abs(lengths_mu[j])]++;
        }
      }

      for (int j = 0; j < 4; j++) {
        if (lengths_mu[j] != 0) {
          wrapped[length]++;
          break;
        }
      }
    }

    cout << "total length " << total_length << endl;

    cout << "lengths: " << endl;

    for (auto it = lengths.cbegin(); it != lengths.cend(); ++it) {
      std::cout << it->first << " " << it->second << "\n";
    }

    cout << "wrapped: " << endl;

    for (auto it = wrapped.cbegin(); it != wrapped.cend(); ++it) {
      std::cout << it->first << " " << it->second << "\n";
    }

    cout << "wrappings: " << endl;

    for (auto it = wrappings.cbegin(); it != wrappings.cend(); ++it) {
      std::cout << it->first << " " << it->second << "\n";
    }
  }

  // std::vector<double> monopole_angles =
  //     make_monopole_angles(angles[0], inverse_laplacian);

  // for (int i = 0; i < 10; i++) {
  //   cout << monopole_angles[i] << endl;
  // }
}