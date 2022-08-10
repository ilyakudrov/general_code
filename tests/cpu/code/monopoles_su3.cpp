#include "../../../lib/cpu/include/abelian_projection_su3.h"
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

using namespace std;

int main(int argc, char *argv[]) {
  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  x_size = 36;
  y_size = 36;
  z_size = 36;
  t_size = 36;

  cout.precision(17);

  int data_size = 4 * x_size * y_size * z_size * t_size;

  string path_abelian =
      "../../confs/su3/mag/36^4/beta6.3/CONFDP_gaugefixed_0001";

  data<su3> conf;
  // conf.read_double_convert_abelian(path_abelian, 8);
  conf.read_double_qc2dstag(path_abelian);
  // conf.read_double_convert_abelian(path_abelian, 0);
  vector<vector<double>> angles = make_angles_SU3(conf.array);

  for (int color = 0; color < angles.size(); color++) {

    vector<double> J = calculate_current(angles[color]);
    vector<loop *> LL = calculate_clusters(J);

    cout << "color = " << color << endl;
    cout << "number of clusters = " << LL.size() << endl;

    int length;

    map<int, int> lengths;
    map<int, int> windings;
    vector<int> lengths_mu;
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
      //     cout << "not multiple of 32: " << i << " " << length << " "
      //               << length_mu_test << " " << lengths_mu[0] << " "
      //               << lengths_mu[1] << " " << lengths_mu[2] << " "
      //               << lengths_mu[3] << endl;
      // }

      // for (int j = 0; j < 4; j++) {
      //   if (lengths_mu[j] != 0) {
      //     cout << "winding occured " << endl;
      //     windings[abs(lengths_mu[j])]++;
      //   }
      // }

      // if (length == 34) {
      //   cout << "length: " << length
      //        << " ;variation: " << cluster_variation(LL[i]) / length << endl;

      //   cout << "length: " << length
      //        << " ;number of sites: " << site_number(LL[i]) << endl;
      // }
    }

    cout << "lengths:" << endl;

    for (auto it = lengths.cbegin(); it != lengths.cend(); ++it) {
      cout << it->first << " " << it->second << "\n";
    }

    cout << endl;
    cout << "windings: " << endl;
    cout << endl;

    for (auto it = windings.begin(); it != windings.end(); ++it) {
      cout << it->first << " " << it->second << "\n";
    }
  }
}