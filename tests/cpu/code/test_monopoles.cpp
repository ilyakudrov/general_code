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
      "../../confs/MA_gauge/su2/su2_suzuki/24^4/beta2.4/CON_fxd_MAG_001.LAT";

  data<abelian> conf;
  // data<su2> conf;
  // conf.read_double_convert_abelian(path_abelian, 0);
  // conf.read_double(path_abelian, 0);
  // conf.read_float(path_abelian, 8);
  conf.read_double_convert_abelian(path_abelian, 8);
  std::vector<double> angles = convert_to_angles(conf.array);
  // std::vector<double> angles
  // = read_angles_float_fortran(path_abelian); std::vector<double> angles =
  // read_double_qc2dstag_convet_abelian(path_abelian);
  // std::vector<double> angles = read_angles_double_fortran(path_abelian);
  // std::vector<double> angles =
  // read_float_fortran_convet_abelian(path_abelian); std::vector<double>
  // angles = read_double_fortran_convet_abelian(path_abelian);

  // std::vector<double> J = calculate_current(angles);
  std::vector<int> J = calculate_current_singular(angles);

  link1 link(x_size, y_size, z_size, t_size);

  for (int y = 0; y < 20; y++) {
    for (int x = 0; x < x_size; x++) {
      link.go_update(x, y, 0, 0);
      for (int mu = 0; mu < 4; mu++) {

        if (J[link.place + mu] > 0.3 || J[link.place + mu] < -0.33) {
          cout << 1 << " " << 1 << " " << y + 1 << " " << x + 1 << " " << mu + 1
               << " " << J[link.place + mu] << endl;
        }
      }
    }
  }

  vector<int> J_number(4);
  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
          link.go_update(x, y, z, t);

          for (int mu = 0; mu < 4; mu++) {

            if (J[link.place + mu] < -1.3 || J[link.place + mu] > 1.3)
              J_number[mu] += 2;

            else if (J[link.place + mu] < -0.3 || J[link.place + mu] > 0.3)
              J_number[mu]++;
          }
        }
      }
    }
  }

  cout << J_number[0] << " " << J_number[1] << " " << J_number[2] << " "
       << J_number[3] << endl;

  std::vector<loop *> LL = calculate_clusters(J);

  std::cout << "number of clusters " << LL.size() << std::endl;

  int length;

  map<int, int> lengths_unwrapped;
  map<int, int> lengths_wrapped;
  map<int, int> space_windings;
  map<int, int> time_windings;
  vector<int> lengths_mu;
  int length_mu_test;
  vector<int> currents;

  int space_currents = 0;
  int time_currents = 0;

  int lengths_sum = 0;

  for (int i = 0; i < LL.size(); i++) {
    length = cluster_length(LL[i]);
    lengths_sum += length;

    lengths_mu = length_mu(LL[i]);

    currents = currents_directions(LL[i]);

    space_currents += currents[0];
    time_currents += currents[1];

    for (int j = 0; j < 3; j++) {
      if (lengths_mu[j] != 0) {
        space_windings[abs(lengths_mu[j]) / x_size]++;
      }
    }

    if (lengths_mu[3] != 0) {
      time_windings[abs(lengths_mu[3]) / t_size]++;
    }

    if (lengths_mu[3] == 0)
      lengths_unwrapped[length]++;
    else if (lengths_mu[3] != 0)
      lengths_wrapped[length]++;
  }

  cout << "unwrapped clusters:" << endl << endl;

  for (auto it = lengths_unwrapped.cbegin(); it != lengths_unwrapped.cend();
       ++it) {
    cout << it->first << "," << it->second << endl;
  }

  cout << endl << "wrapped clusters:" << endl << endl;

  for (auto it = lengths_wrapped.cbegin(); it != lengths_wrapped.cend(); ++it) {
    cout << it->first << "," << it->second << endl;
  }

  cout << endl << "time windings:" << endl << endl;

  for (auto it = time_windings.begin(); it != time_windings.end(); ++it) {
    cout << it->first << "," << it->second << ",time" << endl;
  }

  cout << endl << "space windings:" << endl << endl;

  for (auto it = space_windings.begin(); it != space_windings.end(); ++it) {
    cout << it->first << "," << it->second << ",space" << endl;
  }

  double asymmetry = (space_currents / 3. - time_currents) /
                     (space_currents / 3. + time_currents);

  cout << endl << asymmetry << endl;
}