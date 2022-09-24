#include "../../../lib/cpu/include/abelian_projection_su3.h"
#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/decomposition.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/loop.h"
#include "../../../lib/cpu/include/mag.h"
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

  x_size = 64;
  y_size = 64;
  z_size = 64;
  t_size = 6;

  cout.precision(17);

  int data_size = 4 * x_size * y_size * z_size * t_size;

  // string path_abelian =
  //     "../../confs/MA_gauge/su3/QCD/140MeV/nt6/conf.SP_gaugefixed_0501.ildg";
  string path_abelian =
      "../../confs/MA_gauge/su3/QCD/140MeV/nt6/conf.SP_gaugefixed_0501.ildg";
  // string path_abelian = "../../confs/Landau_U1xU1/gluodynamics/24^4/beta6.0/"
  // "conf_Landau_gaugefixed_0001";
  // string path_abelian = "../../confs/decomposed/monopole/gluodynamics/36^4/"
  //                       "beta6.3/conf_monopole_0001";
  // string path_abelian =
  // "../../confs/decomposed/monopoless/gluodynamics/36^4/"
  //                       "beta6.3/conf_monopoless_0001";

  data<su3> conf;
  // conf.read_double_convert_abelian(path_abelian, 8);
  // conf.read_double_qc2dstag(path_abelian);
  // conf.read_double(path_abelian, 0);
  conf.read_ildg(path_abelian);
  // conf.read_double_convert_abelian(path_abelian, 0);
  // vector<vector<double>> angles = conf.array;
  vector<vector<double>> angles = make_angles_SU3(conf.array);
  // vector<vector<double>> angles = read_double_angles_su3(path_abelian);
  // vector<vector<double>> angles =
  // read_double_su3_convet_angles(path_abelian);

  link1 link(x_size, y_size, z_size, t_size);

  std::vector<std::vector<std::vector<double>>> monopole_plakets =
      make_monopole_plakets(angles);

  for (int color = 0; color < angles.size(); color++) {

    vector<double> J =
        calculate_current_monopole_plakets(monopole_plakets[color]);
    vector<loop *> LL = calculate_clusters(J);

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

    for (int i = 0; i < LL.size(); i++) {
      length = cluster_length(LL[i]);
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

    for (auto it = lengths_wrapped.cbegin(); it != lengths_wrapped.cend();
         ++it) {
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
}