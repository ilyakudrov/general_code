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

  x_size = 16;
  y_size = 16;
  z_size = 16;
  t_size = 16;

  cout.precision(17);

  // string path_abelian =
  //     "../../confs/decomposed/monopole/gluodynamics/16^4/CON_MON_MAG_01001.LAT";
  // string path_abelian =
  // "../../confs/decomposed/monopole/gluodynamics/16^4/old/"
  //                       "CON_MON_MAG_01001.LAT";
  // string path_abelian =
  //     "/home/ilya/soft/lattice/general_code/apps/monopole_decomposition_su3/"
  //     "test/result/conf_monopole_16_1001";
  // string path_abelian =
  //     "/home/ilya/soft/lattice/general_code/apps/monopole_decomposition_su3/"
  //     "test/result/conf_monopole_16_1001_non-unitary";
  // string path_abelian =
  //     "/home/ilya/soft/lattice/general_code/apps/monopole_decomposition_su3/"
  //     "test/result/conf_monopoless_16_1001_non-unitary";
  string path_abelian =
      "/home/ilya/soft/lattice/general_code/apps/monopole_decomposition_su3/"
      "test/result/conf_monopoless_16_1001";
  // string path_abelian = "../../confs/decomposed/monopole/gluodynamics/24^4/"
  //                       "beta6.0/conf_monopole_0001";

  // data<su3_abelian> conf;
  data<su3> conf;
  // conf.read_double_convert_abelian(path_abelian, 8);
  // conf.read_double_qc2dstag(path_abelian);
  conf.read_double(path_abelian, 0);
  // conf.read_double_vitaly(path_abelian, 4);
  // conf.read_ildg(path_abelian);
  // conf.read_double_convert_abelian(path_abelian, 0);
  // vector<vector<double>> angles = conf.array;
  // vector<vector<double>> angles = make_angles_SU3(conf.array);
  // vector<vector<double>> angles = read_double_angles_su3(path_abelian);
  vector<vector<double>> angles = convert_to_angles(conf.array);
  // vector<vector<double>> angles =
  // read_double_su3_convet_angles(path_abelian);

  double sum;
  for (int i = 0; i < 10; i++) {
    sum = 0;
    for (int j = 0; j < 3; j++) {
      sum += angles[j][i];
    }
    cout << "sum = " << sum << endl;
  }

  std::cout << "qc2dstag plaket " << plaket(conf.array) << std::endl;

  cout << "abelian link" << endl << endl;

  for (int i = 0; i < 4; i++) {
    cout << angles[0][i] << " " << angles[1][i] << " " << angles[2][i] << endl;
  }

  std::vector<std::vector<std::vector<double>>> monopole_plakets =
      make_monopole_plakets(angles);

  cout << "abelian plaket" << endl << endl;

  for (int i = 0; i < 6; i++) {
    cout << monopole_plakets[0][i][0] << " " << monopole_plakets[1][i][0] << " "
         << monopole_plakets[2][i][0] << endl;
  }

  vector<double> test(3);
  for (int i = 0; i < 3; i++) {
    test[i] = 0;
    for (int j = 0; j < 6; j++) {
      for (int k = 0; k < monopole_plakets[i][j].size(); k++) {
        test[i] += monopole_plakets[i][j][k];
      }
    }
  }
  cout << "test " << test[0] << " " << test[1] << " " << test[2] << endl;

  link1 link(x_size, y_size, z_size, t_size);

  vector<vector<double>> J_test(3);
  for (int i = 0; i < 3; i++) {
    J_test[i] = calculate_current_monopole_plakets(monopole_plakets[i]);
  }

  for (int y = 0; y < 20; y++) {
    for (int x = 0; x < x_size; x++) {
      link.go_update(x, y, 0, 0);
      for (int mu = 0; mu < 4; mu++) {

        if (J_test[0][link.place + mu] > 0.3 ||
            J_test[0][link.place + mu] < -0.3 ||
            J_test[1][link.place + mu] > 0.3 ||
            J_test[1][link.place + mu] < -0.3 ||
            J_test[2][link.place + mu] > 0.3 ||
            J_test[2][link.place + mu] < -0.3) {
          cout << 1 << " " << 1 << " " << y + 1 << " " << x + 1 << " " << mu + 1
               << " " << J_test[0][link.place + mu] << " "
               << J_test[1][link.place + mu] << " "
               << J_test[2][link.place + mu] << endl;
        }
      }
    }
  }

  vector<vector<int>> J_number(3, vector<int>(4));
  for (int color = 0; color < 3; color++) {
    for (int t = 0; t < t_size; t++) {
      for (int z = 0; z < z_size; z++) {
        for (int y = 0; y < y_size; y++) {
          for (int x = 0; x < x_size; x++) {
            link.go_update(x, y, z, t);

            for (int mu = 0; mu < 4; mu++) {

              if (J_test[color][link.place + mu] < -1.3 ||
                  J_test[color][link.place + mu] > 1.3)
                J_number[color][mu] += 2;

              else if (J_test[color][link.place + mu] < -0.3 ||
                       J_test[color][link.place + mu] > 0.3)
                J_number[color][mu]++;
            }
          }
        }
      }
    }
  }

  for (int color = 0; color < 3; color++) {
    cout << J_number[color][0] << " " << J_number[color][1] << " "
         << J_number[color][2] << " " << J_number[color][3] << endl;
  }
  int J_sum = 0;
  for (int color = 0; color < 3; color++) {
    for (auto a : J_test[color]) {
      J_sum += abs(a);
    }
  }
  int J_sum1 = 0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 4; j++) {
      J_sum1 += J_number[i][j];
    }
  }
  cout << "J_sum = " << J_sum << " J_sum1 = " << J_sum1 << endl;

  for (int color = 0; color < 3; color++) {
    vector<double> J =
        calculate_current_monopole_plakets(monopole_plakets[color]);

    int J_sum = 0;

    for (auto a : J) {
      J_sum += abs(a);
    }
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

    cout << "J_sum = " << J_sum << " lengths sum =  " << lengths_sum << endl;

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