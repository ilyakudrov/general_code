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
int size1;
int size2;

using namespace std;

void J_places(vector<vector<vector<double>>> monopole_plakets) {
  vector<vector<double>> J_test(3);
  for (int i = 0; i < 3; i++) {
    J_test[i] = calculate_current_monopole_plakets(monopole_plakets[i]);
  }

  link1 link(x_size, y_size, z_size, t_size);

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
}

void J_number(vector<vector<vector<double>>> monopole_plakets) {
  vector<vector<double>> J_test(3);
  for (int i = 0; i < 3; i++) {
    J_test[i] = calculate_current_monopole_plakets(monopole_plakets[i]);
  }

  link1 link(x_size, y_size, z_size, t_size);

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
}

template <class T>
tuple<double, T> aver_and_max_plaket(vector<vector<vector<T>>> plakets) {
  double aver = 0;
  T max = 0;
  T sum_tmp;
  int size = plakets[0][0].size();
  for (int j = 0; j < 6; j++) {
    for (int k = 0; k < size; k++) {
      sum_tmp = 0;
      for (int i = 0; i < 3; i++) {
        sum_tmp += plakets[i][j][k];
      }
      aver += sum_tmp;
      if (abs(sum_tmp) > max) {
        max = abs(sum_tmp);
      }
    }
  }
  aver /= size * 6.;

  return tuple<double, T>({aver, max});
}

int main(int argc, char *argv[]) {
  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  x_size = 16;
  y_size = 16;
  z_size = 16;
  t_size = 16;
  size1 = x_size * y_size;
  size2 = x_size * y_size * z_size;

  cout.precision(17);

  // string path_abelian =
  // "../../../apps/monopole_decomposition_su3/test/result/"
  //                       "conf_monopole_16_1001_compensate";
  // string path_abelian =
  // "../../../apps/monopole_decomposition_su3/test/result/"
  //                       "conf_monopole_16_1001_nocompensate";
  // string path_abelian =
  //     "../../confs/SU3_conf/gluodynamics/16^4/su3_mag_u1.01001.lat";
  string path_abelian = "../../../tests/confs/decomposed/monopole/gluodynamics/"
                        "36^4/beta6.3/steps_25/copies=4/conf_monopole_0001";

  // data<su3_abelian> conf;
  Data::data<su3> conf;
  // conf.read_double_convert_abelian(path_abelian, 8);
  // conf.read_double_qc2dstag(path_abelian);
  conf.read_double(path_abelian, 0);
  // conf.read_double_vitaly(path_abelian, 4);
  // conf.read_ildg(path_abelian);
  // conf.read_double_convert_abelian(path_abelian, 0);
  // vector<vector<double>> angles = conf.array;
  vector<vector<double>> angles = make_angles_SU3(conf.array);
  // vector<vector<double>> angles = read_double_angles_su3(path_abelian);
  // vector<vector<double>> angles = convert_to_angles(conf.array);
  // vector<vector<double>> angles =
  // read_double_su3_convet_angles(path_abelian);

  // double sum;
  // for (int i = 0; i < 10; i++) {
  //   sum = 0;
  //   for (int j = 0; j < 3; j++) {
  //     sum += angles[j][i];
  //   }
  //   cout << "sum = " << sum << endl;
  // }

  cout << "qc2dstag plaket " << plaket(conf.array) << endl;

  // vector<vector<vector<double>>> monopole_plakets =
  //     make_monopole_plakets(angles);
  // vector<vector<vector<int>>> dirac_plakets =
  //     make_monopole_plakets_singular(angles);

  vector<vector<vector<double>>> monopole_plakets(3);
  vector<vector<vector<int>>> dirac_plakets(3);

  make_plakets_both(angles, monopole_plakets, dirac_plakets);

  tuple<double, double> abelian_aver = aver_and_max_plaket(monopole_plakets);
  tuple<double, double> dirac_aver = aver_and_max_plaket(dirac_plakets);

  cout << "abelian_aver " << get<0>(abelian_aver) << " " << get<1>(abelian_aver)
       << endl;
  cout << "dirac_aver " << get<0>(dirac_aver) << " " << get<1>(dirac_aver)
       << endl;

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

    cout << endl << "asymmetry " << asymmetry << endl;
  }
}