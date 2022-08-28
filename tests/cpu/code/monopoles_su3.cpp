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

  x_size = 24;
  y_size = 24;
  z_size = 24;
  t_size = 24;

  cout.precision(17);

  int data_size = 4 * x_size * y_size * z_size * t_size;

  // string path_abelian = "../../confs/su3/gluodynamics/36^4/beta6.3/"
  //                       "CONF0001";
  // string path_abelian = "../../confs/su3/mag/gluodynamics/36^4/beta6.3/"
  //                       "CONFDP_gaugefixed_0001";
  // string path_abelian =
  //     "../../confs/su3/Landau_U1xU1/gluodynamics/32^4/beta6.2/"
  //     "conf_Landau_gaugefixed_0001";

  // string path_abelian =
  //     "../../confs/su3/140MeV/nt6/conf.SP_gaugefixed_0501.ildg";
  // string path_abelian = "../../confs/su3/140MeV/nt6/conf.0501";

  // string path_abelian =
  //     "/home/ilya/soft/lattice/general_code/apps/monopole_decomposition_su3/"
  //     "test/result/conf_monopole_24_0001";
  string path_abelian =
      "/home/ilya/soft/lattice/general_code/apps/monopole_decomposition_su3/"
      "test/result/conf_monopoless_24_0001";

  data<su3> conf;
  // conf.read_double_convert_abelian(path_abelian, 8);
  // conf.read_double_qc2dstag(path_abelian);
  // conf.read_double(path_abelian, 0);
  // conf.read_ildg(path_abelian);
  // conf.read_double_convert_abelian(path_abelian, 0);
  // vector<vector<double>> angles = conf.array;
  // vector<vector<double>> angles = make_angles_SU3(conf.array);
  // vector<vector<double>> angles = read_double_angles_su3(path_abelian);
  // angles_project(angles);
  vector<vector<double>> angles = read_double_su3_convet_angles(path_abelian);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      cout << angles[i][j] << endl;
    }
  }

  link1 link(x_size, y_size, z_size, t_size);

  // su3 lambda3;
  // lambda3.matrix[0][0] = complex_t(1, 0);
  // lambda3.matrix[1][1] = complex_t(-1, 0);
  // lambda3.matrix[2][2] = complex_t(0, 0);

  // cout << "lambda3" << endl;
  // cout << lambda3 << endl;

  // su3 lambda8;
  // lambda8.matrix[0][0] = complex_t(1. / (2 * sqrt(3)), 0);
  // lambda8.matrix[1][1] = complex_t(1. / (2 * sqrt(3)), 0);
  // lambda8.matrix[2][2] = complex_t(-2. / (2 * sqrt(3)), 0);

  // cout << "lambda8" << endl;
  // cout << lambda8 << endl;

  // cout << "unity check1 " << conf.array[0] * conf.array[0].conj() << endl;
  // cout << "unity check2 " << (conf.array[0] ^ conf.array[0]) << endl;

  // double sum = 0;
  // for (int i = 0; i < 3; i++) {
  //   sum += conf.array[0].matrix[i][i].real * conf.array[0].matrix[i][i].real
  //   +
  //          conf.array[0].matrix[i][i].imag * conf.array[0].matrix[i][i].imag;
  // }
  // cout << "sum of diag squares " << sum << endl;

  // cout << "lambda sum tr "
  //      << (lambda3 * conf.array[0] * lambda3 * conf.array[0].conj()).tr() +
  //      (lambda8 * conf.array[0] * lambda8 * conf.array[0].conj()).tr()
  //      << endl;
  // cout << "mult " << (lambda8 * conf.array[0]) * lambda8 << endl;

  // cout << "mag_su3 functional " << mag_functional_su3(conf.array) << endl;

  for (int color = 0; color < 3; color++) {

    vector<double> J = calculate_current(angles[color]);
    // vector<int> J = calculate_current_singular(angles[color]);
    vector<loop *> LL = calculate_clusters(J);

    cout << "color = " << color << endl;
    cout << "number of clusters = " << LL.size() << endl;

    int length;

    map<int, int> lengths;
    map<int, int> windings;
    vector<int> lengths_mu;
    int length_mu_test;
    vector<int> currents;

    int space_length = 0;
    int time_length = 0;

    for (int i = 0; i < LL.size(); i++) {
      length = cluster_length(LL[i]);
      lengths[length]++;
      lengths_mu = length_mu(LL[i]);

      currents = currents_directions(LL[i]);

      space_length += currents[0];
      time_length += currents[1];

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

    cout << "length ratio " << 1. * space_length / time_length / 3 << endl;
    cout << endl;
  }
}