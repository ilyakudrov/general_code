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

  x_size = 32;
  y_size = 32;
  z_size = 32;
  t_size = 32;

  cout.precision(17);

  int data_size = 4 * x_size * y_size * z_size * t_size;

  string path_abelian = "../../confs/abelian/mu0.00/CON_32^3x32_0001.LAT";
  // string path_abelian = "../../confs/mon_wl/mu0.00/CON_MON_MAG_031.LAT";
  // string path_abelian = "../../confs/qc2dstag/mu0.05/s0/CONF0201";

  vector<FLOAT> angles = read_float_fortran_convet_abelian(path_abelian);
  // vector<FLOAT> angles = read_angles_float_fortran(path_abelian);
  // vector<FLOAT> angles = read_double_qc2dstag_convet_abelian(path_abelian);

  vector<FLOAT> J = calculate_current(angles);

  vector<loop *> LL = calculate_clusters(J);

  cout << "number of clusters " << LL.size() << endl;

  int length;

  for (int i = 0; i < LL.size(); i++) {
    length = 0;
    cluster_length(LL[i], length);
    vector<int> lengths_mu = {0, 0, 0, 0};
    length_mu(LL[i], lengths_mu);
  }
}