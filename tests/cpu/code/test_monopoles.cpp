#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/link.h"
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

  // string path_abelian = "../../confs/su2/abelian/CON_MON_MAG_031.LAT";
  string path_abelian = "../../confs/abelian/mu0.00/CON_32^3x32_0001.LAT";

  vector<FLOAT> angles = read_angles_float_fortran(path_abelian);

  vector<FLOAT> J = calculate_current(angles);

  // for (int i = 0; i < 500; i++) {
  //   if (J[i] > 0.3 || J[i] < -0.3) {
  //     cout << i << endl;
  //   }
  // }

  // for (int i = 0; i < 100; i++) {
  //   cout << J[i] << endl;
  // }

  vector<loop *> LL;
  vector<FLOAT> J_test(data_size);
  link1 link_test(x_size, y_size, z_size, t_size);

  // first square
  J_test[link_test.place] = 1;
  link_test.move(0, 1);
  J_test[link_test.place + 1] = 1;
  link_test.move(1, 1);
  link_test.move(0, -1);
  J_test[link_test.place] = -1;
  link_test.move(1, -1);
  J_test[link_test.place + 1] = -1;

  link_test.move(0, 1);
  link_test.move(1, 1);

  // second square
  J_test[link_test.place] = 1;
  link_test.move(0, 2);
  J_test[link_test.place + 2] = 1;
  link_test.move(2, 1);
  link_test.move(0, -1);
  J_test[link_test.place] = -1;
  link_test.move(2, -1);
  J_test[link_test.place + 2] = -1;

  // and another one

  link_test.move(0, 5);
  link_test.move(1, 6);
  link_test.move(2, 7);

  // first square
  J_test[link_test.place] = 1;
  link_test.move(0, 1);
  J_test[link_test.place + 1] = 1;
  link_test.move(1, 1);
  link_test.move(0, -1);
  J_test[link_test.place] = -1;
  link_test.move(1, -1);
  J_test[link_test.place + 1] = -1;

  link_test.move(0, 1);
  link_test.move(1, 1);

  // second square
  J_test[link_test.place] = 1;
  link_test.move(0, 2);
  J_test[link_test.place + 2] = 1;
  link_test.move(2, 1);
  link_test.move(0, -1);
  J_test[link_test.place] = -1;
  link_test.move(2, -1);
  J_test[link_test.place + 2] = -1;

  // cout << J[] << endl;

  link1 link_test1(x_size, y_size, z_size, t_size);
  link_test1.go(24, 1, 0, 0);
  link_test1.update(0);
  link_test1.update(1);

  cout << "dir " << find_current(link_test1, J) << endl;
  link_test1.move(1, 1);
  cout << "dir " << find_current(link_test1, J) << endl;

  calculate_clusters(LL, J);

  cout << LL.size() << endl;

  int length;

  // for (int i = 0; i < 10 /*LL.size()*/; i++) {
  //   length = 0;
  //   cluster_length(*LL[i], length);
  //   cout << length << endl;
  // }
}