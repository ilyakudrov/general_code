#include "../include/data.h"
#include "../include/link.h"
#include "../include/matrix.h"
#include "../include/observables.h"
#include "../include/result.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <stdio.h>

using namespace std;

int x_size = 32;
int y_size = 32;
int z_size = 32;
int t_size = 32;

int main(int argc, char *argv[]) {
  FLOAT alpha1 = 1;
  FLOAT alpha2 = 1;
  FLOAT alpha3 = 0.5;
  FLOAT alpha_APE = 0.3;
  FLOAT stout_rho = 0.15;

  vector<vector<su2>> smearing_first(6, vector<su2>(1));
  vector<vector<su2>> smearing_second(4, vector<su2>(1));

  link1<su2> link(x_size, y_size, z_size, t_size);
  link1<abelian> link_double(x_size, y_size, z_size, t_size);
  data<su2> conf;
  data<su2> smeared;
  data<abelian> conf_abelian;
  data<abelian> conf_abelian_smeared;
  char const *path1 = "../../confs/su2/time_32/mu0.00/conf_0001.fl";
  char const *path_abelian = "../../confs/su2/abelian/CON_MON_MAG_031.LAT";
  conf.read_float(path1);
  conf_abelian.read_float_fortran(path_abelian);
  FLOAT aver[2];

  cout.precision(10);
  su2 A;
  unsigned int start_time = clock();

  conf_abelian_smeared.array =
      link_double.smearing_APE(conf_abelian.array, alpha_APE);

  for (int i = 0; i < 10; i++) {
    cout << "smearing_APE_abelian test " << conf_abelian_smeared.array[i]
         << endl;
  }

  cout << "wilson_abelian_smeared aver "
       << wilson(conf_abelian_smeared.array, 1, 1) << endl;

  vector<vector<abelian>> smearing_first_abelian(6, vector<abelian>(1));
  vector<vector<abelian>> smearing_second_abelian(4, vector<abelian>(1));
  smearing_first_abelian =
      link_double.smearing_first_full(conf_abelian.array, alpha3);
  smearing_second_abelian = link_double.smearing_second_full(
      conf_abelian.array, smearing_first_abelian, alpha2);
  conf_abelian_smeared.array = link_double.smearing_HYP(
      conf_abelian.array, smearing_second_abelian, alpha1);

  smearing_first = link.smearing_first_full(conf.array, alpha3);
  smearing_second =
      link.smearing_second_full(conf.array, smearing_first, alpha2);
  smeared.array = link.smearing_HYP(conf.array, smearing_second, alpha1);

  for (int i = 0; i < 4; i++) {
    A = smeared.array[i];
    cout << "smearing_HYP test " << A.a0 << " " << A.a1 << " " << A.a2 << " "
         << A.a3 << endl;
  }
  cout << "right:" << endl;
  cout << "-0.110943 0.629543 -0.583767 0.500582" << endl;
  cout << "0.318657 -0.218438 -0.415944 0.823245" << endl;
  cout << "0.174894 0.667862 0.674797 0.260808" << endl;
  cout << "-0.704066 -0.13241 0.453476 -0.530206" << endl;

  smeared.array = link.smearing_HYP_refresh(conf, alpha1, alpha2, alpha3);

  for (int i = 0; i < 4; i++) {
    A = smeared.array[i];
    cout << "smearing_HYP_refresh test " << A.a0 << " " << A.a1 << " " << A.a2
         << " " << A.a3 << endl;
  }
  cout << "right:" << endl;
  cout << "-0.110943 0.629543 -0.583767 0.500582" << endl;
  cout << "0.318657 -0.218438 -0.415944 0.823245" << endl;
  cout << "0.174894 0.667862 0.674797 0.260808" << endl;
  cout << "-0.704066 -0.13241 0.453476 -0.530206" << endl;

  smeared.array = link.smearing_APE(conf.array, alpha_APE);

  for (int i = 0; i < 4; i++) {
    A = smeared.array[i];
    cout << "smearing_APE test " << A.a0 << " " << A.a1 << " " << A.a2 << " "
         << A.a3 << endl;
  }
  cout << "right:" << endl;
  cout << "-0.0833883 0.641689 -0.604086 0.465147" << endl;
  cout << "0.336007 -0.177589 -0.426239 0.820903" << endl;
  cout << "0.125318 0.626699 0.730398 0.240963" << endl;
  cout << "-0.632818 0.13761 0.745239 0.158818" << endl;

  smeared.array = link.smearing_stout(conf, stout_rho);

  for (int i = 0; i < 4; i++) {
    A = smeared.array[i];
    cout << "smearing_stout test " << A.a0 << " " << A.a1 << " " << A.a2 << " "
         << A.a3 << endl;
  }
  cout << "right:" << endl;
  cout << "-0.0185375 0.684286 -0.660024 0.30948" << endl;
  cout << "0.363087 -0.0991993 -0.575234 0.726246" << endl;
  cout << "0.172742 0.633665 0.73173 0.182208" << endl;
  cout << "-0.816237 -0.0404862 0.571983 -0.0703786" << endl;

  unsigned int end_time = clock();
  unsigned int search_time = end_time - start_time;
  cout << "working time: " << search_time * 1. / CLOCKS_PER_SEC << endl;
}
