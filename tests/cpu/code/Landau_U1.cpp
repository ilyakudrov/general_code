#include "../../../lib/cpu/include/Landau_U1.h"
#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/matrix.h"

#include <ctime>
#include <fstream>
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

  x_size = 40;
  y_size = 40;
  z_size = 40;
  t_size = 40;

  std::cout.precision(17);

  data<su2> conf;

  vector<su2> gauge;
  data<su2> conf1;

  std::vector<abelian> gauge_abelian;
  std::vector<double> gauge_angles;
  std::vector<double> conf_angles;
  std::vector<complex_t> gauge_complex;
  std::vector<complex_t> conf_complex;

  double T_init = 2.5;
  double T_final = 0.01;
  double T_step = 0.01;
  int OR_steps = 4;
  int thermalization_steps = 5;
  int tolerance_digits = 7;
  double tolerance_maximal = 1e-5;
  double tolerance_average = 1e-7;

  string conf_path =
      "../../confs/qc2dstag/40^4/mag/mu0.05/s0/conf_abelian_0201";
  // string conf_path =
  //     "../../confs/qc2dstag/40^4/Landau_U1/mu0.05/s0/conf_Landau_U1";

  conf.read_double(conf_path, 0);
  // conf.read_double(conf_path, 4);

  // for (int i = 0; i < 4; i++) {
  //   cout << conf.array[i] << endl;
  // }

  std::vector<abelian> conf_abelian = convert_to_abelian(conf.array);

  // simulated annealing and maximization

  cout << "initial plaket " << plaket(conf.array) << endl;
  cout << "initial plaket from su2 conf " << Landau_functional(conf.array)
       << endl;

  //   double temperature = 2.5;

  cout << "initial functional " << Landau_functional_abelian(conf_abelian)
       << endl;

  //   gauge_abelian = generate_gauge_abelian_uniform();

  //   start_time = clock();

  //   make_simulated_annealing_test1(gauge_abelian, conf_abelian, T_init,
  //   T_final,
  //                                  T_step, OR_steps, thermalization_steps);

  //   end_time = clock();
  //   search_time = end_time - start_time;
  //   cout << "make_simulated_annealing_test1 time: "
  //        << search_time * 1. / CLOCKS_PER_SEC << endl;

  //   make_maximization_final(gauge_abelian, conf_abelian, OR_steps,
  //                           tolerance_maximal, tolerance_average);

  //   cout << "final functional "
  //        << Landau_functional_gauge_abelian(gauge_abelian, conf_abelian) <<
  //        endl;

  //   gauge_tranformation_abelian(gauge_abelian, conf_abelian);

  //   gauge_tranformation(gauge_abelian, conf.array);

  //   cout << "final plaket " << plaket(conf.array) << endl;

  //   cout << "final functional abelian " <<
  //   Landau_functional_abelian(conf_abelian)
  //        << endl;

  //   cout << "final functional abelian " << Landau_functional(conf.array) <<
  //   endl;

  /*gauge_abelian = generate_gauge_abelian_uniform();

  start_time = clock();

  heat_bath_update_test1(gauge_abelian, conf_abelian, T_init);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "heat_bath_update_test1 time: " << search_time * 1. / CLOCKS_PER_SEC
       << endl;

  cout << "heat_bath_update_test1 functional "
       << Landau_functional_gauge_abelian(gauge_abelian, conf_abelian) << endl;

  gauge_angles = generate_gauge_angles_uniform();

  start_time = clock();

  heat_bath_update_test2(gauge_angles, conf_abelian, T_init);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "heat_bath_update_test2 time: " << search_time * 1. / CLOCKS_PER_SEC
       << endl;

  cout << "heat_bath_update_test2 functional "
       << Landau_functional_gauge_angles(gauge_angles, conf_abelian) << endl;

  gauge_angles = generate_gauge_angles_uniform();
  conf_angles = convert_to_angles(conf.array);

  start_time = clock();

  heat_bath_update_test3(gauge_angles, conf_angles, T_init);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "heat_bath_update_test3 time: " << search_time * 1. / CLOCKS_PER_SEC
       << endl;

  cout << "heat_bath_update_test3 functional "
       << Landau_functional_gauge(gauge_angles, conf_angles) << endl;*/

  //   gauge_complex = generate_gauge_complex_uniform();
  //   conf_complex = convert_to_complex(conf.array);

  //   start_time = clock();

  //   heat_bath_update(gauge_complex, conf_complex, T_init);

  //   end_time = clock();
  //   search_time = end_time - start_time;
  //   cout << "heat_bath_update time: " << search_time * 1. / CLOCKS_PER_SEC
  //        << endl;

  //   cout << "heat_bath_update functional "
  //        << Landau_functional_gauge_complex(gauge_complex, conf_complex) <<
  //        endl;

  gauge_complex = generate_gauge_complex_uniform();
  conf_complex = convert_to_complex(conf.array);

  //   start_time = clock();

  //   make_simulated_annealing(gauge_complex, conf_complex, T_init, T_final,
  //   T_step,
  //                            OR_steps, thermalization_steps);

  //   end_time = clock();
  //   search_time = end_time - start_time;
  //   cout << "simulated_annealing time: " << search_time * 1. / CLOCKS_PER_SEC
  //        << endl;

  make_maximization_final(gauge_complex, conf_complex, OR_steps,
                          tolerance_maximal, tolerance_average);

  cout << "final functional "
       << Landau_functional_gauge_complex(gauge_complex, conf_complex) << endl;
}