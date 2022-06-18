#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/mag.h"
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

  std::vector<spin> spins;

  double T_init = 2.5;
  double T_final = 0.1;
  double T_step = 0.1;
  int OR_steps = 6;
  int thermalization_steps = 20;
  int tolerance_digits = 7;
  double tolerance_maximal = 1e-9;
  double tolerance_average = 1e-11;

  string conf_path = "../../confs/qc2dstag/mu0.05/s0/CONF0201";

  conf.read_double_qc2dstag(conf_path);

  // performance test

  /*cout << "initial plaket " << plaket(conf.array) << endl;

  double temperature = 2.5;

  cout << "initial functional " << MAG_functional_su2(conf.array) << endl;

  spins = generate_spins_uniform();

  start_time = clock();

  heat_bath_update(spins, conf.array, temperature);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "heat_bath_update_tets1 time: " << search_time * 1. /
  CLOCKS_PER_SEC
       << endl;

  cout << "final functional " << MAG_functional_su2_spin(conf.array, spins)
       << endl;

  gauge = make_gauge(spins);
  conf1.array = gauge_tranformation(conf.array, gauge);

  cout << "final plaket " << plaket(conf1.array) << endl;

  spins = generate_spins_uniform();

  vector<int> indices_qube = make_indices_qube(2);

  start_time = clock();

  heat_bath_update_tets2(spins, conf.array, indices_qube, temperature);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "heat_bath_update_tets2 time: " << search_time * 1. /
  CLOCKS_PER_SEC
       << endl;

  cout << "final functional " << MAG_functional_su2_spin(conf.array, spins)
       << endl;

  gauge = make_gauge(spins);
  conf1.array = gauge_tranformation(conf.array, gauge);

  cout << "final plaket " << plaket(conf1.array) << endl;

  spins = generate_spins_uniform();

  start_time = clock();

  heat_bath_update_tets3(spins, conf.array, temperature);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "heat_bath_update_tets3 time: " << search_time * 1. /
  CLOCKS_PER_SEC
       << endl;

  cout << "final functional " << MAG_functional_su2_spin(conf.array, spins)
       << endl;

  gauge = make_gauge(spins);
  conf1.array = gauge_tranformation(conf.array, gauge);

  cout << "final plaket " << plaket(conf1.array) << endl;

  spins = generate_spins_uniform();

  std::vector<int> indices;
  std::vector<char> coordinates;
  indices.reserve(x_size * y_size * z_size * t_size);
  coordinates.reserve(x_size * y_size * z_size * t_size * 4);

  make_indices_qube1(indices, coordinates, 2);

  start_time = clock();

  heat_bath_update_tets4(spins, conf.array, indices, coordinates,
  temperature);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "heat_bath_update_tets4 time: " << search_time * 1. /
  CLOCKS_PER_SEC
       << endl;

  cout << "final functional " << MAG_functional_su2_spin(conf.array, spins)
       << endl;

  gauge = make_gauge(spins);
  conf1.array = gauge_tranformation(conf.array, gauge);

  cout << "final plaket " << plaket(conf1.array) << endl;*/

  // thermalization test

  /*string path_out = "thermalization";
  ofstream stream;
  stream.open(path_out);

  map<double, double> functional;

  spins = generate_spins_uniform();

  for (int i = 0; i < thermalization_steps; i++) {
    heat_bath_update(spins, conf.array, T_init);
  }

  double T = T_init;

  while (T > T_final) {

    heat_bath_update(spins, conf.array, T);

    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(spins, conf.array);
    }

    functional[T] = MAG_functional_su2_spin(conf.array, spins);

    T -= T_step;
  }

  stream << "T,functional" << endl;

  for (auto it = functional.begin(); it != functional.end(); it++) {
    stream << it->first << "," << it->second << endl;
  }

  stream.close();

  cout << "final functional " << MAG_functional_su2_spin(conf.array, spins)
       << endl;

  gauge = make_gauge(spins);
  conf1.array = gauge_tranformation(conf.array, gauge);

  cout << "final functional from su2 conf " << MAG_functional_su2(conf1.array)
       << endl;
  cout << "final plaket " << plaket(conf1.array) << endl;

  conf1.array = conf.array;

  gauge_tranformation_spins(conf1.array, spins);

  cout << "final functional from su2 conf " << MAG_functional_su2(conf1.array)
       << endl;
  cout << "final plaket " << plaket(conf1.array) << endl;*/

  // simulated annealing and maximization

  spins = generate_spins_uniform();

  start_time = clock();

  heat_bath_update(spins, conf.array, 1.);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "heat_bath_update time: " << search_time * 1. / CLOCKS_PER_SEC
       << endl;

  start_time = clock();

  overrelaxation_update(spins, conf.array);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "overrelaxation_update time: " << search_time * 1. / CLOCKS_PER_SEC
       << endl;

  spins = generate_spins_uniform();

  start_time = clock();

  make_simulated_annealing(conf.array, spins, T_init, T_final, T_step, OR_steps,
                           thermalization_steps);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "simulated_annealing time: " << search_time * 1. / CLOCKS_PER_SEC
       << endl;

  make_maximization_approximate(conf.array, spins, OR_steps, tolerance_digits);

  make_maximization_final(conf.array, spins, OR_steps, tolerance_maximal,
                          tolerance_average);
}