#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/mag.h"
#include "../../../lib/cpu/include/matrix.h"

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

  x_size = 40;
  y_size = 40;
  z_size = 40;
  t_size = 40;

  std::cout.precision(17);

  data<su2> conf;

  string conf_path = "../../confs/qc2dstag/mu0.05/s0/CONF0201";

  conf.read_double_qc2dstag(conf_path);

  cout << "initial plaket " << plaket(conf.array) << endl;

  std::vector<spin> spins;

  double temperature = 2.5;

  cout << "initial functional " << MAG_functional_su2(conf.array) << endl;

  spins = generate_spins_uniform();

  start_time = clock();

  heat_bath_update_tets1(spins, conf.array, temperature);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "heat_bath_update_tets1 time: " << search_time * 1. / CLOCKS_PER_SEC
       << endl;

  cout << "final functional " << MAG_functional_su2_spin(conf.array, spins)
       << endl;

  vector<su2> gauge;
  data<su2> conf1;

  gauge = make_gauge(spins);
  conf1.array = gauge_tranformation(conf.array, gauge);

  cout << "final plaket " << plaket(conf1.array) << endl;

  spins = generate_spins_uniform();

  vector<int> indices_qube = make_indices_qube(2);

  start_time = clock();

  heat_bath_update_tets2(spins, conf.array, indices_qube, temperature);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "heat_bath_update_tets2 time: " << search_time * 1. / CLOCKS_PER_SEC
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
  cout << "heat_bath_update_tets3 time: " << search_time * 1. / CLOCKS_PER_SEC
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

  heat_bath_update_tets4(spins, conf.array, indices, coordinates, temperature);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "heat_bath_update_tets4 time: " << search_time * 1. / CLOCKS_PER_SEC
       << endl;

  cout << "final functional " << MAG_functional_su2_spin(conf.array, spins)
       << endl;

  gauge = make_gauge(spins);
  conf1.array = gauge_tranformation(conf.array, gauge);

  cout << "final plaket " << plaket(conf1.array) << endl;
}