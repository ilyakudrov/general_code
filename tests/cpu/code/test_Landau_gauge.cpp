#include "../../../lib/cpu/include/Landau_gauge.h"
#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/flux_tube.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/smearing.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <map>
#include <numeric>
#include <stdio.h>
#include <tuple>
#include <unordered_map>
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

  x_size = 40;
  y_size = 40;
  z_size = 40;
  t_size = 40;

  cout.precision(17);

  data<su2> conf;
  // data<su3_full> conf;
  // data<abelian> conf;
  string conf_path =
      "../../confs/MA_gauge/qc2dstag/40^4/mu0.05/s0/conf_abelian_0201";
  // string conf_path = "../../confs/su2_suzuki/monopoless/CON_OFF_MAG_001.LAT";
  // string conf_path = "../../confs/su2_suzuki/monopole/CON_MON_MAG_001.LAT";
  // conf.read_double(conf_path, 8);
  conf.read_double(conf_path, 0);
  // conf.read_float(conf_path, 8);

  // plakets and polyakov loop
  //   start_time = clock();
  cout << "qc2dstag plaket " << plaket(conf.array) << endl;
  //   cout << "qc2dstag plaket_time " << plaket_time(conf.array) << endl;
  //   cout << "qc2dstag plaket_space " << plaket_space(conf.array) << endl;
  //   cout << "qc2dstag polyakov " << polyakov(conf.array) << endl;
  //   end_time = clock();
  //   search_time = end_time - start_time;
  //   cout << "plaket and staff time: " << search_time * 1. / CLOCKS_PER_SEC
  //        << endl;

  cout << "initial Landau functional: " << functional_Landau(conf.array)
       << endl;

  vector<abelian> angles_conf = convert_su2_abelian(conf.array);
  vector<double> angles_gauge = generate_angles_gauge_random();

  cout << "initial Landau functional abelian: "
       << functional_Landau(angles_conf) << endl;

  double temperature_start = 2.5;
  double temperature_end = 0.01;
  double temperature_step = 0.01;

  int thermalization_steps = 20;

  int n = 5;

  start_time = clock();

  for (int i = 0; i < n; i++) {
    heat_bath_step(angles_gauge, angles_conf, 2);
  }

  end_time = clock();
  search_time = end_time - start_time;
  std::cout << "heat bath step time: "
            << (double)search_time / n / CLOCKS_PER_SEC << std::endl;

  thermalize(angles_gauge, angles_conf, temperature_start,
             thermalization_steps);

  vector<SA_data> SA_data =
      simulated_annealing_test(angles_gauge, angles_conf, temperature_start,
                               temperature_end, temperature_step);

  for (int i = 0; i < 10; i++) {
    relaxation_step(angles_gauge, angles_conf);
  }

  std::ofstream stream_SA;
  stream_SA.precision(17);
  stream_SA.open("SA_test");
  stream_SA << "#temperature,functional" << std::endl;

  for (auto i : SA_data) {
    stream_SA << i.temperature << "," << i.functional << std::endl;
  }

  stream_SA.close();

  cout << "final Landau functional: "
       << functional_Landau(angles_conf, angles_gauge) << endl;
}
