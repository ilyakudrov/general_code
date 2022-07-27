#include "../../../lib/cpu/include/Landau_U1.h"
#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/decomposition.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/monopoles.h"

#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
#include <tuple>
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

  double T_init = 4;
  double T_final = 0.2;
  double T_step = 0.2;
  int OR_steps = 4;
  int thermalization_steps = 30;
  int tolerance_digits = 7;
  double tolerance_maximal = 1e-5;
  double tolerance_average = 1e-7;

  string conf_path =
      "../../confs/qc2dstag/40^4/mag/mu0.05/s0/conf_abelian_0201";
  //   string conf_path =
  //   "../../confs/MA_gauge/su2_suzuki/conf_gaugefixed/24^4/"
  //                      "beta2.4/conf_gaugefixed_0001";
  // string conf_path =
  //     "../../confs/qc2dstag/40^4/Landau_U1/mu0.05/s0/conf_Landau_U1";

  conf.read_double(conf_path, 0);
  // conf.read_double(conf_path, 4);

  std::vector<abelian> conf_abelian = convert_to_abelian(conf.array);

  cout << "initial plaket " << plaket(conf.array) << endl;
  cout << "initial plaket from su2 conf " << Landau_functional(conf.array)
       << endl;

  cout << "initial functional " << Landau_functional_abelian(conf_abelian)
       << endl;

  /*gauge_abelian = generate_gauge_abelian_uniform();

    start_time = clock();

    make_simulated_annealing_test1(gauge_abelian, conf_abelian, T_init,
    T_final,
                                   T_step, OR_steps, thermalization_steps);

    end_time = clock();
    search_time = end_time - start_time;
    cout << "make_simulated_annealing_test1 time: "
         << search_time * 1. / CLOCKS_PER_SEC << endl;

  make_maximization_final(gauge_abelian, conf_abelian, OR_steps,
                          tolerance_maximal, tolerance_average);

  cout << "final functional "
       << Landau_functional_gauge_abelian(gauge_abelian, conf_abelian) << endl;

  gauge_tranformation_abelian(gauge_abelian, conf_abelian);

  gauge_tranformation(gauge_abelian, conf.array);

  cout << "final plaket " << plaket(conf.array) << endl;

  cout << "final functional abelian " << Landau_functional_abelian(conf_abelian)
       << endl;

  cout << "final functional conf " << Landau_functional(conf.array) << endl;

  gauge_abelian = generate_gauge_abelian_uniform();

  start_time = clock();

  heat_bath_update_test1(gauge_abelian, conf_abelian, T_init);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "heat_bath_update_test1 time: " << search_time * 1. / CLOCKS_PER_SEC
       << endl;

  cout << "heat_bath_update_test1 functional "
       << Landau_functional_gauge_abelian(gauge_abelian, conf_abelian) <<
  endl;

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
       << Landau_functional_gauge(gauge_angles, conf_angles) << endl;

  gauge_complex = generate_gauge_complex_uniform();
  conf_complex = convert_to_complex(conf.array);

  start_time = clock();

  heat_bath_update(gauge_complex, conf_complex, T_init);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "heat_bath_update time: " << search_time * 1. / CLOCKS_PER_SEC
       << endl;

  cout << "heat_bath_update functional "
       << Landau_functional_gauge_complex(gauge_complex, conf_complex) <<
  endl;*/

  conf_angles = convert_to_angles(conf.array);

  std::vector<std::vector<int>> monopole_plaket =
      calculate_monopole_plaket_singular(conf_angles);

  std::vector<std::vector<int>> monopole_difference(4, std::vector<int>());
  std::vector<std::vector<int>> monopole_coordinate(4, std::vector<int>());

  monopole_plaket_difference_nonzero(monopole_plaket, monopole_difference,
                                     monopole_coordinate);

  for (int mu = 0; mu < 4; mu++) {
    cout << "monopole difference number " << mu << " "
         << monopole_difference[mu].size() << endl;
  }

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

  start_time = clock();

  cout << "simulated_annealing functional "
       << Landau_functional_gauge_complex(gauge_complex, conf_complex) << endl;

  make_maximization_final(gauge_complex, conf_complex, OR_steps,
                          tolerance_maximal, tolerance_average);

  end_time = clock();
  search_time = end_time - start_time;
  cout << "make_maximization_final time: " << search_time * 1. / CLOCKS_PER_SEC
       << endl;

  normalize_complex(gauge_complex);

  double unitarity_test_aver = 0;
  double unitarity_test_max = 0;
  double tmp = 0;

  for (int i = 0; i < gauge_complex.size(); i++) {
    tmp = 1 - sqrt(gauge_complex[i].imag * gauge_complex[i].imag +
                   gauge_complex[i].real * gauge_complex[i].real);
    if (tmp > unitarity_test_max) {
      unitarity_test_max = tmp;
    }
    unitarity_test_aver += tmp;
  }
  unitarity_test_aver = unitarity_test_aver / gauge_complex.size();
  cout << "unitarity deviation aver = " << unitarity_test_aver
       << ", max = " << unitarity_test_max << endl;

  cout << "final functional "
       << Landau_functional_gauge_complex(gauge_complex, conf_complex) << endl;

  apply_gauge_Landau_complex(gauge_complex, conf_complex);

  cout << "final after applying gauge "
       << Landau_functional_complex(conf_complex) << endl;

  conf_angles = convert_complex_to_angles(conf_complex);

  monopole_plaket = calculate_monopole_plaket_singular(conf_angles);

  for (int mu = 0; mu < 4; mu++) {

    monopole_difference[mu].clear();
    monopole_difference[mu].shrink_to_fit();

    monopole_coordinate[mu].clear();
    monopole_coordinate[mu].shrink_to_fit();
  }

  monopole_plaket_difference_nonzero(monopole_plaket, monopole_difference,
                                     monopole_coordinate);

  for (int mu = 0; mu < 4; mu++) {
    cout << "monopole difference number " << mu << " "
         << monopole_difference[mu].size() << endl;
  }

  // thermalization digram test

  /*string path_out = "thermalization_Landau_U1";

  ofstream stream;
  stream.open(path_out);

  map<tuple<int, double>, double> functional_Landau_U1;

  gauge_complex = generate_gauge_complex_uniform();
  conf_complex = convert_to_complex(conf.array);

  start_time = clock();

  for (int i = 0; i < thermalization_steps; i++) {
    heat_bath_update(gauge_complex, conf_complex, T_init);
  }

  double T = T_init;
  int step = 0;

  while (T > T_final) {

    for (int i = 0; i < 10; i++) {

      heat_bath_update(gauge_complex, conf_complex, T);

      for (int i = 0; i < OR_steps; i++) {
        overrelaxation_update(gauge_complex, conf_complex);
      }
    }

    for (int i = 0; i < 30; i++) {

      heat_bath_update(gauge_complex, conf_complex, T);

      for (int i = 0; i < OR_steps; i++) {
        overrelaxation_update(gauge_complex, conf_complex);
      }

      functional_Landau_U1[tuple<int, double>(i, T)] =
          Landau_functional_gauge_complex(gauge_complex, conf_complex);

      //  step++;
    }

    if (T > 2)
      T -= T_step;
    else
      T -= T_step / 2;
  }

  end_time = clock();
  search_time = end_time - start_time;
  cout << "thermalization time: " << search_time * 1. / CLOCKS_PER_SEC << endl;

  cout << "functional after thermalization "
       << Landau_functional_gauge_complex(gauge_complex, conf_complex) << endl;

  stream << "step,T,functional" << endl;

  for (auto it = functional_Landau_U1.begin(); it != functional_Landau_U1.end();
       it++) {
    stream << get<0>(it->first) << "," << get<1>(it->first) << "," << it->second
           << endl;
  }

  stream.close();

  //   start_time = clock();

  //   make_maximization_final(gauge_complex, conf_complex, OR_steps,
  //                           tolerance_maximal, tolerance_average);

  //   end_time = clock();
  //   search_time = end_time - start_time;
  //   cout << "make_maximization_final time: " << search_time * 1. /
  //   CLOCKS_PER_SEC
  //        << endl;

  //   cout << "final functional "
  //        << Landau_functional_gauge_complex(gauge_complex, conf_complex) <<
  //        endl;*/
}