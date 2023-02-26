#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/flux_tube.h"
#include "../../../lib/cpu/include/matrix.h"

#include <cstring>
#include <ctime>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <vector>

#ifndef MATRIX_PLAKET
#define MATRIX_PLAKET su2
#endif
#ifndef MATRIX_WILSON
#define MATRIX_WILSON su2
#endif

using namespace std;

int x_size;
int y_size;
int z_size;
int t_size;

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double search_time;

  std::string conf_format_plaket;
  std::string conf_format_wilson;
  std::string conf_path_plaket;
  std::string conf_path_wilson;
  std::string output_path_electric_long;
  std::string output_path_magnetic_long;
  std::string output_path_electric_trans;
  std::string output_path_magnetic_trans;
  int L_spat, L_time;
  int x_trans = 0;
  int bytes_skip_plaket = 0;
  int bytes_skip_wilson = 0;
  int T_min, T_max;
  int R_min, R_max;
  bool convert_plaket = 0;
  bool convert_wilson = 0;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-conf_format_plaket") == 0) {
      conf_format_plaket = argv[++i];
    } else if (string(argv[i]) == "-bytes_skip_plaket") {
      bytes_skip_plaket = stoi(string(argv[++i]));
    } else if (std::string(argv[i]) == "-conf_format_wilson") {
      conf_format_wilson = argv[++i];
    } else if (string(argv[i]) == "-bytes_skip_wilson") {
      bytes_skip_wilson = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-convert_wilson") {
      istringstream(string(argv[++i])) >> convert_wilson;
    } else if (string(argv[i]) == "-convert_plaket") {
      istringstream(string(argv[++i])) >> convert_plaket;
    } else if (std::string(argv[i]) == "-conf_path_plaket") {
      conf_path_plaket = argv[++i];
    } else if (std::string(argv[i]) == "-conf_path_wilson") {
      conf_path_wilson = argv[++i];
    } else if (std::string(argv[i]) == "-output_path_electric_long") {
      output_path_electric_long = argv[++i];
    } else if (std::string(argv[i]) == "-output_path_magnetic_long") {
      output_path_magnetic_long = argv[++i];
    } else if (std::string(argv[i]) == "-output_path_electric_trans") {
      output_path_electric_trans = argv[++i];
    } else if (std::string(argv[i]) == "-output_path_magnetic_trans") {
      output_path_magnetic_trans = argv[++i];
    } else if (std::string(argv[i]) == "-L_spat") {
      L_spat = stoi(std::string(argv[++i]));
    } else if (std::string(argv[i]) == "-L_time") {
      L_time = stoi(std::string(argv[++i]));
    } else if (std::string(argv[i]) == "-x_trans") {
      x_trans = stoi(std::string(argv[++i]));
    } else if (std::string(argv[i]) == "-T_min") {
      T_min = stoi(std::string(argv[++i]));
    } else if (std::string(argv[i]) == "-T_max") {
      T_max = stoi(std::string(argv[++i]));
    } else if (std::string(argv[i]) == "-R_min") {
      R_min = stoi(std::string(argv[++i]));
    } else if (std::string(argv[i]) == "-R_max") {
      R_max = stoi(std::string(argv[++i]));
    }
  }

  std::cout << "conf_format_plaket " << conf_format_plaket << std::endl;
  cout << "bytes_skip_plaket " << bytes_skip_plaket << endl;
  cout << "convert_plaket " << convert_plaket << endl;
  std::cout << "conf_format_wilson " << conf_format_wilson << std::endl;
  cout << "bytes_skip_wilson " << bytes_skip_wilson << endl;
  cout << "convert_wilson " << convert_wilson << endl;
  std::cout << "conf_path_plaket " << conf_path_plaket << std::endl;
  std::cout << "conf_path_wilson " << conf_path_wilson << std::endl;
  std::cout << "output_path_electric_long " << output_path_electric_long
            << std::endl;
  std::cout << "output_path_magnetic_long " << output_path_magnetic_long
            << std::endl;
  std::cout << "output_path_electric_trans " << output_path_electric_trans
            << std::endl;
  std::cout << "output_path_magnetic_trans " << output_path_magnetic_trans
            << std::endl;
  std::cout << "L_spat " << L_spat << std::endl;
  std::cout << "L_time " << L_time << std::endl;
  std::cout << "R_min " << R_min << std::endl;
  std::cout << "R_max " << R_max << std::endl;
  std::cout << "T_min " << T_min << std::endl;
  std::cout << "T_max " << T_max << std::endl;
  std::cout << "x_trans " << x_trans << std::endl;

  x_size = L_spat;
  y_size = L_spat;
  z_size = L_spat;
  t_size = L_time;

  data<MATRIX_PLAKET> conf_plaket;
  data<MATRIX_WILSON> conf_wilson;

  get_data(conf_plaket, conf_path_plaket, conf_format_plaket, bytes_skip_plaket,
           convert_plaket);
  get_data(conf_wilson, conf_path_wilson, conf_format_wilson, bytes_skip_wilson,
           convert_wilson);

  double plaket_time_average = plaket_time(conf_plaket.array);
  double plaket_space_average = plaket_space(conf_plaket.array);

  std::cout << "plaket_time " << plaket_time_average << " smeared_plaket_time "
            << plaket_time(conf_wilson.array) << std::endl;
  std::cout << "plaket_space " << plaket_space_average
            << " smeared_plaket_space " << plaket_space(conf_wilson.array)
            << std::endl;

  std::ofstream stream_electric_long;
  std::ofstream stream_magnetic_long;
  std::ofstream stream_electric_trans;
  std::ofstream stream_magnetic_trans;

  stream_electric_long.precision(17);
  stream_magnetic_long.precision(17);
  stream_electric_trans.precision(17);
  stream_magnetic_trans.precision(17);

  stream_electric_long.open(output_path_electric_long);
  stream_magnetic_long.open(output_path_magnetic_long);
  stream_electric_trans.open(output_path_electric_trans);
  stream_magnetic_trans.open(output_path_magnetic_trans);

  stream_electric_long << "T,R,d,correlator,wilson_loop,plaket" << endl;
  stream_magnetic_long << "T,R,d,correlator,wilson_loop,plaket" << endl;
  stream_electric_trans << "T,R,d,correlator,wilson_loop,plaket" << endl;
  stream_magnetic_trans << "T,R,d,correlator,wilson_loop,plaket" << endl;

  vector<vector<MATRIX_PLAKET>> separated_plaket =
      separate_wilson(conf_plaket.array);
  conf_plaket.array.erase(conf_plaket.array.begin(), conf_plaket.array.end());
  vector<vector<MATRIX_WILSON>> separated_wilson =
      separate_wilson(conf_wilson.array);
  conf_wilson.array.erase(conf_wilson.array.begin(), conf_wilson.array.end());

  map<tuple<int, int>, double> wilson_loops =
      wilson_parallel(separated_wilson, R_min, R_max, T_min, T_max);

  map<tuple<int, int, int>, double> flux_tube;

  vector<double> plaket_time_tr = plaket_aver_tr_time(separated_plaket);
  double plaket_time_aver = 0;
  for (int i = 0; i < plaket_time_tr.size(); i++) {
    plaket_time_aver += plaket_time_tr[i];
  }
  plaket_time_aver = plaket_time_aver / plaket_time_tr.size();
  vector<double> plaket_space_tr = plaket_aver_tr_space(separated_plaket);
  double plaket_space_aver = 0;
  for (int i = 0; i < plaket_space_tr.size(); i++) {
    plaket_space_aver += plaket_space_tr[i];
  }
  plaket_space_aver = plaket_space_aver / plaket_space_tr.size();

  start_time = omp_get_wtime();

  flux_tube =
      wilson_plaket_correlator(plaket_time_tr, separated_wilson, T_min, T_max,
                               R_min, R_max, 10, 0, "longitudinal");

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "flux tube longitudinal electric new time: " << search_time << endl;

  for (auto it = flux_tube.begin(); it != flux_tube.end(); it++) {
    stream_electric_long
        << get<0>(it->first) << "," << get<1>(it->first) << ","
        << get<2>(it->first) << "," << it->second << ","
        << wilson_loops[tuple<int, int>(get<0>(it->first), get<1>(it->first))]
        << "," << plaket_time_aver << endl;
  }

  start_time = omp_get_wtime();

  flux_tube =
      wilson_plaket_correlator(plaket_space_tr, separated_wilson, T_min, T_max,
                               R_min, R_max, 10, 0, "longitudinal");

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "flux tube longitudinal magnetic new time: " << search_time << endl;

  for (auto it = flux_tube.begin(); it != flux_tube.end(); it++) {
    stream_magnetic_long
        << get<0>(it->first) << "," << get<1>(it->first) << ","
        << get<2>(it->first) << "," << it->second << ","
        << wilson_loops[tuple<int, int>(get<0>(it->first), get<1>(it->first))]
        << "," << plaket_space_aver << endl;
  }

  start_time = omp_get_wtime();

  flux_tube =
      wilson_plaket_correlator(plaket_time_tr, separated_wilson, T_min, T_max,
                               R_min, R_max, x_size / 2, 0, "transversal");

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "flux tube transversal electric new time: " << search_time << endl;

  for (auto it = flux_tube.begin(); it != flux_tube.end(); it++) {
    stream_electric_long
        << get<0>(it->first) << "," << get<1>(it->first) << ","
        << get<2>(it->first) << "," << it->second << ","
        << wilson_loops[tuple<int, int>(get<0>(it->first), get<1>(it->first))]
        << "," << plaket_time_aver << endl;
  }

  start_time = omp_get_wtime();

  flux_tube =
      wilson_plaket_correlator(plaket_space_tr, separated_wilson, T_min, T_max,
                               R_min, R_max, x_size / 2, 0, "transversal");

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "flux tube transversal magnetic new time: " << search_time << endl;

  for (auto it = flux_tube.begin(); it != flux_tube.end(); it++) {
    stream_magnetic_long
        << get<0>(it->first) << "," << get<1>(it->first) << ","
        << get<2>(it->first) << "," << it->second << ","
        << wilson_loops[tuple<int, int>(get<0>(it->first), get<1>(it->first))]
        << "," << plaket_space_aver << endl;
  }
}
