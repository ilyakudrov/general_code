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
  double observables_time;

  std::string conf_format_plaket;
  std::string conf_format_wilson;
  std::string conf_path_plaket;
  std::string conf_path_wilson;
  std::string output_path_electric;
  std::string output_path_magnetic;
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
    } else if (std::string(argv[i]) == "-output_path_electric") {
      output_path_electric = argv[++i];
    } else if (std::string(argv[i]) == "-output_path_magnetic") {
      output_path_magnetic = argv[++i];
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
  std::cout << "output_path_electric " << output_path_electric << std::endl;
  std::cout << "output_path_magnetic " << output_path_magnetic << std::endl;
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

  vector<double> wilson_loop_trace;
  vector<double> plaket_time_trace =
      calculate_plaket_time_tr(conf_plaket.array);
  vector<double> plaket_space_trace =
      calculate_plaket_space_tr(conf_plaket.array);
  double wilson_loop_average;

  std::cout << "plaket_time " << plaket_time_average << " smeared_plaket_time "
            << plaket_time(conf_wilson.array) << std::endl;
  std::cout << "plaket_space " << plaket_space_average
            << " smeared_plaket_space " << plaket_space(conf_wilson.array)
            << std::endl;

  std::ofstream stream_electric;
  std::ofstream stream_magnetic;

  stream_electric.precision(17);
  stream_magnetic.precision(17);

  stream_electric.open(output_path_electric);
  stream_magnetic.open(output_path_magnetic);

  cout << output_path_electric << endl;
  cout << output_path_magnetic << endl;

  stream_electric << "T,R,d,wilson-plaket-correlator,wilson-loop,plaket"
                  << std::endl;

  stream_magnetic << "T,R,d,wilson-plaket-correlator,wilson-loop,plaket"
                  << std::endl;

  std::map<int, double> correlator_electric;
  std::map<int, double> correlator_magnetic;

  start_time = omp_get_wtime();

  for (int T = T_min; T <= T_max; T += 2) {
    for (int R = R_min; R <= R_max; R += 2) {

      int d_min = -5;
      int d_max = R + 5;

      wilson_loop_trace = calculate_wilson_loop_tr(conf_wilson.array, R, T);

      wilson_loop_average = accumulate(wilson_loop_trace.cbegin(),
                                       wilson_loop_trace.cend(), 0.0) /
                            wilson_loop_trace.size();

      correlator_electric = wilson_plaket_correlator_electric(
          wilson_loop_trace, plaket_time_trace, R, T, x_trans, d_min, d_max);
      correlator_magnetic = wilson_plaket_correlator_magnetic(
          wilson_loop_trace, plaket_space_trace, R, T, x_trans, d_min, d_max);

      for (auto it = correlator_electric.begin();
           it != correlator_electric.end(); ++it) {
        stream_electric << T << "," << R << "," << it->first - R / 2 << ","
                        << it->second << "," << wilson_loop_average << ","
                        << plaket_time_average << std::endl;
      }

      for (auto it = correlator_magnetic.begin();
           it != correlator_magnetic.end(); ++it) {
        stream_magnetic << T << "," << R << "," << it->first - R / 2 << ","
                        << it->second << "," << wilson_loop_average << ","
                        << plaket_space_average << std::endl;
      }
    }
  }

  end_time = omp_get_wtime();
  observables_time = end_time - start_time;
  cout << "flux_tube time: " << observables_time << endl;

  stream_electric.close();
  stream_magnetic.close();
}
