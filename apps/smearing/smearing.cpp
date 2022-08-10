#include "../../lib/cpu/include/smearing.h"
#include "../../lib/cpu/include/basic_observables.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/flux_tube.h"
#include "../../lib/cpu/include/matrix.h"

#include <ctime>
#include <iostream>
#include <omp.h>

using namespace std;

int x_size;
int y_size;
int z_size;
int t_size;

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double smearing_time;
  double observables_time;

  string conf_format_wilson;
  string conf_path_wilson;
  string conf_format_plaket;
  string conf_path_plaket;
  string path_wilson;
  string path_flux;
  double HYP_alpha1, HYP_alpha2, HYP_alpha3;
  double APE_alpha;
  bool APE_enabled, HYP_enabled;
  int HYP_steps, APE_steps;
  int L_spat, L_time;
  bool wilson_enabled, flux_enabled;
  int calculation_APE_start, calculation_step_APE;
  int T_min, T_max, R_min, R_max;
  int bytes_skip_wilson = 0;
  int bytes_skip_plaket = 0;
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-conf_format_wilson") {
      conf_format_wilson = argv[++i];
    } else if (string(argv[i]) == "-bytes_skip_wilson") {
      bytes_skip_wilson = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-conf_path_wilson") {
      conf_path_wilson = argv[++i];
    } else if (string(argv[i]) == "-conf_format_plaket") {
      conf_format_plaket = argv[++i];
    } else if (string(argv[i]) == "-bytes_skip_plaket") {
      bytes_skip_plaket = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-conf_path_plaket") {
      conf_path_plaket = argv[++i];
    } else if (string(argv[i]) == "-HYP_alpha1") {
      HYP_alpha1 = atof(argv[++i]);
    } else if (string(argv[i]) == "-HYP_alpha2") {
      HYP_alpha2 = atof(argv[++i]);
    } else if (string(argv[i]) == "-HYP_alpha3") {
      HYP_alpha3 = atof(argv[++i]);
    } else if (string(argv[i]) == "-APE_alpha") {
      APE_alpha = atof(argv[++i]);
    } else if (string(argv[i]) == "-APE") {
      APE_enabled = stoi(argv[++i]);
    } else if (string(argv[i]) == "-HYP") {
      HYP_enabled = stoi(argv[++i]);
    } else if (string(argv[i]) == "-APE_enabled") {
      istringstream(string(argv[++i])) >> APE_enabled;
    } else if (string(argv[i]) == "-HYP_enabled") {
      istringstream(string(argv[++i])) >> HYP_enabled;
    } else if (string(argv[i]) == "-APE_steps") {
      APE_steps = stoi(argv[++i]);
    } else if (string(argv[i]) == "-HYP_steps") {
      HYP_steps = stoi(argv[++i]);
    } else if (string(argv[i]) == "-L_spat") {
      L_spat = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-L_time") {
      L_time = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-wilson_enabled") {
      istringstream(string(argv[++i])) >> wilson_enabled;
    } else if (string(argv[i]) == "-flux_enabled") {
      istringstream(string(argv[++i])) >> flux_enabled;
    } else if (string(argv[i]) == "-path_wilson") {
      path_wilson = argv[++i];
    } else if (string(argv[i]) == "-path_flux") {
      path_flux = argv[++i];
    } else if (string(argv[i]) == "-T_min") {
      T_min = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-T_max") {
      T_max = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-R_min") {
      R_min = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-R_max") {
      R_max = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-calculation_step_APE") {
      calculation_step_APE = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-calculation_APE_start") {
      calculation_APE_start = stoi(string(argv[++i]));
    }
  }

  x_size = L_spat;
  y_size = L_spat;
  z_size = L_spat;
  t_size = L_time;

  cout << "conf_format_wilson " << conf_format_wilson << endl;
  cout << "conf_path_wilson " << conf_path_wilson << endl;
  cout << "bytes_skip_wilson " << bytes_skip_wilson << endl;
  cout << "conf_format_plaket " << conf_format_plaket << endl;
  cout << "conf_path_plaket " << conf_path_plaket << endl;
  cout << "bytes_skip_plaket " << bytes_skip_plaket << endl;
  cout << "HYP_alpha1 " << HYP_alpha1 << endl;
  cout << "HYP_alpha2 " << HYP_alpha2 << endl;
  cout << "HYP_alpha3 " << HYP_alpha3 << endl;
  cout << "APE_alpha " << APE_alpha << endl;
  cout << "APE_enabled " << APE_enabled << endl;
  cout << "HYP_enabled " << HYP_enabled << endl;
  cout << "APE_steps " << APE_steps << endl;
  cout << "HYP_steps " << HYP_steps << endl;
  cout << "L_spat " << L_spat << endl;
  cout << "L_time " << L_time << endl;
  cout << "path_wilson " << path_wilson << endl;
  cout << "path_flux " << path_flux << endl;
  cout << "wilson_enabled " << wilson_enabled << endl;
  cout << "flux_enabled " << flux_enabled << endl;
  cout << "T_min " << T_min << endl;
  cout << "T_max " << T_max << endl;
  cout << "R_min " << R_min << endl;
  cout << "R_max " << R_max << endl;
  cout << "calculation_step_APE " << calculation_step_APE << endl;
  cout << "calculation_APE_start " << calculation_APE_start << endl;
  cout << endl;

  cout.precision(17);

  data<MATRIX_WILSON> conf_wilson;
  data<MATRIX_PLAKET> conf_plaket;

  vector<double> plaket_time_tr;

  if (flux_enabled) {
    if (string(conf_format_plaket) == "float") {
      conf_plaket.read_float(conf_path_plaket, bytes_skip_plaket);
    } else if (string(conf_format_plaket) == "double") {
      conf_plaket.read_double(conf_path_plaket, bytes_skip_plaket);
    } else if (string(conf_format_plaket) == "double_qc2dstag") {
      conf_plaket.read_double_qc2dstag(conf_path_plaket);
    } else {
      cout << "wrong conf format: " << conf_format_plaket << endl;
      return 0;
    }

    cout << "plaket plaket unsmeared " << plaket(conf_plaket.array) << endl;

    vector<vector<MATRIX_PLAKET>> conf_separated_plaket =
        separate_wilson(conf_plaket.array);

    conf_plaket.array.clear();
    conf_plaket.array.shrink_to_fit();

    plaket_time_tr = plaket_aver_tr_time(conf_separated_plaket);

    for (int i = 0; i < conf_separated_plaket.size(); i++) {
      conf_separated_plaket[i].clear();
      conf_separated_plaket[i].shrink_to_fit();
    }
  }

  if (string(conf_format_wilson) == "float") {
    conf_wilson.read_float(conf_path_wilson, bytes_skip_wilson);
  } else if (string(conf_format_wilson) == "double") {
    conf_wilson.read_double(conf_path_wilson, bytes_skip_wilson);
  } else if (string(conf_format_wilson) == "double_qc2dstag") {
    conf_wilson.read_double_qc2dstag(conf_path_wilson);
  } else {
    cout << "wrong conf format: " << conf_format_wilson << endl;
    return 0;
  }

  cout << "wilson plaket unsmeared " << plaket(conf_wilson.array) << endl;

  ofstream stream_wilson;
  // open file
  if (wilson_enabled) {
    stream_wilson.open(path_wilson);

    stream_wilson.precision(17);

    stream_wilson << "smearing_step,time_size,space_size,wilson_loop" << endl;
  }
  ofstream stream_flux;
  // open file
  if (flux_enabled) {
    stream_flux.open(path_flux);

    stream_flux.precision(17);

    stream_flux << "smearing_step,time_size,space_size,x_long,correlator,"
                   "wilson_loop,plaket"
                << endl;
  }

  map<tuple<int, int>, double> wilson_loops;
  map<tuple<int, int, int>, double> flux_tube;

  vector<vector<MATRIX_WILSON>> conf_separated =
      separate_wilson(conf_wilson.array);

  conf_wilson.array.clear();
  conf_wilson.array.shrink_to_fit();

  if (HYP_enabled == 1) {
    start_time = omp_get_wtime();
    for (int HYP_step = 0; HYP_step < HYP_steps; HYP_step++) {

      smearing_HYP_new(conf_separated, HYP_alpha1, HYP_alpha2, HYP_alpha3);
    }

    end_time = omp_get_wtime();
    smearing_time = end_time - start_time;
    cout << "i=" << HYP_steps << " iterations of HYP time: " << smearing_time
         << endl;
  }

  if (APE_enabled == 1) {

    smearing_time = 0;
    observables_time = 0;

    for (int APE_step = 0; APE_step < APE_steps; APE_step++) {

      start_time = omp_get_wtime();

      smearing_APE_new(conf_separated, APE_alpha);

      end_time = omp_get_wtime();
      smearing_time += end_time - start_time;

      start_time = omp_get_wtime();

      if (APE_step % calculation_step_APE == 0 &&
          APE_step >= calculation_APE_start) {

        if (wilson_enabled) {
          wilson_loops =
              wilson_parallel(conf_separated, R_min, R_max, T_min, T_max);

          for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
            stream_wilson << APE_step + 1 << "," << get<0>(it->first) << ","
                          << get<1>(it->first) << "," << it->second << endl;
          }
        }

        if (flux_enabled) {

          flux_tube = wilson_plaket_correlator(
              plaket_time_tr, conf_separated, 6, L_time / 2, L_spat / 4,
              L_spat / 2, 5, 0, "longitudinal");

          for (auto it = flux_tube.begin(); it != flux_tube.end(); it++) {
            stream_flux << APE_step + 1 << "," << get<0>(it->first) << ","
                        << get<1>(it->first) << "," << get<2>(it->first) << ","
                        << it->second << endl;
          }
        }
      }

      end_time = omp_get_wtime();
      observables_time += end_time - start_time;
    }

    cout << "i=" << APE_steps << " iterations of APE time: " << smearing_time
         << endl;
    cout << "observables time during APE: " << observables_time << endl;
  }

  if (wilson_enabled) {
    stream_wilson.close();
  }
  if (flux_enabled) {
    stream_flux.close();
  }
}
