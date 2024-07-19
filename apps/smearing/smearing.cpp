#include "../../lib/cpu/include/smearing.h"
#include "../../lib/cpu/include/basic_observables.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/decomposition.h"
#include "../../lib/cpu/include/flux_tube.h"
#include "../../lib/cpu/include/matrix.h"

#include <ctime>
#include <iostream>
#include <omp.h>

#include "sys/sysinfo.h"
#include "sys/types.h"

struct sysinfo memInfo;

#ifndef MATRIX_WILSON
#define MATRIX_WILSON su2
#endif

#ifndef MATRIX_PLAKET
#define MATRIX_PLAKET su2
#endif

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
  string conf_path_output;
  string conf_format_plaket;
  string conf_path_plaket;
  string path_wilson;
  string path_flux;
  string path_polyakov_correlator;
  string correlator_type;
  string path_polyakov_loop;
  double HYP_alpha1, HYP_alpha2, HYP_alpha3;
  double APE_alpha;
  bool APE_enabled, HYP_enabled;
  int HYP_steps, APE_steps;
  int L_spat, L_time;
  bool wilson_enabled, flux_enabled, polyakov_correlator_enabled,
      polyakov_loop_enabled;
  int calculation_APE_start, calculation_step_APE;
  int calculation_HYP_start, calculation_step_HYP;
  int T_min, T_max, R_min, R_max;
  int polyakov_correlator_D;
  int bytes_skip_wilson = 0;
  int bytes_skip_plaket = 0;
  bool save_conf = false;
  bool convert_plaket = false;
  bool convert_wilson = false;
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-conf_format_wilson") {
      conf_format_wilson = argv[++i];
    } else if (string(argv[i]) == "-bytes_skip_wilson") {
      bytes_skip_wilson = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-conf_path_wilson") {
      conf_path_wilson = argv[++i];
    } else if (string(argv[i]) == "-conf_path_output") {
      conf_path_output = argv[++i];
    } else if (string(argv[i]) == "-conf_format_plaket") {
      conf_format_plaket = argv[++i];
    } else if (string(argv[i]) == "-bytes_skip_plaket") {
      bytes_skip_plaket = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-conf_path_plaket") {
      conf_path_plaket = argv[++i];
    } else if (string(argv[i]) == "-convert_plaket") {
      istringstream(string(argv[++i])) >> convert_plaket;
    } else if (string(argv[i]) == "-convert_wilson") {
      istringstream(string(argv[++i])) >> convert_wilson;
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
    } else if (string(argv[i]) == "-polyakov_correlator_enabled") {
      istringstream(string(argv[++i])) >> polyakov_correlator_enabled;
    } else if (string(argv[i]) == "-polyakov_loop_enabled") {
      istringstream(string(argv[++i])) >> polyakov_loop_enabled;
    } else if (string(argv[i]) == "-save_conf") {
      istringstream(string(argv[++i])) >> save_conf;
    } else if (string(argv[i]) == "-path_wilson") {
      path_wilson = argv[++i];
    } else if (string(argv[i]) == "-path_flux") {
      path_flux = argv[++i];
    } else if (string(argv[i]) == "-path_polyakov_correlator") {
      path_polyakov_correlator = argv[++i];
    } else if (string(argv[i]) == "-path_polyakov_loop") {
      path_polyakov_loop = argv[++i];
    } else if (string(argv[i]) == "-correlator_type") {
      correlator_type = argv[++i];
    } else if (string(argv[i]) == "-T_min") {
      T_min = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-T_max") {
      T_max = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-R_min") {
      R_min = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-R_max") {
      R_max = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-polyakov_correlator_D") {
      polyakov_correlator_D = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-calculation_step_APE") {
      calculation_step_APE = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-calculation_APE_start") {
      calculation_APE_start = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-calculation_step_HYP") {
      calculation_step_HYP = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-calculation_HYP_start") {
      calculation_HYP_start = stoi(string(argv[++i]));
    }
  }

  x_size = L_spat;
  y_size = L_spat;
  z_size = L_spat;
  t_size = L_time;

  cout << "conf_format_wilson " << conf_format_wilson << endl;
  cout << "conf_path_wilson " << conf_path_wilson << endl;
  cout << "conf_path_output " << conf_path_output << endl;
  cout << "bytes_skip_wilson " << bytes_skip_wilson << endl;
  cout << "conf_format_plaket " << conf_format_plaket << endl;
  cout << "conf_path_plaket " << conf_path_plaket << endl;
  cout << "bytes_skip_plaket " << bytes_skip_plaket << endl;
  cout << "convert_plaket " << convert_plaket << endl;
  cout << "convert_wilson " << convert_wilson << endl;
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
  cout << "path_polyakov_correlator " << path_polyakov_correlator << endl;
  cout << "path_polyakov_loop " << path_polyakov_loop << endl;
  cout << "wilson_enabled " << wilson_enabled << endl;
  cout << "flux_enabled " << flux_enabled << endl;
  cout << "polyakov_correlator_enabled " << polyakov_correlator_enabled << endl;
  cout << "polyakov_loop_enabled " << polyakov_loop_enabled << endl;
  cout << "correlator_type " << correlator_type << endl;
  cout << "save_conf " << save_conf << endl;
  cout << "T_min " << T_min << endl;
  cout << "T_max " << T_max << endl;
  cout << "R_min " << R_min << endl;
  cout << "R_max " << R_max << endl;
  cout << "polyakov_correlator_D " << polyakov_correlator_D << endl;
  cout << "calculation_step_APE " << calculation_step_APE << endl;
  cout << "calculation_APE_start " << calculation_APE_start << endl;
  cout << "calculation_step_HYP " << calculation_step_HYP << endl;
  cout << "calculation_HYP_start " << calculation_HYP_start << endl;
  cout << endl;

  cout.precision(17);

  data<MATRIX_WILSON> conf_wilson;
  data<MATRIX_PLAKET> conf_plaket;

  vector<double> plaket_time_tr;
  double plaket_unsmeared;

  if (flux_enabled) {
    get_data(conf_plaket, conf_path_plaket, conf_format_plaket,
             bytes_skip_plaket, convert_plaket);

    plaket_unsmeared = plaket(conf_plaket.array);

    cout << "plaket plaket unsmeared " << plaket_unsmeared << endl;

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

  get_data(conf_wilson, conf_path_wilson, conf_format_wilson, bytes_skip_wilson,
           convert_wilson);

  cout << "wilson plaket unsmeared " << plaket(conf_wilson.array) << endl;

  ofstream stream_wilson;
  stream_wilson.precision(17);
  // open file
  if (wilson_enabled) {
    stream_wilson.open(path_wilson);

    stream_wilson.precision(17);

    stream_wilson << "smearing_step,time_size,space_size,wilson_loop" << endl;
  }
  ofstream stream_flux;
  stream_flux.precision(17);
  // open file
  if (flux_enabled) {
    stream_flux.open(path_flux);

    stream_flux.precision(17);

    stream_flux << "smearing_step,time_size,space_size,x_long,correlator,"
                   "wilson_loop,plaket"
                << endl;
  }
  ofstream stream_polyakov_correlator;
  stream_polyakov_correlator.precision(17);
  // open file
  if (polyakov_correlator_enabled) {
    stream_polyakov_correlator.open(path_polyakov_correlator);
    stream_polyakov_correlator.precision(17);
    stream_polyakov_correlator << "smearing_step,distance,correlator" << endl;
  }
  ofstream stream_polyakov_loop;
  stream_polyakov_correlator.precision(17);
  if (polyakov_loop_enabled) {
    stream_polyakov_loop.open(path_polyakov_loop);
    stream_polyakov_loop.precision(17);
    stream_polyakov_loop << "HYP_step,polyakov_loop" << endl;
  }

  map<tuple<int, int>, double> wilson_loops;
  map<tuple<int, int, int>, double> flux_tube;
  std::vector<double> polyakov_correlator_vec;
  std ::map<double, double> polyakov_correlator;

  vector<vector<MATRIX_WILSON>> conf_separated =
      separate_wilson(conf_wilson.array);
  conf_wilson.array.clear();
  conf_wilson.array.shrink_to_fit();

  observables_time = 0;
  start_time = omp_get_wtime();
  if (polyakov_correlator_enabled) {
    if (correlator_type == "singlet") {
      polyakov_correlator_vec = polyakov_loop_correlator_singlet(
          conf_separated, polyakov_correlator_D);
      polyakov_correlator = polyakov_average_directions(polyakov_correlator_vec,
                                                        polyakov_correlator_D);
    } else if (correlator_type == "color_average") {
      polyakov_correlator_vec =
          polyakov_loop_correlator(conf_separated, polyakov_correlator_D);
      polyakov_correlator = polyakov_average_directions(polyakov_correlator_vec,
                                                        polyakov_correlator_D);
    } else {
      cout << "invalid correlator_type" << endl;
    }
    for (auto it = polyakov_correlator.begin(); it != polyakov_correlator.end();
         it++) {
      stream_polyakov_correlator << 0 << "," << it->first << "," << it->second
                                 << std::endl;
    }
  }
  if (polyakov_loop_enabled) {
    stream_polyakov_loop << 0 << "," << polyakov_loop_parallel(conf_separated)
                         << std::endl;
  }
  end_time = omp_get_wtime();
  observables_time += end_time - start_time;

  if (HYP_enabled == 1) {
    smearing_time = 0;
    for (int HYP_step = 1; HYP_step <= HYP_steps; HYP_step++) {
      start_time = omp_get_wtime();
      smearing_HYP_parallel(conf_separated, HYP_alpha1, HYP_alpha2, HYP_alpha3);
      end_time = omp_get_wtime();
      smearing_time += end_time - start_time;

      start_time = omp_get_wtime();
      if ((HYP_step - calculation_HYP_start) % calculation_step_HYP == 0 &&
          HYP_step >= calculation_HYP_start) {
        if (polyakov_correlator_enabled) {
          if (correlator_type == "singlet") {
            polyakov_correlator_vec = polyakov_loop_correlator_singlet(
                conf_separated, polyakov_correlator_D);
            polyakov_correlator = polyakov_average_directions(
                polyakov_correlator_vec, polyakov_correlator_D);
          } else if (correlator_type == "color_average") {
            polyakov_correlator_vec =
                polyakov_loop_correlator(conf_separated, polyakov_correlator_D);
            polyakov_correlator = polyakov_average_directions(
                polyakov_correlator_vec, polyakov_correlator_D);
          } else {
            cout << "invalid correlator_type" << endl;
          }
          for (auto it = polyakov_correlator.begin();
               it != polyakov_correlator.end(); it++) {
            stream_polyakov_correlator << HYP_step << "," << it->first << ","
                                       << it->second << std::endl;
          }
        }
        if (polyakov_loop_enabled) {
          stream_polyakov_loop << HYP_step << ","
                               << polyakov_loop_parallel(conf_separated)
                               << std::endl;
        }
      }
      end_time = omp_get_wtime();
      observables_time += end_time - start_time;
    }
    cout << "i=" << HYP_steps << " iterations of HYP time: " << smearing_time
         << endl;
    cout << "HYP observables time: " << observables_time << endl;
  }

  if (wilson_enabled) {
    wilson_loops = wilson_parallel(conf_separated, R_min, R_max, T_min, T_max);

    for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
      stream_wilson << 0 << "," << get<0>(it->first) << "," << get<1>(it->first)
                    << "," << it->second << endl;
    }
  }

  if (APE_enabled == 1) {
    smearing_time = 0;
    for (int APE_step = 1; APE_step <= APE_steps; APE_step++) {

      start_time = omp_get_wtime();

      smearing_APE_parallel(conf_separated, APE_alpha);

      end_time = omp_get_wtime();
      smearing_time += end_time - start_time;

      start_time = omp_get_wtime();

      if ((APE_step - calculation_APE_start) % calculation_step_APE == 0 &&
          APE_step >= calculation_APE_start) {

        if (wilson_enabled) {
          wilson_loops =
              wilson_parallel(conf_separated, R_min, R_max, T_min, T_max);

          for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
            stream_wilson << APE_step << "," << get<0>(it->first) << ","
                          << get<1>(it->first) << "," << it->second << endl;
          }
        }

        if (flux_enabled) {

          flux_tube = wilson_plaket_correlator(
              plaket_time_tr, conf_separated, 2, L_time / 2, L_spat / 4,
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
  if (polyakov_correlator_enabled) {
    stream_polyakov_correlator.close();
  }

  if (save_conf) {
    conf_wilson.array = merge_wilson(conf_separated);
    conf_wilson.write_double(conf_path_output);
  }
}
