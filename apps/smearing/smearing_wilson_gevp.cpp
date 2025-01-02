#include "../../lib/cpu/include/basic_observables.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/matrix.h"
#include "../../lib/cpu/include/smearing.h"

#include <ctime>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>

#ifndef MATRIX_WILSON
#define MATRIX_WILSON su2
#endif

using namespace std;

int x_size;
int y_size;
int z_size;
int t_size;

void write_wilson_loops(map<tuple<int, int>, double> &wilson_tmp,
                        map<tuple<int, int, int, int>, double> &wilson_loops,
                        int smearing1, int smearing2) {
  for (auto it = wilson_tmp.begin(); it != wilson_tmp.end(); it++) {
    wilson_loops[tuple<int, int, int, int>(
        smearing1, smearing2, get<0>(it->first), get<1>(it->first))] +=
        it->second;
  }
}

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double smearing_APE_time;
  double smearing_HYP_time;
  double observables_time;

  string conf_format_wilson;
  string conf_path_wilson;
  string path_wilson;
  string representation;
  double HYP_alpha1, HYP_alpha2, HYP_alpha3;
  double APE_alpha;
  bool HYP_enabled;
  int HYP_steps, APE_steps;
  int L_spat, L_time;
  int calculation_APE_start, calculation_step_APE;
  int T_min, T_max, R_min, R_max;
  int bytes_skip_wilson = 0;
  bool convert_wilson = false;
  int N_dir = 1;
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-conf_format_wilson") {
      conf_format_wilson = argv[++i];
    } else if (string(argv[i]) == "-bytes_skip_wilson") {
      bytes_skip_wilson = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-conf_path_wilson") {
      conf_path_wilson = argv[++i];
    } else if (string(argv[i]) == "-representation") {
      representation = argv[++i];
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
    } else if (string(argv[i]) == "-HYP") {
      HYP_enabled = stoi(argv[++i]);
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
    } else if (string(argv[i]) == "-path_wilson") {
      path_wilson = argv[++i];
    } else if (string(argv[i]) == "-T_min") {
      T_min = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-T_max") {
      T_max = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-R_min") {
      R_min = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-R_max") {
      R_max = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-N_dir") {
      N_dir = stoi(string(argv[++i]));
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
  cout << "representation " << representation << endl;
  cout << "bytes_skip_wilson " << bytes_skip_wilson << endl;
  cout << "convert_wilson " << convert_wilson << endl;
  cout << "HYP_alpha1 " << HYP_alpha1 << endl;
  cout << "HYP_alpha2 " << HYP_alpha2 << endl;
  cout << "HYP_alpha3 " << HYP_alpha3 << endl;
  cout << "APE_alpha " << APE_alpha << endl;
  cout << "HYP_enabled " << HYP_enabled << endl;
  cout << "APE_steps " << APE_steps << endl;
  cout << "HYP_steps " << HYP_steps << endl;
  cout << "L_spat " << L_spat << endl;
  cout << "L_time " << L_time << endl;
  cout << "path_wilson " << path_wilson << endl;
  cout << "T_min " << T_min << endl;
  cout << "T_max " << T_max << endl;
  cout << "R_min " << R_min << endl;
  cout << "R_max " << R_max << endl;
  cout << "N_dir " << N_dir << endl;
  cout << "calculation_step_APE " << calculation_step_APE << endl;
  cout << "calculation_APE_start " << calculation_APE_start << endl;
  cout << endl;

  cout.precision(17);

  Data::data<MATRIX_WILSON> conf_wilson;
  map<tuple<int, int>, double> wilson_tmp;
  map<tuple<int, int, int, int>, double> wilson_loops;
  vector<vector<MATRIX_WILSON>> conf_separated1;
  vector<vector<MATRIX_WILSON>> conf_separated2;

  smearing_APE_time = 0;
  smearing_HYP_time = 0;
  observables_time = 0;

  link1 link(x_size, y_size, z_size, t_size);

  for (int dir = 0; dir < N_dir; dir++) {
    get_data(conf_wilson, conf_path_wilson, conf_format_wilson,
             bytes_skip_wilson, convert_wilson);
    if (dir > 0) {
      conf_wilson.array = swap_directions(conf_wilson.array, dir - 1, 3);
    }
    conf_separated1 = separate_wilson(conf_wilson.array);
    conf_wilson.array.clear();
    conf_wilson.array.shrink_to_fit();

    conf_separated2 = conf_separated1;
    if (HYP_enabled == 1) {
      start_time = omp_get_wtime();
      for (int HYP_step = 1; HYP_step <= HYP_steps; HYP_step++) {
        smearing_HYP_parallel(conf_separated1, HYP_alpha1, HYP_alpha2,
                              HYP_alpha3);
      }
      end_time = omp_get_wtime();
      smearing_HYP_time += end_time - start_time;
      conf_separated2[3] = conf_separated1[3];
    }

    // wilson loops at (0, 0) APE_steps
    start_time = omp_get_wtime();
    if (representation == "fundamental") {
      wilson_tmp = wilson_parallel(conf_separated1, R_min, R_max, T_min, T_max);
    } else if (representation == "adjoint") {
      wilson_tmp =
          wilson_adjoint_parallel(conf_separated1, R_min, R_max, T_min, T_max);
    } else {
      cout << "wrong representation" << endl;
    }
    write_wilson_loops(wilson_tmp, wilson_loops, 0, 0);
    end_time = omp_get_wtime();
    observables_time += end_time - start_time;

    // wilson loops at (0, APE_step2) APE_steps
    for (int APE_step = 1; APE_step <= APE_steps; APE_step++) {
      start_time = omp_get_wtime();
      smearing_APE_parallel(conf_separated2, APE_alpha);
      end_time = omp_get_wtime();
      smearing_APE_time += end_time - start_time;

      if ((APE_step - calculation_APE_start) % calculation_step_APE == 0 &&
          APE_step >= calculation_APE_start) {
        start_time = omp_get_wtime();
        if (representation == "fundamental") {
          wilson_tmp = wilson_gevp_parallel(conf_separated1, conf_separated2,
                                            R_min, R_max, T_min, T_max);
        } else if (representation == "adjoint") {
          wilson_tmp = wilson_gevp_adjoint_parallel(
              conf_separated1, conf_separated2, R_min, R_max, T_min, T_max);
        } else {
          cout << "wrong representation" << endl;
        }
        write_wilson_loops(wilson_tmp, wilson_loops, 0, APE_step);
        end_time = omp_get_wtime();
        observables_time += end_time - start_time;
      }
    }

    // wilson loops at (APE_step1, APE_step2) APE_steps, APE_step1 < APE_step2
    for (int APE_step1 = 1; APE_step1 <= APE_steps; APE_step1++) {
      start_time = omp_get_wtime();
      smearing_APE_parallel(conf_separated1, APE_alpha);
      end_time = omp_get_wtime();
      smearing_APE_time += end_time - start_time;

      if ((APE_step1 - calculation_APE_start) % calculation_step_APE == 0 &&
          APE_step1 >= calculation_APE_start) {
        start_time = omp_get_wtime();
        if (representation == "fundamental") {
          wilson_tmp =
              wilson_parallel(conf_separated1, R_min, R_max, T_min, T_max);
        } else if (representation == "adjoint") {
          wilson_tmp = wilson_adjoint_parallel(conf_separated1, R_min, R_max,
                                               T_min, T_max);
        } else {
          cout << "wrong representation" << endl;
        }
        write_wilson_loops(wilson_tmp, wilson_loops, APE_step1, APE_step1);
        end_time = omp_get_wtime();
        observables_time += end_time - start_time;

        conf_separated2 = conf_separated1;
        for (int APE_step2 = APE_step1 + 1; APE_step2 <= APE_steps;
             APE_step2++) {
          start_time = omp_get_wtime();
          smearing_APE_parallel(conf_separated2, APE_alpha);
          end_time = omp_get_wtime();
          smearing_APE_time += end_time - start_time;

          if ((APE_step2 - APE_step1) % calculation_step_APE == 0) {
            start_time = omp_get_wtime();
            if (representation == "fundamental") {
              wilson_tmp = wilson_gevp_parallel(
                  conf_separated1, conf_separated2, R_min, R_max, T_min, T_max);
            } else if (representation == "adjoint") {
              wilson_tmp = wilson_gevp_adjoint_parallel(
                  conf_separated1, conf_separated2, R_min, R_max, T_min, T_max);
            } else {
              cout << "wrong representation" << endl;
            }
            write_wilson_loops(wilson_tmp, wilson_loops, APE_step1, APE_step2);
            end_time = omp_get_wtime();
            observables_time += end_time - start_time;
          }
        }
      }
    }
  }

  for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
    it->second = it->second / N_dir;
  }

  ofstream stream_wilson;
  stream_wilson.precision(17);
  // open file
  stream_wilson.open(path_wilson);
  stream_wilson.precision(17);
  stream_wilson
      << "smearing_step1,smearing_step2,time_size,space_size,wilson_loop"
      << endl;

  for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
    stream_wilson << get<0>(it->first) << "," << get<1>(it->first) << ","
                  << get<2>(it->first) << "," << get<3>(it->first) << ","
                  << it->second << endl;
  }

  cout << "HYP time: " << smearing_HYP_time << endl;
  cout << "APE time: " << smearing_APE_time << endl;
  cout << "observables time: " << observables_time << endl;

  stream_wilson.close();
}
