#include "../../lib/cpu/include/smearing.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/matrix.h"
#include "../../lib/cpu/include/plaket.h"
#include "../../lib/cpu/include/polyakov_loops.h"
#include "../../lib/cpu/include/wilson_loops.h"

#include <ctime>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>

#ifndef MATRIX
#define MATRIX su3
#endif

using namespace std;

int x_size;
int y_size;
int z_size;
int t_size;
int size1;
int size2;

void polyakov_loop_coorelator_write_result(
    map<tuple<int, double>, double> &result,
    const map<double, double> &correlator, int smearing_step) {
  for (auto it : correlator) {
    result[{smearing_step, it.first}] = it.second;
  }
}

void wilson_loops_write_result(map<tuple<int, int, int>, double> &result,
                               const map<tuple<int, int>, double> &wilson_loops,
                               int smearing_step) {
  for (auto it : wilson_loops) {
    result[{smearing_step, get<0>(it.first), get<1>(it.first)}] = it.second;
  }
}

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double smearing_time;
  double observables_time;

  string conf_format;
  string file_precision;
  string conf_path;
  string conf_path_output;
  string path_wilson;
  string path_polyakov_correlator;
  string correlator_type;
  string path_polyakov_loop;
  double HYP_alpha1, HYP_alpha2, HYP_alpha3;
  double APE_alpha;
  bool APE_enabled, HYP_enabled;
  int HYP_steps, APE_steps;
  int L_spat, L_time;
  bool wilson_enabled, polyakov_correlator_enabled, polyakov_loop_enabled;
  int calculation_APE_start, calculation_step_APE;
  int calculation_HYP_start, calculation_step_HYP;
  int T_min, T_max, R_min, R_max;
  int polyakov_correlator_D;
  int bytes_skip = 0;
  bool save_conf = false;
  bool convert = false;
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "--conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "--file_precision") {
      file_precision = argv[++i];
    } else if (string(argv[i]) == "--bytes_skip") {
      bytes_skip = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--conf_path") {
      conf_path = argv[++i];
    } else if (string(argv[i]) == "--conf_path_output") {
      conf_path_output = argv[++i];
    } else if (string(argv[i]) == "--convert") {
      istringstream(string(argv[++i])) >> convert;
    } else if (string(argv[i]) == "--HYP_alpha1") {
      HYP_alpha1 = atof(argv[++i]);
    } else if (string(argv[i]) == "--HYP_alpha2") {
      HYP_alpha2 = atof(argv[++i]);
    } else if (string(argv[i]) == "--HYP_alpha3") {
      HYP_alpha3 = atof(argv[++i]);
    } else if (string(argv[i]) == "--APE_alpha") {
      APE_alpha = atof(argv[++i]);
    } else if (string(argv[i]) == "--APE") {
      APE_enabled = stoi(argv[++i]);
    } else if (string(argv[i]) == "--HYP") {
      HYP_enabled = stoi(argv[++i]);
    } else if (string(argv[i]) == "--APE_enabled") {
      istringstream(string(argv[++i])) >> APE_enabled;
    } else if (string(argv[i]) == "--HYP_enabled") {
      istringstream(string(argv[++i])) >> HYP_enabled;
    } else if (string(argv[i]) == "--APE_steps") {
      APE_steps = stoi(argv[++i]);
    } else if (string(argv[i]) == "--HYP_steps") {
      HYP_steps = stoi(argv[++i]);
    } else if (string(argv[i]) == "--L_spat") {
      L_spat = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--L_time") {
      L_time = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--wilson_enabled") {
      istringstream(string(argv[++i])) >> wilson_enabled;
    } else if (string(argv[i]) == "--polyakov_correlator_enabled") {
      istringstream(string(argv[++i])) >> polyakov_correlator_enabled;
    } else if (string(argv[i]) == "--polyakov_loop_enabled") {
      istringstream(string(argv[++i])) >> polyakov_loop_enabled;
    } else if (string(argv[i]) == "--save_conf") {
      istringstream(string(argv[++i])) >> save_conf;
    } else if (string(argv[i]) == "--path_wilson") {
      path_wilson = argv[++i];
    } else if (string(argv[i]) == "--path_polyakov_correlator") {
      path_polyakov_correlator = argv[++i];
    } else if (string(argv[i]) == "--path_polyakov_loop") {
      path_polyakov_loop = argv[++i];
    } else if (string(argv[i]) == "--correlator_type") {
      correlator_type = argv[++i];
    } else if (string(argv[i]) == "--T_min") {
      T_min = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--T_max") {
      T_max = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--R_min") {
      R_min = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--R_max") {
      R_max = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--polyakov_correlator_D") {
      polyakov_correlator_D = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--calculation_step_APE") {
      calculation_step_APE = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--calculation_APE_start") {
      calculation_APE_start = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--calculation_step_HYP") {
      calculation_step_HYP = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--calculation_HYP_start") {
      calculation_HYP_start = stoi(string(argv[++i]));
    }
  }

  int x_size1 = L_spat;
  int y_size1 = L_spat;
  int z_size1 = L_spat;
  int t_size1 = L_time;

  cout << "conf_format " << conf_format << endl;
  cout << "file_precision " << file_precision << endl;
  cout << "conf_path " << conf_path << endl;
  cout << "conf_path_output " << conf_path_output << endl;
  cout << "bytes_skip " << bytes_skip << endl;
  cout << "convert " << convert << endl;
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
  cout << "path_polyakov_correlator " << path_polyakov_correlator << endl;
  cout << "path_polyakov_loop " << path_polyakov_loop << endl;
  cout << "wilson_enabled " << wilson_enabled << endl;
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

  Data::LatticeData<DataPatternLexicographical, MATRIX> conf(
      {x_size1, y_size1, z_size1, t_size1});
  Data::read_data_convert(conf, conf_path, conf_format, bytes_skip,
                          file_precision, convert);

  // vector<double> plaket_time_tr;
  // double plaket_unsmeared;

  cout << "plaket unsmeared " << plaket(conf) << endl;

  map<tuple<int, int>, double> wilson_loops;
  map<tuple<int, int, int>, double> wilson_loops_result;
  map<tuple<int, int, int>, double> flux_tube;
  std::vector<double> polyakov_correlator_vec;
  std ::map<double, double> polyakov_correlator;
  std::map<std::tuple<int, double>, double> polyakov_correlator_result;
  std::map<int, double> polyakov_loop_result;

  observables_time = 0;
  start_time = omp_get_wtime();
  if (polyakov_correlator_enabled) {
    if (correlator_type == "singlet") {
      polyakov_correlator_vec =
          polyakov_loop_correlator_singlet(conf, polyakov_correlator_D);
      polyakov_correlator = polyakov_average_directions(polyakov_correlator_vec,
                                                        polyakov_correlator_D);
    } else if (correlator_type == "color_average") {
      polyakov_correlator_vec =
          polyakov_loop_correlator_color_average(conf, polyakov_correlator_D);
      polyakov_correlator = polyakov_average_directions(polyakov_correlator_vec,
                                                        polyakov_correlator_D);
    } else {
      cout << "invalid correlator_type" << endl;
    }
    polyakov_loop_coorelator_write_result(polyakov_correlator_result,
                                          polyakov_correlator, 0);
  }
  if (polyakov_loop_enabled) {
    polyakov_loop_result[0] = polyakov_loop(conf);
  }
  end_time = omp_get_wtime();
  observables_time += end_time - start_time;
  if (HYP_enabled == 1) {
    smearing_time = 0;
    for (int HYP_step = 1; HYP_step <= HYP_steps; HYP_step++) {
      start_time = omp_get_wtime();
      smearing_HYP(conf, HYP_alpha1, HYP_alpha2, HYP_alpha3);
      end_time = omp_get_wtime();
      smearing_time += end_time - start_time;

      start_time = omp_get_wtime();
      if ((HYP_step - calculation_HYP_start) % calculation_step_HYP == 0 &&
          HYP_step >= calculation_HYP_start) {
        if (polyakov_correlator_enabled) {
          if (correlator_type == "singlet") {
            polyakov_correlator_vec =
                polyakov_loop_correlator_singlet(conf, polyakov_correlator_D);
            polyakov_correlator = polyakov_average_directions(
                polyakov_correlator_vec, polyakov_correlator_D);
          } else if (correlator_type == "color_average") {
            polyakov_correlator_vec = polyakov_loop_correlator_color_average(
                conf, polyakov_correlator_D);
            polyakov_correlator = polyakov_average_directions(
                polyakov_correlator_vec, polyakov_correlator_D);
          } else {
            cout << "invalid correlator_type" << endl;
          }
          polyakov_loop_coorelator_write_result(polyakov_correlator_result,
                                                polyakov_correlator, HYP_step);
          if (polyakov_loop_enabled) {
            polyakov_loop_result[HYP_step] = polyakov_loop(conf);
          }
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
    wilson_loops = wilson_loop(conf, R_min, R_max, T_min, T_max);
    wilson_loops_write_result(wilson_loops_result, wilson_loops, 0);
  }

  if (APE_enabled == 1) {
    smearing_time = 0;
    for (int APE_step = 1; APE_step <= APE_steps; APE_step++) {

      start_time = omp_get_wtime();

      smearing_APE(conf, APE_alpha);

      end_time = omp_get_wtime();
      smearing_time += end_time - start_time;

      start_time = omp_get_wtime();

      if ((APE_step - calculation_APE_start) % calculation_step_APE == 0 &&
          APE_step >= calculation_APE_start) {

        if (wilson_enabled) {
          wilson_loops = wilson_loop(conf, R_min, R_max, T_min, T_max);
          wilson_loops_write_result(wilson_loops_result, wilson_loops,
                                    APE_step);
        }
      }
      end_time = omp_get_wtime();
      observables_time += end_time - start_time;
    }
    cout << "i=" << APE_steps << " iterations of APE time: " << smearing_time
         << endl;
    cout << "observables time during APE: " << observables_time << endl;
  }

  ofstream stream_wilson;
  if (wilson_enabled) {
    stream_wilson.open(path_wilson);
    stream_wilson.precision(17);
    stream_wilson << "smearing_step,time_size,space_size,wilson_loop" << endl;
    for (auto it : wilson_loops_result) {
      stream_wilson << get<0>(it.first) << "," << get<1>(it.first) << ","
                    << get<2>(it.first) << "," << it.second << endl;
    }
    stream_wilson.close();
  }
  ofstream stream_polyakov_correlator;
  if (polyakov_correlator_enabled) {
    stream_polyakov_correlator.open(path_polyakov_correlator);
    stream_polyakov_correlator.precision(17);
    stream_polyakov_correlator << "smearing_step,distance,correlator" << endl;
    for (auto it : polyakov_correlator_result) {
      stream_polyakov_correlator << std::get<0>(it.first) << ","
                                 << std::get<1>(it.first) << "," << it.second
                                 << endl;
    }
    stream_polyakov_correlator.close();
  }
  ofstream stream_polyakov_loop;
  if (polyakov_loop_enabled) {
    stream_polyakov_loop.open(path_polyakov_loop);
    stream_polyakov_loop.precision(17);
    stream_polyakov_loop << "HYP_step,polyakov_loop" << endl;
    for (auto it : polyakov_loop_result) {
      stream_polyakov_correlator << it.first << "," << it.second << endl;
    }
    stream_polyakov_loop.close();
  }
  if (save_conf) {
    FilePatternLexicographical<4, MATRIX> file_pattern_lexicographical;
    conf.write_data(conf_path_output, file_pattern_lexicographical);
  }
}
