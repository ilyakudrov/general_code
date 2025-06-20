#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/plaket.h"

#include <ctime>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>

#ifndef MATRIX
#define MATRIX su2
#endif

using namespace std;

int x_size;
int y_size;
int z_size;
int t_size;
int size1;
int size2;

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double observables_time;

  string conf_format;
  string file_precision;
  string conf_path;
  string path_wilson;
  string representation;
  string axis;
  int L_spat, L_time;
  int T_min, T_max, R_min, R_max;
  int bytes_skip = 0;
  bool convert = 0;
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "--conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "--file_precision") {
      file_precision = argv[++i];
    } else if (string(argv[i]) == "--bytes_skip") {
      bytes_skip = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--conf_path") {
      conf_path = argv[++i];
    } else if (string(argv[i]) == "--convert") {
      istringstream(string(argv[++i])) >> convert;
    } else if (string(argv[i]) == "--L_spat") {
      L_spat = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--L_time") {
      L_time = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--path_wilson") {
      path_wilson = argv[++i];
    } else if (string(argv[i]) == "--representation") {
      representation = argv[++i];
    } else if (string(argv[i]) == "--axis") {
      axis = argv[++i];
    } else if (string(argv[i]) == "--T_min") {
      T_min = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--T_max") {
      T_max = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--R_min") {
      R_min = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--R_max") {
      R_max = stoi(string(argv[++i]));
    }
  }

  int x_size1 = L_spat;
  int y_size1 = L_spat;
  int z_size1 = L_spat;
  int t_size1 = L_time;

  cout << "conf_format " << conf_format << endl;
  cout << "file_precision " << file_precision << endl;
  cout << "conf_path " << conf_path << endl;
  cout << "bytes_skip " << bytes_skip << endl;
  cout << "convert " << convert << endl;
  cout << "L_spat " << L_spat << endl;
  cout << "L_time " << L_time << endl;
  cout << "path_wilson " << path_wilson << endl;
  cout << "representation " << representation << endl;
  cout << "axis " << axis << endl;
  cout << "T_min " << T_min << endl;
  cout << "T_max " << T_max << endl;
  cout << "R_min " << R_min << endl;
  cout << "R_max " << R_max << endl;
  cout << endl;

  cout.precision(17);

  Data::LatticeData<DataPatternLexicographical, MATRIX> conf(
      {x_size1, y_size1, z_size1, t_size1});
  Data::read_data_convert(conf, conf_path, conf_format, bytes_skip,
                          file_precision, convert);

  cout << "plaket " << plaket(conf) << endl;

  vector<vector<MATRIX>> conf_separated;

  start_time = omp_get_wtime();

  map<tuple<int, int>, double> wilson_loops_onaxis;
  map<tuple<int, double>, double> wilson_loops_offaxis;

  if (axis.compare("on-axis") == 0) {
    cout << "on-axis" << endl;
    if (representation.compare("fundamental") == 0) {
      cout << "fundamental" << endl;
      wilson_loops_onaxis = wilson_loop(conf.array, R_min, R_max, T_min, T_max);
    } else if (representation.compare("adjoint") == 0) {
      cout << "adjoint" << endl;
      wilson_loops_onaxis =
          wilson_loop_adjoint(conf.array, R_min, R_max, T_min, T_max);
    } else {
      cout << "wrong representation" << endl;
    }
    ofstream stream_wilson;
    stream_wilson.precision(17);
    stream_wilson.open(path_wilson);
    stream_wilson << "time_size,space_size,wilson_loop" << endl;
    stream_wilson.close();
    for (auto it = wilson_loops_onaxis.begin(); it != wilson_loops_onaxis.end();
         it++) {
      stream_wilson << get<0>(it->first) << "," << get<1>(it->first) << ","
                    << it->second << endl;
    }
  } else if (axis.compare("off-axis") == 0) {
    cout << "off-axis" << endl;
    if (representation.compare("fundamental") == 0) {
      cout << "fundamental" << endl;
      wilson_loops_offaxis = wilson_offaxis_result(conf.array, R_min - 0.1,
                                                   R_max + 0.1, T_min, T_max);
    } else if (representation.compare("adjoint") == 0) {
      cout << "adjoint" << endl;
      wilson_loops_offaxis = wilson_offaxis_adjoint_result(
          conf.array, R_min - 0.1, R_max + 0.1, T_min, T_max);
    } else {
      cout << "wrong representation" << endl;
    }
    ofstream stream_wilson;
    stream_wilson.precision(17);
    stream_wilson.open(path_wilson);
    stream_wilson << "time_size,space_size,wilson_loop" << endl;
    stream_wilson.close();
    for (auto it = wilson_loops_offaxis.begin();
         it != wilson_loops_offaxis.end(); it++) {
      stream_wilson << get<0>(it->first) << "," << get<1>(it->first) << ","
                    << it->second << endl;
    }
  } else {
    cout << "wrong axis" << endl;
  }

  end_time = omp_get_wtime();
  observables_time = end_time - start_time;
  cout << "wilson loops time: " << observables_time << endl;
}
