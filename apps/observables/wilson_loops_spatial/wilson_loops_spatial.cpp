#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/smearing.h"

#include <ctime>
#include <iostream>
#include <omp.h>

#ifndef MATRIX
#define MATRIX su2
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

  string conf_format;
  string conf_path;
  string path_wilson;
  string representation;
  int L_spat, L_time;
  int APE_steps;
  double alpha;
  int T_min, T_max, R_min, R_max;
  int bytes_skip = 0;
  bool convert = 0;
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "-bytes_skip") {
      bytes_skip = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-conf_path") {
      conf_path = argv[++i];
    } else if (string(argv[i]) == "-convert") {
      istringstream(string(argv[++i])) >> convert;
    } else if (string(argv[i]) == "-L_spat") {
      L_spat = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-L_time") {
      L_time = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-path_wilson") {
      path_wilson = argv[++i];
    } else if (string(argv[i]) == "-representation") {
      representation = argv[++i];
    } else if (string(argv[i]) == "-APE_steps") {
      APE_steps = stoi(argv[++i]);
    } else if (string(argv[i]) == "-alpha") {
      alpha = atof(argv[++i]);
    } else if (string(argv[i]) == "-T_min") {
      T_min = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-T_max") {
      T_max = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-R_min") {
      R_min = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-R_max") {
      R_max = stoi(string(argv[++i]));
    }
  }

  x_size = L_spat;
  y_size = L_spat;
  z_size = L_spat;
  t_size = L_time;

  cout << "conf_format " << conf_format << endl;
  cout << "conf_path " << conf_path << endl;
  cout << "bytes_skip " << bytes_skip << endl;
  cout << "convert " << convert << endl;
  cout << "alpha " << alpha << endl;
  cout << "APE_steps " << APE_steps << endl;
  cout << "L_spat " << L_spat << endl;
  cout << "L_time " << L_time << endl;
  cout << "path_wilson " << path_wilson << endl;
  cout << "representation " << representation << endl;
  cout << "T_min " << T_min << endl;
  cout << "T_max " << T_max << endl;
  cout << "R_min " << R_min << endl;
  cout << "R_max " << R_max << endl;
  cout << endl;

  cout.precision(17);

  data<MATRIX> conf;

  get_data(conf, conf_path, conf_format, bytes_skip, convert);

  cout << "plaket " << plaket(conf.array) << endl;

  ofstream stream_wilson;
  stream_wilson.precision(17);
  // open file
  stream_wilson.open(path_wilson);

  stream_wilson << "time_size,space_size,wilson_loop" << endl;

  vector<vector<MATRIX>> conf_separated = separate_wilson(conf.array);

  conf.array.clear();
  conf.array.shrink_to_fit();
  cout << "plaket space " << plaket_space_parallel(conf_separated) << endl;

  start_time = omp_get_wtime();

  vector<vector<MATRIX>> smeared =
      smearing_APE_2d_initial(conf_separated, alpha);

  for (int step = 1; step < APE_steps; step++) {
    smearing_APE_2d(smeared, alpha);
  }

  map<tuple<int, int>, double> wilson_loops;

  wilson_loops = wilson_spatial_parallel(conf_separated, smeared, R_min, R_max,
                                         T_min, T_max);

  for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
    stream_wilson << get<0>(it->first) << "," << get<1>(it->first) << ","
                  << it->second << endl;
  }

  end_time = omp_get_wtime();
  observables_time = end_time - start_time;
  cout << "wilson loops time: " << observables_time << endl;

  stream_wilson.close();
}
