#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/matrix.h"

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

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double observables_time;

  string conf_format;
  string conf_path;
  string path;
  int L_spat, L_time;
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
    } else if (string(argv[i]) == "-path") {
      path = argv[++i];
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
  cout << "L_spat " << L_spat << endl;
  cout << "L_time " << L_time << endl;
  cout << "path " << path << endl;
  cout << endl;

  cout.precision(17);

  data<MATRIX> conf;

  get_data(conf, conf_path, conf_format, bytes_skip, convert);

  ofstream stream;
  stream.precision(17);
  // open file
  stream.open(path);

  stream << "plaket,plaket_space,plaket_time" << endl;

  // vector<vector<MATRIX>> conf_separated;
  // conf.array.clear();
  // conf.array.shrink_to_fit();

  start_time = omp_get_wtime();

  // double plaket = plaket_parallel(conf_separated);
  // double plaket_space = plaket_space_parallel(conf_separated);
  // double plaket_time = plaket_time_parallel(conf_separated);
  double plaket1 = plaket(conf.array);
  double plaket_space1 = plaket_space(conf.array);
  double plaket_time1 = plaket_time(conf.array);

  stream << plaket1 << "," << plaket_space1 << "," << plaket_time1 << endl;

  end_time = omp_get_wtime();
  observables_time = end_time - start_time;
  cout << "plaket time: " << observables_time << endl;

  stream.close();
}
