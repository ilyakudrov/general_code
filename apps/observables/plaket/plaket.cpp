#include "../../../lib/cpu/include/plaket.h"
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
int size1;
int size2;

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double observables_time;

  string conf_format;
  string file_precision;
  string conf_path;
  string path;
  int L_spat, L_time;
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
    } else if (string(argv[i]) == "--path") {
      path = argv[++i];
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
  cout << "path " << path << endl;
  cout << endl;

  cout.precision(17);

  Data::LatticeData<DataPatternLexicographical, MATRIX> conf(
      {x_size1, y_size1, z_size1, t_size1});
  Data::read_data_convert(conf, conf_path, conf_format, bytes_skip,
                          file_precision, convert);

  start_time = omp_get_wtime();
  double plaket1 = plaket(conf);
  double plaket_space1 = plaket_space(conf);
  double plaket_time1 = plaket_time(conf);

  end_time = omp_get_wtime();
  observables_time = end_time - start_time;
  cout << "plaket time: " << observables_time << endl;

  ofstream stream;
  stream.precision(17);
  stream.open(path);
  stream << "plaket,plaket_space,plaket_time" << endl;
  stream << plaket1 << "," << plaket_space1 << "," << plaket_time1 << endl;
  stream.close();
}
