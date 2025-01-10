#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/link.h"
#include "../../lib/cpu/include/matrix.h"

#include <iostream>
#include <math.h>
#include <sstream>
#include <vector>

#define SPACE_ITER_START                                                       \
  for (int t = 0; t < t_size; t++) {                                           \
    for (int z = 0; z < z_size; z++) {                                         \
      for (int y = 0; y < y_size; y++) {                                       \
        for (int x = 0; x < x_size; x++) {
#define SPACE_ITER_END                                                         \
  }                                                                            \
  }                                                                            \
  }                                                                            \
  }

#define PLACE_QC2DSTAG                                                         \
  dir *x_size *y_size *z_size *t_size + (t) * x_size *y_size *z_size +         \
      (z) * x_size *y_size + (y) * x_size + x

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

template <class T> vector<su2> conf_transform_qc2dstag(vector<T> &conf);

template <> vector<su2> conf_transform_qc2dstag(vector<su2> &conf) {
  int lattice_size = 4 * x_size * y_size * z_size * t_size;
  vector<su2> conf_qc2dstag(lattice_size);
  link1 link(x_size, y_size, z_size, t_size);
  int dir;
  SPACE_ITER_START
  link.go_update(x, y, z, t);
  for (int dir1 = 0; dir1 < 4; dir1++) {
    if (dir1 == 3)
      dir = 0;
    else
      dir = dir1 + 1;
    conf_qc2dstag[PLACE_QC2DSTAG].a0 = conf[link.place + dir1].a0;
    conf_qc2dstag[PLACE_QC2DSTAG].a1 = conf[link.place + dir1].a3;
    conf_qc2dstag[PLACE_QC2DSTAG].a2 = conf[link.place + dir1].a2;
    conf_qc2dstag[PLACE_QC2DSTAG].a3 = conf[link.place + dir1].a1;
  }
  SPACE_ITER_END
  return conf_qc2dstag;
}

su2 make_su2_from_abelian(abelian &A) {
  su2 B;
  B.a0 = cos(A.phi);
  B.a1 = sin(A.phi);
  B.a2 = 0;
  B.a3 = 0;
  return B;
}

template <> vector<su2> conf_transform_qc2dstag(vector<abelian> &conf) {
  int lattice_size = 4 * x_size * y_size * z_size * t_size;
  vector<su2> conf_qc2dstag(lattice_size);
  link1 link(x_size, y_size, z_size, t_size);
  int dir;
  su2 A;
  SPACE_ITER_START
  link.go_update(x, y, z, t);
  for (int dir1 = 0; dir1 < 4; dir1++) {
    if (dir1 == 3)
      dir = 0;
    else
      dir = dir1 + 1;
    A = make_su2_from_abelian(conf[link.place + dir1]);
    conf_qc2dstag[PLACE_QC2DSTAG].a0 = A.a0;
    conf_qc2dstag[PLACE_QC2DSTAG].a1 = A.a3;
    conf_qc2dstag[PLACE_QC2DSTAG].a2 = A.a2;
    conf_qc2dstag[PLACE_QC2DSTAG].a3 = A.a1;
  }
  SPACE_ITER_END
  return conf_qc2dstag;
}

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double observables_time;

  string conf_format;
  string conf_path;
  string conf_path_output;
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
    } else if (string(argv[i]) == "-conf_path_output") {
      conf_path_output = argv[++i];
    } else if (string(argv[i]) == "-convert") {
      istringstream(string(argv[++i])) >> convert;
    } else if (string(argv[i]) == "-L_spat") {
      L_spat = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-L_time") {
      L_time = stoi(string(argv[++i]));
    }
  }

  x_size = L_spat;
  y_size = L_spat;
  z_size = L_spat;
  t_size = L_time;
  size1 = x_size * y_size;
  size2 = x_size * y_size * z_size;

  cout << "conf_format " << conf_format << endl;
  cout << "conf_path " << conf_path << endl;
  cout << "bytes_skip " << bytes_skip << endl;
  cout << "convert " << convert << endl;
  cout << "L_spat " << L_spat << endl;
  cout << "L_time " << L_time << endl;
  cout << endl;

  Data::data<MATRIX> conf;
  get_data(conf, conf_path, conf_format, bytes_skip, convert);
  Data::data<su2> conf_qc2dstag;
  conf_qc2dstag.array = conf_transform_qc2dstag(conf.array);
  conf_qc2dstag.write_double(conf_path_output);
}