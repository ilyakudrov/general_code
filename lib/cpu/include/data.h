#pragma once

#ifdef DOUBLE
#define FLOAT double
#else
#define FLOAT float
#endif

#include "math.h"
#include "matrix.h"
#include <fstream>
#include <iostream>
#include <vector>

// Class which contains array of matrices of type T
// Order of matrices used in program is (usual ordering)
// D * (t * Nx*Ny*Nz + z * Nx*Ny + y * Nx + x) + mu,
// where Nx, Ny, Nz, Nt - lattice size
// x, y, z, t - lattice coordinates
// D - dimension of the lattice (4)
// mu - link direction
template <class T> class data {
public:
  // array of matrices
  std::vector<T> array;
  data();

  // read conf file of floats, usual ordering
  void read_float(char const *file_name);

  // read conf file of floats, ml5 ordering
  // takes vector of floats, which is obtained by read_full_ml5 function, and number of a configuration
  void read_float_ml5(const vector<float> &array_ml5, int conf_num);

  // read conf file of floats, fortran ordering
  void read_float_fortran(char const *file_name);

  // read conf file of doubles, usual ordering
  void read_double(char const *file_name);

  // read conf file of doubles, qc2dstag ordering
  void read_double_qc2dstag(char const *file_name);

  // read conf file of doubles, fortran ordering
  void read_double_fortran(char const *file_name);

  // writes conf in file, usual ordering, double
  void write_double(char const *file_name);

  // writes conf in file, usual ordering, float
  void write_float(char const *file_name);

  // writes conf in file, fortran ordering, float
  void write_float_fortran(char const *file_name);

  void read_float_fortran_convert_abelian(char const *file_name);
  void read_float_convert_abelian(char const *file_name);
};
// read conf_num configurations from ml5 file and write them to vector in order
vector<float> read_full_ml5(char const *file_name, int conf_num);