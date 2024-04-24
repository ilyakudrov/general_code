#pragma once

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
  void read_float(std::string &file_name, int bytes_skip);

  // read conf file of doubles, usual ordering
  void read_double(std::string &file_name, int bytes_skip);

  void read_double_abelian(std::string &file_name, int bytes_skip);

  void read_double_offdiagonal(std::string &file_name, int bytes_skip);

  void read_double_vitaly(std::string &file_name, int bytes_skip);

  void read_double_vitaly_abelian(std::string &file_name, int bytes_skip);

  void read_ildg(std::string &file_name);

  // read conf file of floats, ml5 ordering
  // takes vector of floats, which is obtained by read_full_ml5 function, and
  // number of a configuration
  void read_float_ml5(const std::vector<float> &array_ml5, int conf_num);

  // read conf file of doubles, qc2dstag ordering
  void read_double_qc2dstag(std::string &file_name);

  // writes conf in file, usual ordering, double
  void write_double(std::string &file_name);

  // writes conf in file, usual ordering, float
  void write_float(std::string &file_name);

  // writes conf in file, fortran ordering, float
  void write_float_fortran(std::string &file_name);

  void read_float_convert_abelian(std::string &file_name, int bytes_skip);
  void read_double_convert_abelian(std::string &file_name, int bytes_skip);
  void read_double_qc2dstag_convert_abelian(std::string &file_name);
};
// read conf_num configurations from ml5 file and write them to vector in order
std::vector<float> read_full_ml5(std::string &file_name, int conf_num);

template <class T>
void read_file(data<T> &conf_data, std::string file_path,
               std::string file_format, int bytes_skip);

template <class T>
void get_data(data<T> &conf_data, std::string file_path,
              std::string file_format, int bytes_skip, bool convert);

template <class T>
std::vector<T> swap_directions(std::vector<T> &conf, int dir1, int dir2);