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

template <class T> class data {
public:
  std::vector<T> array;
  data();
  void read_float(char const *file_name); // read conf file of floats
  void read_float_fortran(char const *file_name);
  void read_double(char const *file_name); // read conf file of double
  void read_double_qc2dstag(char const *file_name);
  void read_double_fortran(char const *file_name);
  void write_double(char const *file_name); // writes in file
  void write_float(char const *file_name);

  void write_float_fortran(char const *file_name);
  void read_float_fortran_convert_abelian(char const *file_name);
  void read_float_convert_abelian(char const *file_name);
};