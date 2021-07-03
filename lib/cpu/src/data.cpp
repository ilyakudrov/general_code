#include "../include/data.h"

#define PLACE_DATA                                                             \
  (t) * 4 * x_size *y_size *z_size + (z)*4 * x_size *y_size + (y)*4 * x_size + \
      (x)*4 + dir - 1

#define PLACE_QC2DSTAG                                                         \
  (dir - 1) * x_size *y_size *z_size *t_size * 4 +                             \
      (t)*x_size *y_size *z_size * 4 + (z)*x_size *y_size * 4 +                \
      (y)*x_size * 4 + (x)*4

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

// global variables for lattice size (should be removed)
extern int x_size;
extern int y_size;
extern int z_size;
extern int t_size;

template <class T> data<T>::data() {
  array.reserve(4 * x_size * y_size * z_size * t_size);
}

template <> void data<su2>::read_float(std::string &file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  std::ifstream stream(file_name);
  std::vector<float> v(data_size1 * 4);
  if (!stream.read((char *)&v[0], data_size1 * 4 * sizeof(float)))
    std::cout << "read_float<su2> error: " << file_name << std::endl;
  su2 A;
  for (int i = 0; i < data_size1; i++) {
    A.a0 = (FLOAT)v[i * 4];
    A.a1 = (FLOAT)v[i * 4 + 1];
    A.a2 = (FLOAT)v[i * 4 + 2];
    A.a3 = (FLOAT)v[i * 4 + 3];
    array.push_back(A);
  }
  stream.close();
}

std::vector<float> read_full_ml5(std::string &file_name, int conf_num) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  std::ifstream stream(file_name);
  std::vector<float> v(conf_num * data_size1 * 4);
  if (!stream.read((char *)&v[0], conf_num * data_size1 * 4 * sizeof(float)))
    std::cout << "read_full_ml5<su2> error: " << file_name << std::endl;
  return v;
  stream.close();
}

template <>
void data<su2>::read_float_ml5(const std::vector<float> &array_ml5,
                               int conf_num) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  su2 A;
  int i;
  SPACE_ITER_START
  for (int mu = 1; mu <= 4; mu++) {
    if (mu == 4) {
      i = 4 *
          (t * x_size * y_size * z_size + x + y * x_size + z * x_size * y_size);
    } else {
      i = 4 * (t * x_size * y_size * z_size + x + y * x_size +
               z * x_size * y_size) +
          mu;
    }
    A.a0 = (FLOAT)array_ml5[data_size1 * 4 * conf_num + i * 4];
    A.a3 = (FLOAT)array_ml5[data_size1 * 4 * conf_num + i * 4 + 1];
    A.a2 = (FLOAT)array_ml5[data_size1 * 4 * conf_num + i * 4 + 2];
    A.a1 = (FLOAT)array_ml5[data_size1 * 4 * conf_num + i * 4 + 3];
    array.push_back(A);
  }
  SPACE_ITER_END
}

template <> void data<abelian>::read_float(std::string &file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  std::ifstream stream(file_name);
  std::vector<float> v(data_size1 * 4);
  if (!stream.read((char *)&v[0], data_size1 * sizeof(float)))
    std::cout << "read_float<abelian> error: " << file_name << std::endl;
  for (int i = 0; i < data_size1; i++) {
    array.push_back(abelian(1, (FLOAT)v[i]));
  }
  stream.close();
}

template <> void data<su2>::read_float_fortran(std::string &file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  std::ifstream stream(file_name);
  std::vector<float> v(data_size1 * 4 + 2);
  if (!stream.read((char *)&v[0], (data_size1 * 4 + 2) * sizeof(float)))
    std::cout << "read_float_fortran<su2> error: " << file_name << std::endl;
  su2 A;
  for (int i = 0; i < data_size1; i++) {
    A.a0 = (FLOAT)v[i * 4 + 1];
    A.a1 = (FLOAT)v[i * 4 + 2];
    A.a2 = (FLOAT)v[i * 4 + 3];
    A.a3 = (FLOAT)v[i * 4 + 4];
    array.push_back(A);
  }
  stream.close();
}

template <> void data<abelian>::read_float_fortran(std::string &file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  std::ifstream stream(file_name);
  std::vector<float> v(data_size1 + 2);
  if (!stream.read((char *)&v[0], (data_size1 + 2) * sizeof(float)))
    std::cout << "read_float_fortran<abelian> error: " << file_name
              << std::endl;
  for (int i = 0; i < data_size1; i++) {
    array.push_back(abelian(1, (FLOAT)v[i + 1]));
  }
  stream.close();
}

template <>
void data<abelian>::read_float_fortran_convert_abelian(std::string &file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  std::ifstream stream(file_name);
  std::vector<float> v(data_size1 * 4 + 2);
  if (!stream.read((char *)&v[0], (data_size1 * 4 + 2) * sizeof(float)))
    std::cout << "read_float_fortran_convert_abelian<abelian> error: "
              << file_name << std::endl;
  for (int i = 0; i < data_size1; i++) {
    array.push_back(abelian(1, (FLOAT)atan2(v[i * 4 + 4], v[i * 4 + 1])));
  }
  stream.close();
}

template <>
void data<abelian>::read_float_convert_abelian(std::string &file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  std::ifstream stream(file_name);
  std::vector<float> v(data_size1 * 4);
  if (!stream.read((char *)&v[0], (data_size1 * 4 + 1) * sizeof(float)))
    std::cout << "read_float_convert_abelian<abelian> error: " << file_name
              << std::endl;
  for (int i = 0; i < data_size1; i++) {
    array.push_back(abelian(1, (FLOAT)atan2(v[i * 4 + 3], v[i * 4 + 0])));
  }
  stream.close();
}

template <> void data<su2>::read_double(std::string &file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  std::ifstream stream(file_name);
  std::vector<double> v(data_size1 * 4);
  if (!stream.read((char *)&v[0], data_size1 * 4 * sizeof(double)))
    std::cout << "read_double<su2> error: " << file_name << std::endl;
  su2 A;
  for (int i = 0; i < data_size1; i++) {
    A.a0 = (FLOAT)v[i * 4];
    A.a1 = (FLOAT)v[i * 4 + 1];
    A.a2 = (FLOAT)v[i * 4 + 2];
    A.a3 = (FLOAT)v[i * 4 + 3];
    array.push_back(A);
  }
  stream.close();
}

template <> void data<abelian>::read_double(std::string &file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  std::ifstream stream(file_name);
  std::vector<double> v(data_size1 * 4);
  if (!stream.read((char *)&v[0], data_size1 * sizeof(double)))
    std::cout << "read_double<abelian> error: " << file_name << std::endl;
  for (int i = 0; i < data_size1; i++) {
    array.push_back(abelian(1, (FLOAT)v[i]));
  }
  stream.close();
}

template <> void data<su2>::read_double_qc2dstag(std::string &file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  std::ifstream stream(file_name);
  std::vector<double> v(data_size1 * 4);
  if (!stream.read((char *)&v[0], data_size1 * 4 * sizeof(double)))
    std::cout << "read_double_qc2dstag<su2> error: " << file_name << std::endl;
  su2 A;
  int dir;
  SPACE_ITER_START
  for (int dir1 = 1; dir1 <= 4; dir1++) {
    if (dir1 == 4)
      dir = 1;
    else
      dir = dir1 + 1;
    A.a0 = v[PLACE_QC2DSTAG + 0];
    A.a3 = v[PLACE_QC2DSTAG + 1];
    A.a2 = v[PLACE_QC2DSTAG + 2];
    A.a1 = v[PLACE_QC2DSTAG + 3];
    array.push_back(A);
  }
  SPACE_ITER_END
  stream.close();
}

template <> void data<abelian>::read_double_qc2dstag(std::string &file_name) {
  std::cout << "there's no reason for implementation" << std::endl;
}

template <> void data<su2>::read_double_fortran(std::string &file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  std::ifstream stream(file_name);
  std::vector<double> v(data_size1 * 4);
  stream.ignore(4);
  if (!stream.read((char *)&v[0], (data_size1 * 4) * sizeof(double)))
    std::cout << "read_double_fortran<su2> error: " << file_name << std::endl;
  su2 A;
  for (int i = 0; i < data_size1; i++) {
    A.a0 = (FLOAT)v[i * 4];
    A.a1 = (FLOAT)v[i * 4 + 1];
    A.a2 = (FLOAT)v[i * 4 + 2];
    A.a3 = (FLOAT)v[i * 4 + 3];
    array.push_back(A);
  }
  stream.close();
}

template <> void data<abelian>::read_double_fortran(std::string &file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  std::ifstream stream(file_name);
  std::vector<double> v(data_size1);
  stream.ignore(4);
  if (!stream.read((char *)&v[0], (data_size1) * sizeof(double)))
    std::cout << "read_double_fortran<abelian> error: " << file_name
              << std::endl;
  for (int i = 0; i < data_size1; i++) {
    array.push_back(abelian(1, (FLOAT)v[i]));
  }
  stream.close();
}

template <> void data<su2>::write_float(std::string &file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  std::ofstream stream(file_name);
  std::vector<float> v(data_size1 * 4);
  for (int i = 0; i < data_size1; i++) {
    v[i * 4] = (float)array[i].a0;
    v[i * 4 + 1] = (float)array[i].a1;
    v[i * 4 + 2] = (float)array[i].a2;
    v[i * 4 + 3] = (float)array[i].a3;
  }
  if (!stream.write((char *)&v[0], data_size1 * 4 * sizeof(float)))
    std::cout << "write_float<su2> error: " << file_name << std::endl;
  stream.close();
}

template <> void data<abelian>::write_float(std::string &file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  std::ofstream stream(file_name);
  std::vector<float> v(data_size1);
  for (int i = 0; i < data_size1; i++) {
    v[i] = (float)array[i].phi;
  }
  if (!stream.write((char *)&v[0], (data_size1) * sizeof(float)))
    std::cout << "write_float<abelian> error: " << file_name << std::endl;
  stream.close();
}

template <> void data<abelian>::write_float_fortran(std::string &file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  std::ofstream stream(file_name);
  std::vector<float> v(data_size1 + 2);
  for (int i = 0; i < data_size1; i++) {
    v[i + 1] = (float)array[i].phi;
  }
  if (!stream.write((char *)&v[0], (data_size1 + 2) * sizeof(float)))
    std::cout << "write_float_fortran<abelian> error: " << file_name
              << std::endl;
  stream.close();
}

template <> void data<su2>::write_double(std::string &file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  std::ofstream stream(file_name);
  std::vector<FLOAT> v(data_size1 * 4);
  for (int i = 0; i < data_size1; i++) {
    v[i * 4] = (FLOAT)array[i].a0;
    v[i * 4 + 1] = (FLOAT)array[i].a1;
    v[i * 4 + 2] = (FLOAT)array[i].a2;
    v[i * 4 + 3] = (FLOAT)array[i].a3;
  }
  if (!stream.write((char *)&v[0], data_size1 * 4 * sizeof(double)))
    std::cout << "write_double<su2> error: " << file_name << std::endl;
  stream.close();
}

template <> void data<abelian>::write_double(std::string &file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  std::ofstream stream(file_name);
  std::vector<double> v(data_size1);
  for (int i = 0; i < data_size1; i++) {
    v[i] = (double)array[i].phi;
  }
  if (!stream.write((char *)&v[0], (data_size1) * sizeof(double)))
    std::cout << "write_double<abelian> error: " << file_name << std::endl;
  stream.close();
}

template class data<su2>;
template class data<abelian>;