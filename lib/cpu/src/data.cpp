#include "data.h"

extern int x_size;
extern int y_size;
extern int z_size;
extern int t_size;

using namespace std;

template <class T> data<T>::data() {
  array.reserve(4 * x_size * y_size * z_size * t_size);
}

template <> void data<su2>::read_float(char const *file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  ifstream stream(file_name);
  vector<float> v(data_size1 * 4);
  if (!stream.read((char *)&v[0], data_size1 * 4 * sizeof(float)))
    cout << "read_float error: " << file_name << endl;
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

template <> void data<abelian>::read_float(char const *file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  ifstream stream(file_name);
  vector<float> v(data_size1 * 4);
  if (!stream.read((char *)&v[0], data_size1 * sizeof(float)))
    cout << "read_float error: " << file_name << endl;
  for (int i = 0; i < data_size1; i++) {
    array.push_back(abelian(1, (FLOAT)v[i]));
  }
  stream.close();
}

template <> void data<su2>::read_float_fortran(char const *file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  ifstream stream(file_name);
  vector<float> v(data_size1 * 4 + 2);
  if (!stream.read((char *)&v[0], (data_size1 * 4 + 2) * sizeof(float)))
    cout << "read_float error: " << file_name << endl;
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

template <> void data<abelian>::read_float_fortran(char const *file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  ifstream stream(file_name);
  vector<float> v(data_size1 + 2);
  if (!stream.read((char *)&v[0], (data_size1 + 2) * sizeof(float)))
    cout << "read_float error: " << file_name << endl;
  for (int i = 0; i < data_size1; i++) {
    array.push_back(abelian(1, (FLOAT)v[i + 1]));
  }
  stream.close();
}

template <>
void data<abelian>::read_float_fortran_convert_abelian(char const *file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  ifstream stream(file_name);
  vector<float> v(data_size1 * 4 + 2);
  if (!stream.read((char *)&v[0], (data_size1 * 4 + 2) * sizeof(float)))
    cout << "read_float error: " << file_name << endl;
  for (int i = 0; i < data_size1; i++) {
    array.push_back(abelian(1, (FLOAT)atan2(v[i * 4 + 4], v[i * 4 + 1])));
  }
  stream.close();
}

template <>
void data<abelian>::read_float_convert_abelian(char const *file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  ifstream stream(file_name);
  vector<float> v(data_size1 * 4);
  if (!stream.read((char *)&v[0], (data_size1 * 4) * sizeof(float)))
    cout << "read_float error: " << file_name << endl;
  for (int i = 0; i < data_size1; i++) {
    array.push_back(abelian(1, (FLOAT)atan2(v[i * 4 + 3], v[i * 4 + 0])));
  }
  stream.close();
}

template <> void data<su2>::read_double(char const *file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  ifstream stream(file_name);
  vector<double> v(data_size1 * 4);
  if (!stream.read((char *)&v[0], data_size1 * 4 * sizeof(double)))
    cout << "read_double error: " << file_name << endl;
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

template <> void data<su2>::read_double_fortran(char const *file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  ifstream stream(file_name);
  vector<double> v(data_size1 * 4);
  stream.ignore(4);
  if (!stream.read((char *)&v[0], (data_size1 * 4) * sizeof(double)))
    cout << "read_float error: " << file_name << endl;
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

template <> void data<su2>::write_float(char const *file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  ofstream stream(file_name);
  vector<float> v(data_size1 * 4);
  for (int i = 0; i < data_size1; i++) {
    v[i * 4] = (float)array[i].a0;
    v[i * 4 + 1] = (float)array[i].a1;
    v[i * 4 + 2] = (float)array[i].a2;
    v[i * 4 + 3] = (float)array[i].a3;
  }
  if (!stream.write((char *)&v[0], data_size1 * 4 * sizeof(float)))
    cout << "write_float error: " << file_name << endl;
  stream.close();
}

template <> void data<abelian>::write_float_fortran(char const *file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  ofstream stream(file_name);
  vector<float> v(data_size1 + 2);
  for (int i = 0; i < data_size1; i++) {
    v[i + 1] = (float)array[i].phi;
  }
  if (!stream.write((char *)&v[0], (data_size1 + 2) * sizeof(float)))
    cout << "write_double error: " << file_name << endl;
  stream.close();
}

template <> void data<su2>::write_double(char const *file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  ofstream stream(file_name);
  vector<FLOAT> v(data_size1 * 4);
  for (int i = 0; i < data_size1; i++) {
    v[i * 4] = (FLOAT)array[i].a0;
    v[i * 4 + 1] = (FLOAT)array[i].a1;
    v[i * 4 + 2] = (FLOAT)array[i].a2;
    v[i * 4 + 3] = (FLOAT)array[i].a3;
  }
  if (!stream.write((char *)&v[0], data_size1 * 4 * sizeof(double)))
    cout << "write_double error: " << file_name << endl;
  stream.close();
}

/*
data_float::data_float() {
  array.reserve(4 * x_size * y_size * z_size * t_size);
}

void data_float::read_float(char const *file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  ifstream stream(file_name);
  vector<float> v(data_size1 * 4);
  if (!stream.read((char *)&v[0], data_size1 * sizeof(float)))
    cout << "read_float error: " << file_name << endl;
  for (int i = 0; i < data_size1; i++) {
    array.push_back((FLOAT)v[i]);
  }
  stream.close();
}

void data_float::read_float_fortran(char const *file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  ifstream stream(file_name);
  vector<float> v(data_size1 + 2);
  if (!stream.read((char *)&v[0], (data_size1 + 2) * sizeof(float)))
    cout << "read_float error: " << file_name << endl;
  for (int i = 0; i < data_size1; i++) {
    array.push_back((FLOAT)v[i + 1]);
  }
  stream.close();
}

void data_float::write_float_fortran(char const *file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  ofstream stream(file_name);
  vector<float> v(data_size1 + 2);
  for (int i = 0; i < data_size1; i++) {
    v[i + 1] = (float)array[i];
  }
  if (!stream.write((char *)&v[0], (data_size1 + 2) * sizeof(float)))
    cout << "write_double error: " << file_name << endl;
  stream.close();
}

void data_float::read_float_fortran_convert_abelian(char const *file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  ifstream stream(file_name);
  vector<float> v(data_size1 * 4 + 2);
  if (!stream.read((char *)&v[0], (data_size1 * 4 + 2) * sizeof(float)))
    cout << "read_float error: " << file_name << endl;
  for (int i = 0; i < data_size1; i++) {
    array.push_back((FLOAT)atan2(v[i * 4 + 4], v[i * 4 + 1]));
  }
  stream.close();
}

void data_float::read_float_convert_abelian(char const *file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  ifstream stream(file_name);
  vector<float> v(data_size1 * 4);
  if (!stream.read((char *)&v[0], (data_size1 * 4) * sizeof(float)))
    cout << "read_float error: " << file_name << endl;
  for (int i = 0; i < data_size1; i++) {
    array.push_back((FLOAT)atan2(v[i * 4 + 3], v[i * 4 + 0]));
  }
  stream.close();
}
*/

template class data<su2>;
template class data<abelian>;