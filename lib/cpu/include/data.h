#pragma once

#include "math.h"

#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

// Class which contains array of matrices of type T
// Order of matrices used in program is (usual ordering)
// D * (t * Nx*Ny*Nz + z * Nx*Ny + y * Nx + x) + mu,
// where Nx, Ny, Nz, Nt - lattice size
// x, y, z, t - lattice coordinates
// D - dimension of the lattice (4)
// mu - link direction

namespace Data {
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

template <class DataPattern, class MatrixType> class Data1 {
public:
  std::vector<MatrixType> array;
  DataPattern data_pattern;
  Data1(const DataPattern &data_pattern) : data_pattern(data_pattern) {
    array = std::vector<MatrixType>(data_pattern.get_data_size());
  }

  template <class FilePattern, class precision>
  void fill_array(std::vector<precision> &data, FilePattern file_pattern) {
    MatrixType A;
    std::vector<double> matrix_data(MatrixType::data_size);
    for (auto lat_coord : data_pattern.get_multi_index()) {
      for (int mu = 0; mu < data_pattern.lat_dim.size(); mu++) {
        for (int element_num = 0; element_num < MatrixType::data_size;
             element_num++) {
          matrix_data[element_num] = data[file_pattern.get_index_matrix_data(
              data_pattern.lat_dim, lat_coord, mu, element_num,
              MatrixType::data_size)];
        }
        data_pattern.lat_coord = lat_coord;
        array[data_pattern.get_index_link(mu)] = MatrixType(matrix_data);
      }
    }
  }

  template <class FilePattern>
  void read_data(std::string file_path, FilePattern file_pattern,
                 int bytes_skip, std::string file_precision) {
    std::ifstream stream(file_path);
    std::vector<float> float_v;
    std::vector<double> double_v;
    if (file_precision == std::string_view("float")) {
      float_v.resize(data_pattern.get_data_size() * MatrixType::data_size);
      if (!stream.read((char *)&float_v[0],
                       (data_pattern.get_data_size() * MatrixType::data_size) *
                           sizeof(float))) {
        std::cout << "read_data float error: " << file_path << std::endl;
      }
      fill_array(float_v, file_pattern);
    } else if (file_precision == std::string_view("double")) {
      double_v.resize(data_pattern.get_data_size() * MatrixType::data_size);
      if (!stream.read((char *)&double_v[0],
                       (data_pattern.get_data_size() * MatrixType::data_size) *
                           sizeof(double))) {
        std::cout << "read_data double error: " << file_path << std::endl;
      }
      fill_array(double_v, file_pattern);
    } else {
      std::cout << "wrong precision in read_data" << std::endl;
    }
    stream.close();
  }

  template <class FilePattern>
  void fill_data_to_write(std::vector<double> &data, FilePattern file_pattern) {
    std::vector<double> matrix_data(MatrixType::data_size);
    for (auto lat_coord : data_pattern.get_multi_index()) {
      for (int mu = 0; mu < data_pattern.lat_dim.size(); mu++) {
        data_pattern.lat_coord = lat_coord;
        matrix_data = array[data_pattern.get_index_link(mu)].get_data();
        for (int element_num = 0; element_num < MatrixType::data_size;
             element_num++) {
          data[file_pattern.get_index_matrix_data(
              data_pattern.lat_dim, lat_coord, mu, element_num,
              MatrixType::data_size)] = matrix_data[element_num];
        }
      }
    }
  }

  template <class FilePattern>
  void write_data(std::string file_path, FilePattern file_pattern) {
    std::vector<double> data(data_pattern.get_data_size() *
                             MatrixType::data_size);
    fill_data_to_write(data, file_pattern);
    std::ofstream stream(file_path);
    if (!stream.write((char *)&data[0], data_pattern.get_data_size() *
                                            MatrixType::data_size *
                                            sizeof(double)))
      std::cout << "write_data error: " << file_path << std::endl;
    stream.close();
  }
};

} // namespace Data

// read conf_num configurations from ml5 file and write them to vector in order
std::vector<float> read_full_ml5(std::string &file_name, int conf_num);

template <class T>
void read_file(Data::data<T> &conf_data, std::string file_path,
               std::string file_format, int bytes_skip);

template <class T>
void get_data(Data::data<T> &conf_data, std::string file_path,
              std::string file_format, int bytes_skip, bool convert);

// only for symmetric lattice
template <class T>
std::vector<T> swap_directions(const std::vector<T> &conf, int dir1, int dir2);