#include "../include/data.h"
#include "../include/c-lime/lime_fixed_types.h"
#include "../include/c-lime/lime_reader.h"
#include "../include/decomposition.h"
#include "../include/link.h"

#include <cstring>
#include <fstream>

#define PLACE_DATA                                                             \
  (t) * 4 * x_size *y_size *z_size + (z) * 4 * x_size *y_size +                \
      (y) * 4 * x_size + (x) * 4 + dir - 1

#define PLACE_QC2DSTAG                                                         \
  (dir - 1) * x_size *y_size *z_size *t_size * 4 +                             \
      (t) * x_size *y_size *z_size * 4 + (z) * x_size *y_size * 4 +            \
      (y) * x_size * 4 + (x) * 4

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

template <class T> data<T>::data() {}

template <> void data<su2>::read_float(std::string &file_name, int bytes_skip) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  array.reserve(4 * x_size * y_size * z_size * t_size);
  std::ifstream stream(file_name);
  std::vector<float> v(data_size * 4);
  stream.ignore(bytes_skip);
  if (!stream.read((char *)&v[0], (data_size * 4) * sizeof(float)))
    std::cout << "read_float<su2> error: " << file_name << std::endl;
  su2 A;
  for (int i = 0; i < data_size; i++) {
    A.a0 = (double)v[i * 4];
    A.a1 = (double)v[i * 4 + 1];
    A.a2 = (double)v[i * 4 + 2];
    A.a3 = (double)v[i * 4 + 3];
    array.push_back(A);
  }
  stream.close();
}

template <>
void data<abelian>::read_float(std::string &file_name, int bytes_skip) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  array.reserve(4 * x_size * y_size * z_size * t_size);
  std::ifstream stream(file_name);
  std::vector<float> v(data_size);
  stream.ignore(bytes_skip);
  if (!stream.read((char *)&v[0], (data_size) * sizeof(float)))
    std::cout << "read_float<abelian> error: " << file_name << std::endl;
  for (int i = 0; i < data_size; i++) {
    array.push_back(abelian(1, (double)v[i]));
  }
  stream.close();
}

template <> void data<su3>::read_float(std::string &file_name, int bytes_skip) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  array.reserve(4 * x_size * y_size * z_size * t_size);
  std::ifstream stream(file_name);
  std::vector<float> v(data_size * 18);
  stream.ignore(bytes_skip);
  if (!stream.read((char *)&v[0], (data_size * 18) * sizeof(float)))
    std::cout << "read_float<su3> error: " << file_name << std::endl;
  su3 A;
  for (int i = 0; i < data_size; i++) {

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {

        A.matrix[i][j].real = v[i * 18 + (i * 3 + j) * 2];
        A.matrix[i][j].imag = v[i * 18 + (i * 3 + j) * 2 + 1];
      }
    }
    array.push_back(A);
  }
  stream.close();
}

template <>
void data<su3_abelian>::read_float(std::string &file_name, int bytes_skip) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  std::ifstream stream(file_name);
  std::vector<float> v(data_size * 4);
  array = std::vector<su3_abelian>(data_size);
  stream.ignore(bytes_skip);
  for (int i = 0; i < 3; i++) {
    if (!stream.read((char *)&v[0], (data_size) * sizeof(float)))
      std::cout << "read_float<su3_abelian> error: " << file_name << std::endl;
    for (int j = 0; j < data_size; j++) {
      array[j].matrix[i] = complex_t(cos(v[j]), sin(v[j]));
    }
  }
}

template <>
void data<su3_angles>::read_float(std::string &file_name, int bytes_skip) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  std::ifstream stream(file_name);
  std::vector<float> v(data_size * 4);
  array = std::vector<su3_angles>(data_size);
  stream.ignore(bytes_skip);
  for (int i = 0; i < 3; i++) {
    if (!stream.read((char *)&v[0], (data_size) * sizeof(float)))
      std::cout << "read_float<su3_angles> error: " << file_name << std::endl;
    for (int j = 0; j < data_size; j++) {
      array[j].matrix[i] = v[j];
    }
  }
}

template <>
void data<abelian>::read_double(std::string &file_name, int bytes_skip) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  array.reserve(4 * x_size * y_size * z_size * t_size);
  std::ifstream stream(file_name);
  std::vector<double> v(data_size);
  stream.ignore(bytes_skip);
  if (!stream.read((char *)&v[0], (data_size) * sizeof(double)))
    std::cout << "read_double<abelian> error: " << file_name << std::endl;
  for (int i = 0; i < data_size; i++) {
    array.push_back(abelian(1, (double)v[i]));
  }
  stream.close();
}

template <>
void data<su2>::read_double(std::string &file_name, int bytes_skip) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  array.reserve(data_size);
  std::ifstream stream(file_name);
  std::vector<double> v(data_size * 4);
  stream.ignore(bytes_skip);
  if (!stream.read((char *)&v[0], (data_size * 4) * sizeof(double)))
    std::cout << "read_double<su2> error: " << file_name << std::endl;
  su2 A;
  for (int i = 0; i < data_size; i++) {
    A.a0 = (double)v[i * 4];
    A.a1 = (double)v[i * 4 + 1];
    A.a2 = (double)v[i * 4 + 2];
    A.a3 = (double)v[i * 4 + 3];
    array.push_back(A);
  }
  stream.close();
}

template <>
void data<su3>::read_double(std::string &file_name, int bytes_skip) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  array = std::vector<su3>(data_size);
  std::ifstream stream(file_name);
  std::vector<double> v(data_size * 18);
  stream.ignore(bytes_skip);
  if (!stream.read((char *)&v[0], data_size * 18 * sizeof(double)))
    std::cout << "data<su3>::read_double error: " << file_name << std::endl;
  long int index = 0;
  link1 link(x_size, y_size, z_size, t_size);
  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
          for (int mu = 0; mu < 4; mu++) {
            for (int k = 0; k < 3; k++) {
              for (int j = 0; j < 3; j++) {
                link.go_update(x, y, z, t);
                array[link.place + mu].matrix[k][j] =
                    complex_t(v[index], v[index + 1]);
                index += 2;
              }
            }
          }
        }
      }
    }
  }
  stream.close();
}

template <>
void data<su3_abelian>::read_double(std::string &file_name, int bytes_skip) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  std::ifstream stream(file_name);
  std::vector<double> v(data_size);
  array = std::vector<su3_abelian>(data_size);
  stream.ignore(bytes_skip);
  for (int i = 0; i < 3; i++) {
    if (!stream.read((char *)&v[0], (data_size) * sizeof(double)))
      std::cout << "read_double<su3_abelian> error: " << file_name << std::endl;
    for (int j = 0; j < data_size; j++) {
      array[j].matrix[i] = complex_t(cos(v[j]), sin(v[j]));
    }
  }
}

template <>
void data<su3_angles>::read_double(std::string &file_name, int bytes_skip) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  std::ifstream stream(file_name);
  std::vector<double> v(data_size);
  array = std::vector<su3_angles>(data_size);
  stream.ignore(bytes_skip);
  for (int i = 0; i < 3; i++) {
    if (!stream.read((char *)&v[0], (data_size) * sizeof(double)))
      std::cout << "read_double<su3_angles> error: " << file_name << std::endl;
    for (int j = 0; j < data_size; j++) {
      array[j].matrix[i] = v[j];
    }
  }
}

template <>
void data<abelian>::read_double_abelian(std::string &file_name,
                                        int bytes_skip) {
  std::cout << "read_double_abelian<abelian> is not implemented" << std::endl;
}

template <>
void data<su2>::read_double_abelian(std::string &file_name, int bytes_skip) {
  std::cout << "read_double_abelian<su2> is not implemented" << std::endl;
}

template <>
void data<su3>::read_double_abelian(std::string &file_name, int bytes_skip) {
  std::cout << "read_double_abelian<su3> is not implemented" << std::endl;
}

template <>
void data<su3_abelian>::read_double_abelian(std::string &file_name,
                                            int bytes_skip) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  array = std::vector<su3_abelian>(data_size);
  std::ifstream stream(file_name);
  std::vector<double> v(data_size * 18);
  stream.ignore(bytes_skip);
  if (!stream.read((char *)&v[0], data_size * 18 * sizeof(double)))
    std::cout << "data<su3_abelian>::read_double_abelian error: " << file_name
              << std::endl;
  long int index = 0;
  complex_t complex_tmp;
  link1 link(x_size, y_size, z_size, t_size);
  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
          for (int mu = 0; mu < 4; mu++) {
            for (int k = 0; k < 3; k++) {
              for (int j = 0; j < 3; j++) {
                link.go_update(x, y, z, t);
                if (k == j) {
                  complex_tmp = complex_t(v[index], v[index + 1]);
                  array[link.place + mu].matrix[j] =
                      complex_tmp / complex_tmp.module();
                }
                index += 2;
              }
            }
          }
        }
      }
    }
  }
  stream.close();
}

template <>
void data<abelian>::read_double_offdiagonal(std::string &file_name,
                                            int bytes_skip) {
  std::cout << "read_double_offdiagonal<abelian> is not implemented"
            << std::endl;
}

template <>
void data<su2>::read_double_offdiagonal(std::string &file_name,
                                        int bytes_skip) {
  std::cout << "read_double_offdiagonal<su2> is not implemented" << std::endl;
}

template <>
void data<su3>::read_double_offdiagonal(std::string &file_name,
                                        int bytes_skip) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  array = std::vector<su3>(data_size);
  std::ifstream stream(file_name);
  std::vector<double> v(data_size * 18);
  stream.ignore(bytes_skip);
  if (!stream.read((char *)&v[0], data_size * 18 * sizeof(double)))
    std::cout << "data<su3>::read_double_offdiagonal error: " << file_name
              << std::endl;
  long int index = 0;
  complex_t complex_tmp;
  link1 link(x_size, y_size, z_size, t_size);
  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
          for (int mu = 0; mu < 4; mu++) {
            for (int k = 0; k < 3; k++) {
              for (int j = 0; j < 3; j++) {
                link.go_update(x, y, z, t);
                array[link.place + mu].matrix[k][j] =
                    complex_t(v[index], v[index + 1]);
                index += 2;
              }
            }
          }
        }
      }
    }
  }
  array = get_offdiagonal(array);
  stream.close();
}

template <>
void data<su3_abelian>::read_double_offdiagonal(std::string &file_name,
                                                int bytes_skip) {
  std::cout << "read_double_offdiagonal<su3_abelian> is not implemented"
            << std::endl;
}

template <>
void data<su3_abelian>::read_double_vitaly(std::string &file_name,
                                           int bytes_skip) {
  int lattice_size = x_size * y_size * z_size * t_size;
  std::ifstream stream(file_name);
  std::vector<double> v(lattice_size * 4 * 3);
  array = std::vector<su3_abelian>(lattice_size * 4);
  stream.ignore(bytes_skip);
  if (!stream.read((char *)&v[0], (lattice_size * 4 * 3) * sizeof(double)))
    std::cout << "read_double<su3_abelian> error: " << file_name << std::endl;
  for (int i = 0; i < lattice_size; i++) {
    for (int mu = 0; mu < 4; mu++) {
      for (int j = 0; j < 3; j++) {
        array[i * 4 + mu].matrix[j] =
            complex_t(cos(v[mu * lattice_size * 3 + j * lattice_size + i]),
                      sin(v[mu * lattice_size * 3 + j * lattice_size + i]));
      }
    }
  }
}

template <>
void data<su3_angles>::read_double_vitaly(std::string &file_name,
                                          int bytes_skip) {
  int lattice_size = x_size * y_size * z_size * t_size;
  std::ifstream stream(file_name);
  std::vector<double> v(lattice_size * 4 * 3);
  array = std::vector<su3_angles>(lattice_size * 4);
  stream.ignore(bytes_skip);
  if (!stream.read((char *)&v[0], (lattice_size * 4 * 3) * sizeof(double)))
    std::cout << "read_double<su3_angles> error: " << file_name << std::endl;
  for (int i = 0; i < lattice_size; i++) {
    for (int mu = 0; mu < 4; mu++) {
      for (int j = 0; j < 3; j++) {
        array[i * 4 + mu].matrix[j] =
            v[mu * lattice_size * 3 + j * lattice_size + i];
      }
    }
  }
}

template <>
void data<su2>::read_double_vitaly(std::string &file_name, int bytes_skip) {
  std::cout << "read_double_vitaly<su2> is not implemented" << std::endl;
}

template <>
void data<su3>::read_double_vitaly(std::string &file_name, int bytes_skip) {
  int lattice_size = x_size * y_size * z_size * t_size;
  std::ifstream stream(file_name);
  std::vector<double> v(lattice_size * 4 * 18);
  array = std::vector<su3>(lattice_size * 4);
  stream.ignore(bytes_skip);
  if (!stream.read((char *)&v[0], (lattice_size * 4 * 18) * sizeof(double)))
    std::cout << "data<su3>::read_double_vitaly error: " << file_name
              << std::endl;
  for (int i = 0; i < lattice_size; i++) {
    for (int mu = 0; mu < 4; mu++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          array[i * 4 + mu].matrix[j][k] =
              complex_t(v[mu * lattice_size * 18 + k * lattice_size * 6 +
                          j * lattice_size * 2 + i * 2],
                        v[mu * lattice_size * 18 + k * lattice_size * 6 +
                          j * lattice_size * 2 + i * 2 + 1]);
        }
      }
    }
  }
}

template <>
void data<abelian>::read_double_vitaly(std::string &file_name, int bytes_skip) {
  std::cout << "read_double_vitaly<abelian> is not implemented" << std::endl;
}

template <>
void data<su3_abelian>::read_double_vitaly_abelian(std::string &file_name,
                                                   int bytes_skip) {
  int lattice_size = x_size * y_size * z_size * t_size;
  std::ifstream stream(file_name);
  std::vector<double> v(lattice_size * 4 * 18);
  array = std::vector<su3_abelian>(lattice_size * 4);
  stream.ignore(bytes_skip);
  if (!stream.read((char *)&v[0], (lattice_size * 4 * 18) * sizeof(double)))
    std::cout << "data<su3>::read_double_vitaly error: " << file_name
              << std::endl;
  complex_t complex_tmp;
  for (int i = 0; i < lattice_size; i++) {
    for (int mu = 0; mu < 4; mu++) {
      for (int j = 0; j < 3; j++) {
        complex_tmp =
            complex_t(v[mu * lattice_size * 18 + j * lattice_size * 6 +
                        j * lattice_size * 2 + i * 2],
                      v[mu * lattice_size * 18 + j * lattice_size * 6 +
                        j * lattice_size * 2 + i * 2 + 1]);
        array[i * 4 + mu].matrix[j] = complex_tmp / complex_tmp.module();
      }
    }
  }
}

template <>
void data<su2>::read_double_vitaly_abelian(std::string &file_name,
                                           int bytes_skip) {
  std::cout << "read_double_vitaly_abelian<su2> is not implemented"
            << std::endl;
}

template <>
void data<su3>::read_double_vitaly_abelian(std::string &file_name,
                                           int bytes_skip) {
  std::cout << "read_double_vitaly_abelian<su3> is not implemented"
            << std::endl;
}

template <>
void data<abelian>::read_double_vitaly_abelian(std::string &file_name,
                                               int bytes_skip) {
  std::cout << "read_double_vitaly_abelian<abelian> is not implemented"
            << std::endl;
}

double reverseValue(const char *data) {
  double result;

  char *dest = (char *)&result;

  for (int i = 0; i < sizeof(double); i++) {
    dest[i] = data[sizeof(double) - i - 1];
  }
  return result;
}

template <> void data<su3>::read_ildg(std::string &file_name) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  array.reserve(4 * x_size * y_size * z_size * t_size);

  FILE *fp;
  fp = fopen(file_name.c_str(), "r");

  LimeReader *reader;
  reader = limeCreateReader(fp);

  int status;
  char *lime_type;
  n_uint64_t nbytes;
  while ((status = limeReaderNextRecord(reader)) != LIME_EOF) {

    if (status != LIME_SUCCESS) {
      fprintf(stderr, "limeReaderNextRecord returned status = %d\n", status);
    }

    lime_type = limeReaderType(reader);
    nbytes = limeReaderBytes(reader);

    if (strcmp(lime_type, "ildg-binary-data") == 0) {

      std::vector<double> v(nbytes / sizeof(double));

      status = limeReaderReadData((char *)&v[0], &nbytes, reader);

      for (int i = 0; i < nbytes / sizeof(double); i++) {
        v[i] = reverseValue((char *)&v[i]);
      }

      if (status != LIME_SUCCESS) {
        fprintf(stderr, "limeReaderReadData returned status = %d\n", status);
      }

      su3 A;
      int place;
      for (int i = 0; i < data_size1; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            place = i * 18 + j * 6 + k * 2;
            A.matrix[j][k] = complex_t(v[place], v[place + 1]);
          }
        }
        array.push_back(A);
      }
    }
  }
}

template <> void data<su2>::read_ildg(std::string &file_name) {
  std::cout << "read_ildg<su2> is not implemented" << std::endl;
}

template <> void data<abelian>::read_ildg(std::string &file_name) {
  std::cout << "read_ildg<abelian> is not implemented" << std::endl;
}

template <> void data<su3_abelian>::read_ildg(std::string &file_name) {
  std::cout << "read_ildg<su3_abelian> is not implemented" << std::endl;
}

template <> void data<su3_angles>::read_ildg(std::string &file_name) {
  std::cout << "read_ildg<su3_angles> is not implemented" << std::endl;
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
  array.reserve(4 * x_size * y_size * z_size * t_size);
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
    A.a0 = (double)array_ml5[data_size1 * 4 * conf_num + i * 4];
    A.a3 = (double)array_ml5[data_size1 * 4 * conf_num + i * 4 + 1];
    A.a2 = (double)array_ml5[data_size1 * 4 * conf_num + i * 4 + 2];
    A.a1 = (double)array_ml5[data_size1 * 4 * conf_num + i * 4 + 3];
    array.push_back(A);
  }
  SPACE_ITER_END
}

template <>
void data<abelian>::read_float_convert_abelian(std::string &file_name,
                                               int bytes_skip) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  array.reserve(4 * x_size * y_size * z_size * t_size);
  std::ifstream stream(file_name);
  std::vector<float> v(data_size * 4);
  stream.ignore(bytes_skip);
  if (!stream.read((char *)&v[0], (data_size * 4) * sizeof(float)))
    std::cout << "read_float_convert_abelian<abelian> error: " << file_name
              << std::endl;
  for (int i = 0; i < data_size; i++) {
    array.push_back(
        abelian(1, atan2((double)v[i * 4 + 3], (double)v[i * 4 + 0])));
  }
  stream.close();
}

template <>
void data<su2>::read_float_convert_abelian(std::string &file_name,
                                           int bytes_skip) {
  std::cout << "wrong type for read_float_convert_abelian" << std::endl;
}

template <>
void data<abelian>::read_double_convert_abelian(std::string &file_name,
                                                int bytes_skip) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  array.reserve(data_size);
  std::ifstream stream(file_name);
  std::vector<double> v(data_size * 4);
  stream.ignore(bytes_skip);
  if (!stream.read((char *)&v[0], (data_size * 4) * sizeof(double)))
    std::cout << "read_double_convert_abelian<abelian> error: " << file_name
              << std::endl;
  for (int i = 0; i < data_size; i++) {
    array.push_back(abelian(1, atan2(v[i * 4 + 3], v[i * 4 + 0])));
  }
  stream.close();
}

template <>
void data<su2>::read_double_convert_abelian(std::string &file_name,
                                            int bytes_skip) {
  std::cout << "wrong type for read_double_convert_abelian" << std::endl;
}

template <>
void data<abelian>::read_double_qc2dstag_convert_abelian(
    std::string &file_name) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  array.reserve(data_size);
  std::ifstream stream(file_name);
  std::vector<double> v(data_size * 4);
  if (!stream.read((char *)&v[0], data_size * 4 * sizeof(double)))
    std::cout << "read_double_qc2dstag<su2> error: " << file_name << std::endl;
  int dir;
  SPACE_ITER_START
  for (int dir1 = 1; dir1 <= 4; dir1++) {
    if (dir1 == 4)
      dir = 1;
    else
      dir = dir1 + 1;
    array.push_back(abelian(
        1, (double)atan2(v[PLACE_QC2DSTAG + 3], v[PLACE_QC2DSTAG + 0])));
  }
  SPACE_ITER_END
  stream.close();
}

template <>
void data<su2>::read_double_qc2dstag_convert_abelian(std::string &file_name) {
  std::cout << "wrong type for read_double_qc2dstag_convert_abelian"
            << std::endl;
}

template <> void data<su2>::read_double_qc2dstag(std::string &file_name) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  array.reserve(data_size);
  std::ifstream stream(file_name);
  std::vector<double> v(data_size * 4);
  if (!stream.read((char *)&v[0], data_size * 4 * sizeof(double)))
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

template <> void data<su3>::read_double_qc2dstag(std::string &file_name) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  array.clear();
  array.reserve(data_size);
  std::ifstream stream(file_name);
  std::vector<double> v(data_size * 18);
  if (!stream.read((char *)&v[0], data_size * 18 * sizeof(double)))
    std::cout << "read_double_qc2dstag<su3> error: " << file_name << std::endl;
  su3 A;
  int dir;
  int place;
  SPACE_ITER_START
  for (int dir1 = 1; dir1 <= 4; dir1++) {
    if (dir1 == 4)
      dir = 1;
    else
      dir = dir1 + 1;

    int place =
        ((dir - 1) * x_size * y_size * z_size * t_size +
         t * x_size * y_size * z_size + z * x_size * y_size + y * x_size + x) *
        18;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {

        A.matrix[i][j].real = v[place + (i * 3 + j) * 2];
        A.matrix[i][j].imag = v[place + (i * 3 + j) * 2 + 1];
      }
    }
    array.push_back(A);
  }
  SPACE_ITER_END
  stream.close();
}

template <> void data<abelian>::read_double_qc2dstag(std::string &file_name) {
  std::cout << "there's no implementation for abelian" << std::endl;
}

template <>
void data<su3_abelian>::read_double_qc2dstag(std::string &file_name) {
  std::cout << "there's no implementation for su3_abelian" << std::endl;
}

template <>
void data<su3_angles>::read_double_qc2dstag(std::string &file_name) {
  std::cout << "there's no implementation for su3_angles" << std::endl;
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
  std::vector<double> v(data_size1 * 4);
  for (int i = 0; i < data_size1; i++) {
    v[i * 4] = (double)array[i].a0;
    v[i * 4 + 1] = (double)array[i].a1;
    v[i * 4 + 2] = (double)array[i].a2;
    v[i * 4 + 3] = (double)array[i].a3;
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

template <> void data<su3>::write_double(std::string &file_name) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  std::ofstream stream(file_name);
  if (!stream.write((char *)&array[0], data_size * 18 * sizeof(double)))
    std::cout << "write_double<su3> error: " << file_name << std::endl;
  stream.close();
}

template <> void data<su3_abelian>::write_double(std::string &file_name) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  std::ofstream stream(file_name);

  for (int i = 0; i < 3; i++) {
    std::vector<double> angles(data_size);
    for (int j = 0; j < data_size; j++) {
      angles[j] = atan2(array[j].matrix[i].imag, array[j].matrix[i].real);
    }

    if (!stream.write((char *)&angles[0], data_size * sizeof(double)))
      std::cout << "write_double<su3_abelian> error: " << file_name
                << std::endl;
  }

  stream.close();
}

template class data<su2>;
template class data<abelian>;
template class data<su3>;
template class data<su3_abelian>;
template class data<su3_angles>;

template <class T>
void read_file(data<T> &conf_data, std::string file_path,
               std::string file_format, int bytes_skip) {
  if (std::string(file_format) == "float") {
    conf_data.read_float(file_path, bytes_skip);
  } else if (std::string(file_format) == "double") {
    conf_data.read_double(file_path, bytes_skip);
  } else if (std::string(file_format) == "double_vitaly") {
    conf_data.read_double_vitaly(file_path, bytes_skip);
  } else if (std::string(file_format) == "double_qc2dstag") {
    conf_data.read_double_qc2dstag(file_path);
  } else if (std::string(file_format) == "ildg") {
    conf_data.read_ildg(file_path);
  } else {
    std::cout << "wrong conf format: " << file_format << std::endl;
  }
}

template <>
void get_data(data<su3_abelian> &conf_data, std::string file_path,
              std::string file_format, int bytes_skip, bool convert) {
  if (!convert) {
    read_file(conf_data, file_path, file_format, bytes_skip);
  } else {
    data<su3> data_tmp;
    read_file(data_tmp, file_path, file_format, bytes_skip);
    conf_data.array = get_abelian(data_tmp.array);
  }
}

template <>
void get_data(data<su3_angles> &conf_data, std::string file_path,
              std::string file_format, int bytes_skip, bool convert) {
  if (!convert) {
    read_file(conf_data, file_path, file_format, bytes_skip);
  } else {
    data<su3> data_tmp;
    read_file(data_tmp, file_path, file_format, bytes_skip);
    conf_data.array = su3_to_su3_angles(data_tmp.array);
  }
}

template <>
void get_data(data<su3> &conf_data, std::string file_path,
              std::string file_format, int bytes_skip, bool convert) {
  if (!convert) {
    read_file(conf_data, file_path, file_format, bytes_skip);
  } else {
    data<su3> data_tmp;
    read_file(data_tmp, file_path, file_format, bytes_skip);
    conf_data.array = get_offdiagonal(data_tmp.array);
  }
}

template <>
void get_data(data<abelian> &conf_data, std::string file_path,
              std::string file_format, int bytes_skip, bool convert) {
  if (!convert) {
    read_file(conf_data, file_path, file_format, bytes_skip);
  } else {
    data<su2> data_tmp;
    read_file(data_tmp, file_path, file_format, bytes_skip);
    conf_data.array = get_abelian(data_tmp.array);
  }
}

template <>
void get_data(data<su2> &conf_data, std::string file_path,
              std::string file_format, int bytes_skip, bool convert) {
  if (!convert) {
    read_file(conf_data, file_path, file_format, bytes_skip);
  } else {
    data<su2> data_tmp;
    read_file(data_tmp, file_path, file_format, bytes_skip);
    conf_data.array = get_offdiagonal(data_tmp.array);
  }
}

template <class T>
std::vector<T> swap_directions(const std::vector<T> &conf, int dir1, int dir2) {
  std::vector<T> conf1 = conf;
  link1 link2(x_size, y_size, z_size, t_size);
  link1 link3(x_size, y_size, z_size, t_size);
  std::vector<int> coordinates(4);
  std::vector<int> coordinates1(4);
  int tmp;
  for (coordinates[3] = 0; coordinates[3] < t_size; coordinates[3]++) {
    for (coordinates[2] = 0; coordinates[2] < z_size; coordinates[2]++) {
      for (coordinates[1] = 0; coordinates[1] < y_size; coordinates[1]++) {
        for (coordinates[0] = 0; coordinates[0] < x_size; coordinates[0]++) {
          link2.go_update(coordinates);
          coordinates1 = coordinates;
          tmp = coordinates1[dir1];
          coordinates1[dir1] = coordinates1[dir2];
          coordinates1[dir2] = tmp;
          link3.go_update(coordinates1);
          conf1[link3.place + dir1] = conf[link2.place + dir2];
          conf1[link3.place + dir2] = conf[link2.place + dir1];
          for (int mu = 0; mu < 4; mu++) {
            if (mu != dir1 && mu != dir2) {
              conf1[link3.place + mu] = conf[link2.place + mu];
            }
          }
        }
      }
    }
  }
  return conf1;
}

// su2
template std::vector<su2> swap_directions(const std::vector<su2> &conf,
                                          int dir1, int dir2);

// abelian
template std::vector<abelian> swap_directions(const std::vector<abelian> &conf,
                                              int dir1, int dir2);

// su3
template std::vector<su3> swap_directions(const std::vector<su3> &conf,
                                          int dir1, int dir2);

// su3_abelian
template std::vector<su3_abelian>
swap_directions(const std::vector<su3_abelian> &conf, int dir1, int dir2);