#include "../include/monopoles.h"
#include "../include/data.h"
#include "../include/indexing.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <link.h>
#include <numeric>

#define SPACE_ITER_START                                                       \
  for (int t = 0; t < t_size; t++) {                                           \
    for (int z = 0; z < z_size; z++) {                                         \
      for (int y = 0; y < y_size; y++) {                                       \
        for (int x = 0; x < x_size; x++) {                                     \
          link.go(x, y, z, t);                                                 \
          link.update(0);                                                      \
          link.update(1);                                                      \
          link.update(2);                                                      \
          link.update(3);

#define SPACE_ITER_END                                                         \
  }                                                                            \
  }                                                                            \
  }                                                                            \
  }

#define SPACE_ITER_START_SIMPLE                                                \
  for (int t = 0; t < t_size; t++) {                                           \
    for (int z = 0; z < z_size; z++) {                                         \
      for (int y = 0; y < y_size; y++) {                                       \
        for (int x = 0; x < x_size; x++) {
#define SPACE_ITER_END_SIMPLE                                                  \
  }                                                                            \
  }                                                                            \
  }                                                                            \
  }

#define PLACE_QC2DSTAG                                                         \
  (dir - 1) * x_size *y_size *z_size *t_size * 4 +                             \
      (t) * x_size *y_size *z_size * 4 + (z) * x_size *y_size * 4 +            \
      (y) * x_size * 4 + (x) * 4

void get_current(std::vector<std::vector<double>> &monopole_plaket, double *J,
                 DataPatternLexicographical &data_pattern) {
  double j0, j1, j2, j3;
  DataPatternLexicographical data_pattern_x(data_pattern);
  data_pattern_x.move_forward(1, 0);
  DataPatternLexicographical data_pattern_y(data_pattern);
  data_pattern_y.move_forward(1, 1);
  DataPatternLexicographical data_pattern_z(data_pattern);
  data_pattern_z.move_forward(1, 2);
  DataPatternLexicographical data_pattern_t(data_pattern);
  data_pattern_t.move_forward(1, 3);

  int index;
  data_pattern.move_forward(1, 3);
  data_pattern_x.move_forward(1, 3);
  data_pattern_y.move_forward(1, 3);
  data_pattern_z.move_forward(1, 3);
  index = data_pattern.get_index_site();
  j3 = monopole_plaket[3][data_pattern_x.get_index_site()] -
       monopole_plaket[3][index] -
       (monopole_plaket[1][data_pattern_y.get_index_site()] -
        monopole_plaket[1][index]) +
       monopole_plaket[0][data_pattern_z.get_index_site()] -
       monopole_plaket[0][index];
  data_pattern.move_backward(1, 3);
  data_pattern_x.move_backward(1, 3);
  data_pattern_y.move_backward(1, 3);
  data_pattern_z.move_backward(1, 3);
  data_pattern.move_forward(1, 0);
  data_pattern_t.move_forward(1, 0);
  data_pattern_y.move_forward(1, 0);
  data_pattern_z.move_forward(1, 0);
  index = data_pattern.get_index_site();
  j0 = -(monopole_plaket[3][data_pattern_t.get_index_site()] -
         monopole_plaket[3][index]) -
       (monopole_plaket[5][data_pattern_y.get_index_site()] -
        monopole_plaket[5][index]) +
       monopole_plaket[4][data_pattern_z.get_index_site()] -
       monopole_plaket[4][index];
  data_pattern.move_backward(1, 0);
  data_pattern_t.move_backward(1, 0);
  data_pattern_y.move_backward(1, 0);
  data_pattern_z.move_backward(1, 0);
  data_pattern.move_forward(1, 1);
  data_pattern_t.move_forward(1, 1);
  data_pattern_x.move_forward(1, 1);
  data_pattern_z.move_forward(1, 1);
  index = data_pattern.get_index_site();
  j1 = monopole_plaket[1][data_pattern_t.get_index_site()] -
       monopole_plaket[1][index] +
       monopole_plaket[5][data_pattern_x.get_index_site()] -
       monopole_plaket[5][index] -
       (monopole_plaket[2][data_pattern_z.get_index_site()] -
        monopole_plaket[2][index]);
  data_pattern.move_backward(1, 1);
  data_pattern_t.move_backward(1, 1);
  data_pattern_x.move_backward(1, 1);
  data_pattern_z.move_backward(1, 1);
  data_pattern.move_forward(1, 2);
  data_pattern_t.move_forward(1, 2);
  data_pattern_x.move_forward(1, 2);
  data_pattern_y.move_forward(1, 2);
  index = data_pattern.get_index_site();
  j2 = -(monopole_plaket[0][data_pattern_t.get_index_site()] -
         monopole_plaket[0][index]) -
       (monopole_plaket[4][data_pattern_x.get_index_site()] -
        monopole_plaket[4][index]) +
       monopole_plaket[2][data_pattern_y.get_index_site()] -
       monopole_plaket[2][index];
  data_pattern.move_backward(1, 2);
  data_pattern_t.move_backward(1, 2);
  data_pattern_x.move_backward(1, 2);
  data_pattern_y.move_backward(1, 2);

  J[0] = j0 / 2 / M_PI;
  J[1] = j1 / 2 / M_PI;
  J[2] = j2 / 2 / M_PI;
  J[3] = j3 / 2 / M_PI;
}

void get_current_singular(std::vector<std::vector<int>> &monopole_plaket,
                          double *J, DataPatternLexicographical &data_pattern) {
  double j0, j1, j2, j3;
  DataPatternLexicographical data_pattern_x(data_pattern);
  data_pattern_x.move_forward(1, 0);
  DataPatternLexicographical data_pattern_y(data_pattern);
  data_pattern_y.move_forward(1, 1);
  DataPatternLexicographical data_pattern_z(data_pattern);
  data_pattern_z.move_forward(1, 2);
  DataPatternLexicographical data_pattern_t(data_pattern);
  data_pattern_t.move_forward(1, 3);

  int index;
  data_pattern.move_forward(1, 3);
  data_pattern_x.move_forward(1, 3);
  data_pattern_y.move_forward(1, 3);
  data_pattern_z.move_forward(1, 3);
  index = data_pattern.get_index_site();
  j3 = -(monopole_plaket[3][data_pattern_x.get_index_site()] -
         monopole_plaket[3][index] -
         (monopole_plaket[1][data_pattern_y.get_index_site()] -
          monopole_plaket[1][index]) +
         monopole_plaket[0][data_pattern_z.get_index_site()] -
         monopole_plaket[0][index]);
  data_pattern.move_backward(1, 3);
  data_pattern_x.move_backward(1, 3);
  data_pattern_y.move_backward(1, 3);
  data_pattern_z.move_backward(1, 3);
  data_pattern.move_forward(1, 0);
  data_pattern_t.move_forward(1, 0);
  data_pattern_y.move_forward(1, 0);
  data_pattern_z.move_forward(1, 0);
  index = data_pattern.get_index_site();
  j0 = -(-(monopole_plaket[3][data_pattern_t.get_index_site()] -
           monopole_plaket[3][index]) -
         (monopole_plaket[5][data_pattern_y.get_index_site()] -
          monopole_plaket[5][index]) +
         monopole_plaket[4][data_pattern_z.get_index_site()] -
         monopole_plaket[4][index]);
  data_pattern.move_backward(1, 0);
  data_pattern_t.move_backward(1, 0);
  data_pattern_y.move_backward(1, 0);
  data_pattern_z.move_backward(1, 0);
  data_pattern.move_forward(1, 1);
  data_pattern_t.move_forward(1, 1);
  data_pattern_x.move_forward(1, 1);
  data_pattern_z.move_forward(1, 1);
  index = data_pattern.get_index_site();
  j1 = -(monopole_plaket[1][data_pattern_t.get_index_site()] -
         monopole_plaket[1][index] +
         monopole_plaket[5][data_pattern_x.get_index_site()] -
         monopole_plaket[5][index] -
         (monopole_plaket[2][data_pattern_z.get_index_site()] -
          monopole_plaket[2][index]));
  data_pattern.move_backward(1, 1);
  data_pattern_t.move_backward(1, 1);
  data_pattern_x.move_backward(1, 1);
  data_pattern_z.move_backward(1, 1);
  data_pattern.move_forward(1, 2);
  data_pattern_t.move_forward(1, 2);
  data_pattern_x.move_forward(1, 2);
  data_pattern_y.move_forward(1, 2);
  index = data_pattern.get_index_site();
  j2 = -(-(monopole_plaket[0][data_pattern_t.get_index_site()] -
           monopole_plaket[0][index]) -
         (monopole_plaket[4][data_pattern_x.get_index_site()] -
          monopole_plaket[4][index]) +
         monopole_plaket[2][data_pattern_y.get_index_site()] -
         monopole_plaket[2][index]);
  data_pattern.move_backward(1, 2);
  data_pattern_t.move_backward(1, 2);
  data_pattern_x.move_backward(1, 2);
  data_pattern_y.move_backward(1, 2);
}

// read configuration of angles
std::vector<double> read_angles_float_fortran(std::string &file_path) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  std::vector<double> angles;
  angles.reserve(data_size1);
  std::ifstream stream(file_path);
  std::vector<float> v(data_size1 + 2);
  if (!stream.read((char *)&v[0], (data_size1 + 2) * sizeof(float)))
    std::cout << "read_angles_float_fortran error: " << file_path << std::endl;
  for (int i = 0; i < data_size1; i++) {
    angles.push_back((double)v[i + 1]);
  }
  stream.close();
  return angles;
}

std::vector<double> read_angles_double_fortran(std::string &file_path) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  std::vector<double> angles;
  angles.reserve(data_size1);
  std::ifstream stream(file_path);
  std::vector<double> v(data_size1);
  stream.ignore(4);
  if (!stream.read((char *)&v[0], (data_size1) * sizeof(double)))
    std::cout << "read_angles_double_fortran error: " << file_path << std::endl;
  for (int i = 0; i < data_size1; i++) {
    angles.push_back((double)v[i]);
  }
  stream.close();
  return angles;
}

// read configuration of su2 matrices and extract abelian degrees of freedom
std::vector<double> read_float_fortran_convet_abelian(std::string &file_path) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  std::vector<double> angles;
  angles.reserve(data_size1 * 4);
  std::ifstream stream(file_path);
  std::vector<float> v(data_size1 * 4 + 1);
  if (!stream.read((char *)&v[0], (data_size1 * 4 + 1) * sizeof(float)))
    std::cout << "read_float_fortran_convet_abelian error: " << file_path
              << std::endl;
  for (int i = 0; i < data_size1; i++) {
    angles.push_back((double)atan2(v[i * 4 + 4], v[i * 4 + 1]));
  }
  stream.close();
  return angles;
}

std::vector<double> read_double_fortran_convet_abelian(std::string &file_path) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  std::vector<double> angles;
  angles.reserve(data_size1 * 4);
  std::ifstream stream(file_path);
  std::vector<double> v(data_size1 * 4);
  stream.ignore(4);
  if (!stream.read((char *)&v[0], (data_size1 * 4) * sizeof(double)))
    std::cout << "read_double_fortran_convet_abelian error: " << file_path
              << std::endl;
  for (int i = 0; i < data_size1; i++) {
    angles.push_back((double)atan2(v[i * 4 + 3], v[i * 4]));
  }
  stream.close();
  return angles;
}

std::vector<std::vector<double>>
read_double_su3_convet_angles(std::string &file_path) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  std::vector<std::vector<double>> angles(3, std::vector<double>(data_size));

  std::ifstream stream(file_path);
  std::vector<double> v(data_size * 18);

  if (!stream.read((char *)&v[0], (data_size * 18) * sizeof(double)))
    std::cout << "read_double_su3_convet_angles error: " << file_path
              << std::endl;

  for (int j = 0; j < data_size; j++) {
    for (int i = 0; i < 3; i++) {
      angles[i][j] = atan2(v[j * 18 + i * 8 + 1], v[j * 18 + i * 8]);
    }
  }
  stream.close();
  return angles;
}

std::vector<double> convert_to_angles(std::vector<abelian> &conf_abelian) {
  std::vector<double> angles(conf_abelian.size());

  for (int i = 0; i < conf_abelian.size(); i++) {
    angles[i] = conf_abelian[i].phi;
  }

  return angles;
}

std::vector<std::vector<double>>
convert_to_angles(std::vector<su3_abelian> &conf) {
  std::vector<std::vector<double>> angles(3, std::vector<double>(conf.size()));

  for (int i = 0; i < conf.size(); i++) {
    for (int j = 0; j < 3; j++) {
      angles[j][i] = atan2(conf[i].matrix[j].imag(), conf[i].matrix[j].real());
    }
  }

  return angles;
}

std::vector<std::vector<double>> convert_to_angles(std::vector<su3> &conf) {
  std::vector<std::vector<double>> angles(3, std::vector<double>(conf.size()));

  for (int i = 0; i < conf.size(); i++) {
    for (int j = 0; j < 3; j++) {
      angles[j][i] =
          atan2(conf[i].matrix(j, j).imag(), conf[i].matrix(j, j).real());
    }
  }

  return angles;
}

std::vector<double>
read_double_qc2dstag_convet_abelian(std::string &file_path) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  std::vector<double> angles;
  angles.reserve(data_size1 * 4);
  std::ifstream stream(file_path);
  std::vector<double> v(data_size1 * 4);
  if (!stream.read((char *)&v[0], (data_size1 * 4) * sizeof(double)))
    std::cout << "read_float_convert_abelian<abelian> error: " << file_path
              << std::endl;
  int dir;
  SPACE_ITER_START_SIMPLE
  for (int dir1 = 1; dir1 <= 4; dir1++) {
    if (dir1 == 4)
      dir = 1;
    else
      dir = dir1 + 1;
    angles.push_back((double)atan2(v[PLACE_QC2DSTAG + 3], v[PLACE_QC2DSTAG]));
  }
  SPACE_ITER_END_SIMPLE
  stream.close();
  return angles;
}

// calculate monopole_plaketes on the lattice
std::vector<std::vector<double>>
calculate_monopole_plaket(std::vector<double> &angles) {
  int data_size = x_size * y_size * z_size * t_size;
  std::vector<std::vector<double>> plakets(6, std::vector<double>(data_size));
  link1 link(x_size, y_size, z_size, t_size);
  int count = 0;
  for (int mu = 0; mu < 4; mu++) {
    link.move_dir(mu);
    for (int nu = mu + 1; nu < 4; nu++) {
      SPACE_ITER_START

      plakets[count][link.place / 4] = link.monopole_plaket_mu(angles, nu);

      SPACE_ITER_END

      count++;
    }
  }
  return plakets;
}

std::vector<std::vector<std::vector<double>>> calculate_monopole_plaket(
    const Data::LatticeData<DataPatternLexicographical, su3_angles> &conf) {
  DataPatternLexicographical data_pattern(conf.lat_dim);
  std::vector<std::vector<std::vector<double>>> plakets(
      3, std::vector<std::vector<double>>(
             6, std::vector<double>(data_pattern.get_lattice_size())));
  int count = 0;
  int index;
  double plaket;
  su3_angles A;
  for (int mu = 0; mu < 4; mu++) {
    for (int nu = mu + 1; nu < 4; nu++) {
      for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
        for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
          for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
            for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
              data_pattern.lat_coord = {x, y, z, t};
              index = data_pattern.get_index_site();
              double plaket = 0;
              A = conf[data_pattern.get_index_link(mu)];
              data_pattern.move_forward(1, mu);
              A = A * conf[data_pattern.get_index_link(nu)];
              data_pattern.move_backward(1, mu);
              data_pattern.move_forward(1, nu);
              A = A ^ conf[data_pattern.get_index_link(mu)];
              data_pattern.move_backward(1, nu);
              A = A ^ conf[data_pattern.get_index_link(nu)];
              for (int i = 0; i < 3; i++) {
                plaket = A.matrix[i];
                while (plaket > M_PI) {
                  plaket -= 2 * M_PI;
                }
                while (plaket < -M_PI) {
                  plaket += 2 * M_PI;
                }
                plakets[i][count][index] = plaket;
              }
            }
          }
        }
      }
      count++;
    }
  }
  return plakets;
}

std::vector<std::vector<double>> calculate_monopole_plaket(
    const Data::LatticeData<DataPatternLexicographical, abelian> &conf) {
  DataPatternLexicographical data_pattern(conf.lat_dim);
  std::vector<std::vector<double>> plakets(
      6, std::vector<double>(data_pattern.get_lattice_size()));
  int count = 0;
  int index;
  double plaket;
  abelian A;
  for (int mu = 0; mu < 4; mu++) {
    for (int nu = mu + 1; nu < 4; nu++) {
      for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
        for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
          for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
            for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
              data_pattern.lat_coord = {x, y, z, t};
              index = data_pattern.get_index_site();
              double plaket = 0;
              A = conf[data_pattern.get_index_link(mu)];
              data_pattern.move_forward(1, mu);
              A = A * conf[data_pattern.get_index_link(nu)];
              data_pattern.move_backward(1, mu);
              data_pattern.move_forward(1, nu);
              A = A ^ conf[data_pattern.get_index_link(mu)];
              data_pattern.move_backward(1, nu);
              A = A ^ conf[data_pattern.get_index_link(nu)];
              plaket = A.phi;
              while (plaket > M_PI) {
                plaket -= 2 * M_PI;
              }
              while (plaket < -M_PI) {
                plaket += 2 * M_PI;
              }
              plakets[count][index] = plaket;
            }
          }
        }
      }
      count++;
    }
  }
  return plakets;
}

std::vector<std::vector<std::vector<int>>> calculate_monopole_plaket_singular(
    const Data::LatticeData<DataPatternLexicographical, su3_angles> &conf) {
  DataPatternLexicographical data_pattern(conf.lat_dim);
  std::vector<std::vector<std::vector<int>>> plakets(
      3, std::vector<std::vector<int>>(
             6, std::vector<int>(data_pattern.get_lattice_size())));
  int count = 0;
  int index;
  double plaket;
  int plaket_singular;
  su3_angles A;
  for (int mu = 0; mu < 4; mu++) {
    for (int nu = mu + 1; nu < 4; nu++) {
      for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
        for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
          for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
            for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
              data_pattern.lat_coord = {x, y, z, t};
              index = data_pattern.get_index_site();
              double plaket = 0;
              A = conf[data_pattern.get_index_link(mu)];
              data_pattern.move_forward(1, mu);
              A = A * conf[data_pattern.get_index_link(nu)];
              data_pattern.move_backward(1, mu);
              data_pattern.move_forward(1, nu);
              A = A ^ conf[data_pattern.get_index_link(mu)];
              data_pattern.move_backward(1, nu);
              A = A ^ conf[data_pattern.get_index_link(nu)];
              for (int i = 0; i < 3; i++) {
                plaket = A.matrix[i];
                plaket_singular = 0;
                while (plaket > M_PI) {
                  plaket -= 2 * M_PI;
                  plaket_singular++;
                }
                while (plaket < -M_PI) {
                  plaket += 2 * M_PI;
                  plaket_singular--;
                }
                plakets[i][count][index] = plaket_singular;
              }
            }
          }
        }
      }
      count++;
    }
  }
  return plakets;
}

std::vector<std::vector<int>> calculate_monopole_plaket_singular(
    const Data::LatticeData<DataPatternLexicographical, abelian> &conf) {
  DataPatternLexicographical data_pattern(conf.lat_dim);
  std::vector<std::vector<int>> plakets(
      6, std::vector<int>(data_pattern.get_lattice_size()));
  int count = 0;
  int index;
  double plaket;
  int plaket_singular;
  abelian A;
  for (int mu = 0; mu < 4; mu++) {
    for (int nu = mu + 1; nu < 4; nu++) {
      for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
        for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
          for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
            for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
              data_pattern.lat_coord = {x, y, z, t};
              index = data_pattern.get_index_site();
              double plaket = 0;
              A = conf[data_pattern.get_index_link(mu)];
              data_pattern.move_forward(1, mu);
              A = A * conf[data_pattern.get_index_link(nu)];
              data_pattern.move_backward(1, mu);
              data_pattern.move_forward(1, nu);
              A = A ^ conf[data_pattern.get_index_link(mu)];
              data_pattern.move_backward(1, nu);
              A = A ^ conf[data_pattern.get_index_link(nu)];
              plaket = A.phi;
              plaket_singular = 0;
              while (plaket > M_PI) {
                plaket -= 2 * M_PI;
                plaket_singular++;
              }
              while (plaket < -M_PI) {
                plaket += 2 * M_PI;
                plaket_singular--;
              }
              plakets[count][index] = plaket_singular;
            }
          }
        }
      }
      count++;
    }
  }
  return plakets;
}

std::vector<std::vector<int>>
calculate_monopole_plaket_singular(std::vector<double> &angles) {
  int data_size = x_size * y_size * z_size * t_size;
  std::vector<std::vector<int>> singular(6, std::vector<int>(data_size));
  link1 link(x_size, y_size, z_size, t_size);

  int count = 0;
  for (int mu = 0; mu < 4; mu++) {
    link.move_dir(mu);
    for (int nu = mu + 1; nu < 4; nu++) {
      SPACE_ITER_START

      singular[count][link.place / 4] =
          link.monopole_plaket_singular_mu(angles, nu);

      SPACE_ITER_END

      count++;
    }
  }
  return singular;
}

std::vector<std::vector<int>>
calculate_monopole_plaket_singular(std::vector<double> &angles,
                                   DataPatternLexicographical &data_pattern) {
  int data_size = data_pattern.get_lattice_size();
  std::vector<std::vector<int>> singular(6, std::vector<int>(data_size));
  int count = 0;
  int index;
  double plaket;
  int plaket_singular;
  for (int mu = 0; mu < 4; mu++) {
    for (int nu = mu + 1; nu < 4; nu++) {
      for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
        for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
          for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
            for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
              data_pattern.lat_coord = {x, y, z, t};
              index = data_pattern.get_index_site();
              plaket = 0;
              plaket = angles[data_pattern.get_index_link(mu)];
              data_pattern.move_forward(1, mu);
              plaket += angles[data_pattern.get_index_link(nu)];
              data_pattern.move_backward(1, mu);
              data_pattern.move_forward(1, nu);
              plaket -= angles[data_pattern.get_index_link(mu)];
              data_pattern.move_backward(1, nu);
              plaket -= angles[data_pattern.get_index_link(nu)];
              plaket_singular = 0;
              while (plaket > M_PI) {
                plaket -= 2 * M_PI;
                plaket_singular++;
              }
              while (plaket < -M_PI) {
                plaket += 2 * M_PI;
                plaket_singular--;
              }
              singular[count][index] = plaket_singular;
            }
          }
        }
      }

      count++;
    }
  }
  return singular;
}

std::vector<double> calculate_current(std::vector<double> &angles) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  link1 link(x_size, y_size, z_size, t_size);

  std::vector<double> J(data_size);

  std::vector<std::vector<double>> monopole_plaket =
      calculate_monopole_plaket(angles);

  SPACE_ITER_START

  link.get_current(monopole_plaket, &J[link.place]);

  SPACE_ITER_END

  return J;
}

std::vector<std::vector<double>> calculate_current(
    const Data::LatticeData<DataPatternLexicographical, su3_angles> &conf) {
  DataPatternLexicographical data_pattern(conf.lat_dim);
  std::vector<std::vector<std::vector<double>>> monopole_plakets =
      make_monopole_plakets(conf);
  std::vector<std::vector<double>> J(
      3, std::vector<double>(data_pattern.get_data_size()));
  for (int i = 0; i < 3; i++) {
    for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
      for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
        for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
          for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
            data_pattern.lat_coord = {x, y, z, t};
            get_current(monopole_plakets[i],
                        &J[i][data_pattern.get_index_link(0)], data_pattern);
          }
        }
      }
    }
  }
  return J;
}

std::vector<double> calculate_current(
    const Data::LatticeData<DataPatternLexicographical, abelian> &conf) {
  DataPatternLexicographical data_pattern(conf.lat_dim);
  std::vector<std::vector<double>> monopole_plakets =
      calculate_monopole_plaket(conf);
  std::vector<double> J(data_pattern.get_data_size());
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          get_current(monopole_plakets, &J[data_pattern.get_index_link(0)],
                      data_pattern);
        }
      }
    }
  }
  return J;
}

std::vector<std::vector<double>> calculate_current_singular(
    const Data::LatticeData<DataPatternLexicographical, su3_angles> &conf) {
  DataPatternLexicographical data_pattern(conf.lat_dim);
  std::vector<std::vector<std::vector<int>>> monopole_plakets_singular =
      make_monopole_plakets_singular(conf);
  std::vector<std::vector<double>> J(
      3, std::vector<double>(data_pattern.get_data_size()));
  for (int i = 0; i < 3; i++) {
    for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
      for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
        for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
          for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
            data_pattern.lat_coord = {x, y, z, t};
            get_current_singular(monopole_plakets_singular[i],
                                 &J[i][data_pattern.get_index_link(0)],
                                 data_pattern);
          }
        }
      }
    }
  }
  return J;
}

std::vector<double> calculate_current_singular(
    const Data::LatticeData<DataPatternLexicographical, abelian> &conf) {
  DataPatternLexicographical data_pattern(conf.lat_dim);
  std::vector<std::vector<int>> monopole_plaket_singular =
      calculate_monopole_plaket_singular(conf);
  std::vector<double> J(data_pattern.get_data_size());
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          get_current_singular(monopole_plaket_singular,
                               &J[data_pattern.get_index_link(0)],
                               data_pattern);
        }
      }
    }
  }
  return J;
}

std::vector<int> calculate_current_singular(std::vector<double> &angles) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  link1 link(x_size, y_size, z_size, t_size);

  std::vector<int> J(data_size);

  std::vector<std::vector<int>> monopole_plaket_singular =
      calculate_monopole_plaket_singular(angles);

  SPACE_ITER_START

  link.get_current_singular(monopole_plaket_singular, &J[link.place]);

  SPACE_ITER_END

  return J;
}

std::vector<double> calculate_current_monopole_plakets(
    std::vector<std::vector<double>> monopole_plakets) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  link1 link(x_size, y_size, z_size, t_size);

  std::vector<double> J(data_size);

  SPACE_ITER_START

  link.get_current(monopole_plakets, &J[link.place]);

  SPACE_ITER_END

  return J;
}

std::vector<double> calculate_current_monopole_plakets(
    std::vector<std::vector<double>> monopole_plakets,
    DataPatternLexicographical data_pattern) {
  std::vector<double> J(data_pattern.get_data_size());
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          get_current(monopole_plakets, &J[data_pattern.get_index_link(0)],
                      data_pattern);
        }
      }
    }
  }
  return J;
}

std::vector<std::vector<std::vector<double>>>
make_monopole_plakets(std::vector<std::vector<double>> &angles) {
  link1 link(x_size, y_size, z_size, t_size);

  std::vector<std::vector<std::vector<double>>> monopole_plakets(3);

  for (int i = 0; i < 3; i++) {
    monopole_plakets[i] = calculate_monopole_plaket(angles[i]);
  }

  double sum;
  int place;
  double extremal;
  for (int j = 0; j < monopole_plakets[0].size(); j++) {
    for (int k = 0; k < monopole_plakets[0][0].size(); k++) {

      sum = 0;
      for (int i = 0; i < 3; i++) {
        sum += monopole_plakets[i][j][k];
      }

      if (sum >= M_PI) {
        place = 0;
        extremal = monopole_plakets[0][j][k];

        if (monopole_plakets[1][j][k] >= extremal) {
          extremal = monopole_plakets[1][j][k];
          place = 1;
        }

        if (monopole_plakets[2][j][k] >= extremal) {
          place = 2;
        }

        monopole_plakets[place][j][k] -= 2 * M_PI;
      }

      if (sum <= -M_PI) {
        place = 0;
        extremal = monopole_plakets[0][j][k];

        if (monopole_plakets[1][j][k] <= extremal) {
          extremal = monopole_plakets[1][j][k];
          place = 1;
        }

        if (monopole_plakets[2][j][k] <= extremal) {
          place = 2;
        }

        monopole_plakets[place][j][k] += 2 * M_PI;
      }
    }
  }

  return monopole_plakets;
}

std::vector<std::vector<std::vector<double>>> make_monopole_plakets(
    const Data::LatticeData<DataPatternLexicographical, su3_angles> &conf) {
  std::vector<std::vector<std::vector<double>>> monopole_plakets =
      calculate_monopole_plaket(conf);
  double sum;
  int place;
  double extremal;
  for (int j = 0; j < monopole_plakets[0].size(); j++) {
    for (int k = 0; k < monopole_plakets[0][0].size(); k++) {
      sum = 0;
      for (int i = 0; i < 3; i++) {
        sum += monopole_plakets[i][j][k];
      }
      if (sum >= M_PI) {
        place = 0;
        extremal = monopole_plakets[0][j][k];
        if (monopole_plakets[1][j][k] >= extremal) {
          extremal = monopole_plakets[1][j][k];
          place = 1;
        }
        if (monopole_plakets[2][j][k] >= extremal) {
          place = 2;
        }
        monopole_plakets[place][j][k] -= 2 * M_PI;
      }
      if (sum <= -M_PI) {
        place = 0;
        extremal = monopole_plakets[0][j][k];
        if (monopole_plakets[1][j][k] <= extremal) {
          extremal = monopole_plakets[1][j][k];
          place = 1;
        }
        if (monopole_plakets[2][j][k] <= extremal) {
          place = 2;
        }
        monopole_plakets[place][j][k] += 2 * M_PI;
      }
    }
  }
  return monopole_plakets;
}

std::vector<std::vector<std::vector<int>>>
make_monopole_plakets_singular(std::vector<std::vector<double>> &angles) {
  link1 link(x_size, y_size, z_size, t_size);

  std::vector<std::vector<std::vector<int>>> dirac_plakets(3);

  for (int i = 0; i < 3; i++) {
    dirac_plakets[i] = calculate_monopole_plaket_singular(angles[i]);
  }

  double sum;
  int place;
  double extremal;
  for (int j = 0; j < dirac_plakets[0].size(); j++) {
    for (int k = 0; k < dirac_plakets[0][0].size(); k++) {

      sum = 0;
      for (int i = 0; i < 3; i++) {
        sum += dirac_plakets[i][j][k];
      }

      if (sum >= 1) {
        place = 0;
        extremal = dirac_plakets[0][j][k];

        if (dirac_plakets[1][j][k] >= extremal) {
          extremal = dirac_plakets[1][j][k];
          place = 1;
        }

        if (dirac_plakets[2][j][k] >= extremal) {
          place = 2;
        }

        dirac_plakets[place][j][k] -= 1;
      }

      if (sum <= -1) {
        place = 0;
        extremal = dirac_plakets[0][j][k];

        if (dirac_plakets[1][j][k] <= extremal) {
          extremal = dirac_plakets[1][j][k];
          place = 1;
        }

        if (dirac_plakets[2][j][k] <= extremal) {
          place = 2;
        }

        dirac_plakets[place][j][k] += 1;
      }
    }
  }

  return dirac_plakets;
}

std::vector<std::vector<std::vector<int>>> make_monopole_plakets_singular(
    const Data::LatticeData<DataPatternLexicographical, su3_angles> &conf) {
  std::vector<std::vector<std::vector<int>>> dirac_plakets =
      calculate_monopole_plaket_singular(conf);
  double sum;
  int place;
  double extremal;
  for (int j = 0; j < dirac_plakets[0].size(); j++) {
    for (int k = 0; k < dirac_plakets[0][0].size(); k++) {
      sum = 0;
      for (int i = 0; i < 3; i++) {
        sum += dirac_plakets[i][j][k];
      }
      if (sum >= 1) {
        place = 0;
        extremal = dirac_plakets[0][j][k];
        if (dirac_plakets[1][j][k] >= extremal) {
          extremal = dirac_plakets[1][j][k];
          place = 1;
        }
        if (dirac_plakets[2][j][k] >= extremal) {
          place = 2;
        }
        dirac_plakets[place][j][k] -= 1;
      }
      if (sum <= -1) {
        place = 0;
        extremal = dirac_plakets[0][j][k];
        if (dirac_plakets[1][j][k] <= extremal) {
          extremal = dirac_plakets[1][j][k];
          place = 1;
        }
        if (dirac_plakets[2][j][k] <= extremal) {
          place = 2;
        }
        dirac_plakets[place][j][k] += 1;
      }
    }
  }
  return dirac_plakets;
}

void make_plakets_both(
    std::vector<std::vector<double>> &angles,
    std::vector<std::vector<std::vector<double>>> &monopole_plakets,
    std::vector<std::vector<std::vector<int>>> &dirac_plakets) {
  link1 link(x_size, y_size, z_size, t_size);

  std::vector<std::vector<std::vector<double>>> monopole_plakets_tmp(3);

  for (int i = 0; i < 3; i++) {
    monopole_plakets_tmp[i] = calculate_monopole_plaket(angles[i]);
  }

  std::vector<std::vector<std::vector<int>>> dirac_plakets_tmp(3);

  for (int i = 0; i < 3; i++) {
    dirac_plakets_tmp[i] = calculate_monopole_plaket_singular(angles[i]);
  }

  double sum;
  int place;
  double extremal;
  for (int j = 0; j < monopole_plakets_tmp[0].size(); j++) {
    for (int k = 0; k < monopole_plakets_tmp[0][0].size(); k++) {

      sum = 0;
      for (int i = 0; i < 3; i++) {
        sum += monopole_plakets_tmp[i][j][k];
      }

      if (sum >= M_PI) {
        place = 0;
        extremal = monopole_plakets_tmp[0][j][k];

        if (monopole_plakets_tmp[1][j][k] >= extremal) {
          extremal = monopole_plakets_tmp[1][j][k];
          place = 1;
        }

        if (monopole_plakets_tmp[2][j][k] >= extremal) {
          place = 2;
        }

        monopole_plakets_tmp[place][j][k] -= 2 * M_PI;
        dirac_plakets_tmp[place][j][k] += 1;
      }

      if (sum <= -M_PI) {
        place = 0;
        extremal = monopole_plakets_tmp[0][j][k];

        if (monopole_plakets_tmp[1][j][k] <= extremal) {
          extremal = monopole_plakets_tmp[1][j][k];
          place = 1;
        }

        if (monopole_plakets_tmp[2][j][k] <= extremal) {
          place = 2;
        }

        monopole_plakets_tmp[place][j][k] += 2 * M_PI;
        dirac_plakets_tmp[place][j][k] -= 1;
      }
    }
  }

  monopole_plakets = std::move(monopole_plakets_tmp);
  dirac_plakets = std::move(dirac_plakets_tmp);
}

void make_plakets_both(
    const Data::LatticeData<DataPatternLexicographical, su3_angles> &conf,
    std::vector<std::vector<std::vector<double>>> &monopole_plakets,
    std::vector<std::vector<std::vector<int>>> &dirac_plakets) {
  std::vector<std::vector<std::vector<double>>> monopole_plakets_tmp =
      calculate_monopole_plaket(conf);
  std::vector<std::vector<std::vector<int>>> dirac_plakets_tmp =
      calculate_monopole_plaket_singular(conf);
  double sum;
  int place;
  double extremal;
  for (int j = 0; j < monopole_plakets_tmp[0].size(); j++) {
    for (int k = 0; k < monopole_plakets_tmp[0][0].size(); k++) {
      sum = 0;
      for (int i = 0; i < 3; i++) {
        sum += monopole_plakets_tmp[i][j][k];
      }
      if (sum >= M_PI) {
        place = 0;
        extremal = monopole_plakets_tmp[0][j][k];
        if (monopole_plakets_tmp[1][j][k] >= extremal) {
          extremal = monopole_plakets_tmp[1][j][k];
          place = 1;
        }
        if (monopole_plakets_tmp[2][j][k] >= extremal) {
          place = 2;
        }
        monopole_plakets_tmp[place][j][k] -= 2 * M_PI;
        dirac_plakets_tmp[place][j][k] += 1;
      }
      if (sum <= -M_PI) {
        place = 0;
        extremal = monopole_plakets_tmp[0][j][k];
        if (monopole_plakets_tmp[1][j][k] <= extremal) {
          extremal = monopole_plakets_tmp[1][j][k];
          place = 1;
        }
        if (monopole_plakets_tmp[2][j][k] <= extremal) {
          place = 2;
        }
        monopole_plakets_tmp[place][j][k] += 2 * M_PI;
        dirac_plakets_tmp[place][j][k] -= 1;
      }
    }
  }
  monopole_plakets = std::move(monopole_plakets_tmp);
  dirac_plakets = std::move(dirac_plakets_tmp);
}

// returns 0 if no current has been found, direction +-1..4 if it has
template <class T> int find_current(link1 &link, std::vector<T> &J) {
  for (int mu = 0; mu < 4; mu++) {
    if ((J[link.place + mu] > 0.3) || (J[link.place + mu] < -0.3))
      return mu + 1;
  }
  for (int mu = 0; mu < 4; mu++) {
    link.move(mu, -1);
    if ((J[link.place + mu] > 0.3) || (J[link.place + mu] < -0.3))
      return -mu - 1;
    link.move(mu, 1);
  }
  return 0;
}

template <class T>
int find_current(const std::vector<T> &J,
                 DataPatternLexicographical &data_pattern) {
  int index;
  for (int mu = 0; mu < 4; mu++) {
    index = data_pattern.get_index_link(mu);
    if ((J[index] > 0.3) || (J[index] < -0.3))
      return mu + 1;
  }
  for (int mu = 0; mu < 4; mu++) {
    data_pattern.move_backward(1, mu);
    index = data_pattern.get_index_link(mu);
    data_pattern.move_forward(1, mu);
    if ((J[index] > 0.3) || (J[index] < -0.3))
      return -mu - 1;
  }
  return 0;
}

// find all directions with current and add loops with that neighbours
template <class T>
std::vector<loop *> find_paths(std::vector<loop *> &neighbours,
                               std::vector<T> &J) {
  // std::vector for new sites with current
  std::vector<loop *> neighbours_new;
  loop *loop_tmp;
  link1 link(x_size, y_size, z_size, t_size);
  int J_tmp;
  // for all previously found sites of cluster find new sites with current
  for (int i = 0; i < neighbours.size(); i++) {
    // std::cout << i << std::endl;
    // go to site
    link.go_update(neighbours[i]->coordinate[0], neighbours[i]->coordinate[1],
                   neighbours[i]->coordinate[2], neighbours[i]->coordinate[3]);
    // check all directions for a current
    // if a current is found, add it to cluster and to std::vector of new sites
    for (int mu = 0; mu < 4; mu++) {
      if (J[link.place + mu] > 0.3) {
        J_tmp = std::lround(J[link.place + mu]);
        J[link.place + mu] = 0.;
        link.move(mu, 1);
        loop_tmp = new loop(link);
        neighbours[i]->link.push_back(loop_tmp);
        neighbours[i]->charge.push_back(J_tmp);
        neighbours_new.push_back(loop_tmp);
        link.move(mu, -1);
      }
    }
    for (int mu = 0; mu < 4; mu++) {
      link.move(mu, -1);
      if (J[link.place + mu] < -0.3) {
        J_tmp = std::lround(J[link.place + mu]);
        J[link.place + mu] = 0.;
        loop_tmp = new loop(link);
        neighbours[i]->link.push_back(loop_tmp);
        neighbours[i]->charge.push_back(J_tmp);
        neighbours_new.push_back(loop_tmp);
      }
      link.move(mu, 1);
    }
  }
  return neighbours_new;
}

// find all directions with current and add loops with that neighbours
template <class T>
std::vector<loop_new *> find_paths(std::vector<loop_new *> &neighbours,
                                   std::vector<T> &J,
                                   DataPatternLexicographical &data_pattern) {
  // std::vector for new sites with current
  std::vector<loop_new *> neighbours_new;
  loop_new *loop_tmp;
  int J_tmp;
  int index;
  // for all previously found sites of cluster find new sites with current
  for (int i = 0; i < neighbours.size(); i++) {
    // go to site
    data_pattern.lat_coord = neighbours[i]->coordinate;
    // check all directions for a current
    // if a current is found, add it to cluster and to std::vector of new sites
    for (int mu = 0; mu < 4; mu++) {
      index = data_pattern.get_index_link(mu);
      if (J[index] > 0.3) {
        J_tmp = std::lround(J[index]);
        J[index] = 0.;
        data_pattern.move_forward(1, mu);
        loop_tmp = new loop_new(data_pattern.lat_coord);
        neighbours[i]->nodes.push_back(loop_tmp);
        neighbours[i]->charge.push_back(J_tmp);
        neighbours_new.push_back(loop_tmp);
        data_pattern.move_backward(1, mu);
      }
    }
    for (int mu = 0; mu < 4; mu++) {
      data_pattern.move_backward(1, mu);
      index = data_pattern.get_index_link(mu);
      if (J[index] < -0.3) {
        J_tmp = std::lround(J[index]);
        J[index] = 0.;
        loop_tmp = new loop_new(data_pattern.lat_coord);
        neighbours[i]->nodes.push_back(loop_tmp);
        neighbours[i]->charge.push_back(J_tmp);
        neighbours_new.push_back(loop_tmp);
      }
      data_pattern.move_forward(1, mu);
    }
  }
  return neighbours_new;
}

// find cluster which has site ll
template <class T> void find_cluster(loop *ll, std::vector<T> &J) {
  std::vector<loop *> neighbours = {ll};
  // while find_path finds new sites with current
  int count = 0;
  do {
    // find new sites with current and add them to loops
    neighbours = find_paths(neighbours, J);
  } while (neighbours.size() > 0);
}

template <class T>
void find_cluster(loop_new *ll, std::vector<T> &J,
                  DataPatternLexicographical &data_pattern) {
  std::vector<loop_new *> neighbours = {ll};
  // while find_path finds new sites with current
  do {
    // find new sites with current and add them to loops
    neighbours = find_paths(neighbours, J, data_pattern);
  } while (neighbours.size() > 0);
}

// find all clusters on a lattice not using recurrence
template <class T> std::vector<loop *> calculate_clusters(std::vector<T> &J) {
  int dir1;
  std::vector<loop *> LL;

  link1 link(x_size, y_size, z_size, t_size);

  SPACE_ITER_START

  // find site with at least one non-zero current
  dir1 = find_current(link, J);
  // if there's a current
  if (dir1 != 0) {

    // create a new loop and add it to cluster std::vector
    LL.push_back(new loop(link));

    // start finding cluster starting at this point
    find_cluster(LL[LL.size() - 1], J);
  }

  SPACE_ITER_END

  return LL;
}

template <class T>
std::vector<loop_new *>
calculate_clusters(std::vector<T> &J,
                   DataPatternLexicographical &data_pattern) {
  int dir1;
  std::vector<loop_new *> LL;
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          // find site with at least one non-zero current
          dir1 = find_current(J, data_pattern);
          // if there's a current
          if (dir1 != 0) {
            // create a new loop and add it to cluster std::vector
            LL.push_back(new loop_new(data_pattern.lat_coord));
            // start finding cluster starting at this point
            find_cluster(LL[LL.size() - 1], J, data_pattern);
          }
        }
      }
    }
  }
  return LL;
}

// functions for obtaining information about clusters for testing
void print_currents(loop *ll) {
  std::cout << ll->coordinate[0] << " " << ll->coordinate[1] << " "
            << ll->coordinate[2] << " " << ll->coordinate[3] << std::endl;
  for (int i = 0; i < ll->link.size(); i++) {
    print_currents(ll->link[i]);
  }
}

void print_currents(loop_new *ll) {
  std::cout << ll->coordinate[0] << " " << ll->coordinate[1] << " "
            << ll->coordinate[2] << " " << ll->coordinate[3] << std::endl;
  for (int i = 0; i < ll->nodes.size(); i++) {
    print_currents(ll->nodes[i]);
  }
}

void check_for_coordinate(loop *loop, int coordinate[4], bool &include) {
  if (loop->coordinate[0] == coordinate[0] &&
      loop->coordinate[1] == coordinate[1] &&
      loop->coordinate[2] == coordinate[2] &&
      loop->coordinate[3] == coordinate[3])
    include = true;
  for (int i = 0; i < loop->link.size(); i++) {
    check_for_coordinate(loop->link[i], coordinate, include);
  }
}

// monopole observables

int cluster_length(loop *ll) {
  int length = 0;
  link1 link(x_size, y_size, z_size, t_size);

  cluster_length_recurrent(ll, length);

  return length;
}

void cluster_length_recurrent(loop_new *ll, int &length) {
  for (int i = 0; i < ll->nodes.size(); i++) {
    cluster_length_recurrent(ll->nodes[i], length);
    length += abs(ll->charge[i]);
  }
}

int cluster_length(loop_new *ll) {
  int length = 0;
  cluster_length_recurrent(ll, length);
  return length;
}

void cluster_length_recurrent(loop *ll, int &length) {
  for (int i = 0; i < ll->link.size(); i++) {
    cluster_length_recurrent(ll->link[i], length);
    length += abs(ll->charge[i]);
  }
}

void cluster_sites(loop *ll) {
  std::cout << "x: " << ll->coordinate[0] << " y: " << ll->coordinate[1]
            << " z: " << ll->coordinate[2] << " t: " << ll->coordinate[3]
            << std::endl;
  for (int i = 0; i < ll->link.size(); i++) {
    cluster_sites(ll->link[i]);
  }
}

std::vector<int> length_mu(loop *ll) {
  std::vector<int> lengths_mu = {0, 0, 0, 0};

  length_mu_recurrent(ll, lengths_mu);

  return lengths_mu;
}

void length_mu_recurrent(loop_new *ll, std::vector<int> &lengths_mu) {
  for (int i = 0; i < ll->nodes.size(); i++) {
    length_mu_recurrent(ll->nodes[i], lengths_mu);
    int mu = 0;
    int difference;
    do {
      difference = ll->nodes[i]->coordinate[mu] - ll->coordinate[mu];
      mu++;
    } while (difference == 0);
    mu--;
    lengths_mu[mu] += ll->charge[i];
  }
}

std::vector<int> length_mu(loop_new *ll) {
  std::vector<int> lengths_mu = {0, 0, 0, 0};
  length_mu_recurrent(ll, lengths_mu);
  return lengths_mu;
}

void length_mu_recurrent(loop *ll, std::vector<int> &lengths_mu) {
  for (int i = 0; i < ll->link.size(); i++) {
    length_mu_recurrent(ll->link[i], lengths_mu);
    int mu = 0;
    int difference;
    do {
      difference = ll->link[i]->coordinate[mu] - ll->coordinate[mu];
      mu++;
    } while (difference == 0);
    mu--;
    lengths_mu[mu] += ll->charge[i];
  }
}

std::vector<int> currents_directions(loop *ll) {
  // [0] is spatial, [1] is temporal
  std::vector<int> directions = {0, 0};

  currents_directions_recurrent(ll, directions);

  return directions;
}

void currents_directions_recurrent(loop_new *ll,
                                   std::tuple<int, int> &directions) {
  for (int i = 0; i < ll->nodes.size(); i++) {
    currents_directions_recurrent(ll->nodes[i], directions);
    int mu = 0;
    int difference;
    do {
      difference = ll->nodes[i]->coordinate[mu] - ll->coordinate[mu];
      mu++;
    } while (difference == 0);
    mu--;
    if (mu == 3) {
      std::get<1>(directions) += abs(ll->charge[i]);
    } else {
      std::get<0>(directions) += abs(ll->charge[i]);
    }
  }
}

std::tuple<int, int> currents_directions(loop_new *ll) {
  // [0] is spatial, [1] is temporal
  std::tuple<int, int> directions = {0, 0};
  currents_directions_recurrent(ll, directions);
  return directions;
}

void currents_directions_recurrent(loop *ll, std::vector<int> &directions) {
  for (int i = 0; i < ll->link.size(); i++) {
    currents_directions_recurrent(ll->link[i], directions);
    int mu = 0;
    int difference;
    do {
      difference = ll->link[i]->coordinate[mu] - ll->coordinate[mu];
      mu++;
    } while (difference == 0);
    mu--;
    if (mu == 3) {
      directions[1] += abs(ll->charge[i]);
    } else {
      directions[0] += abs(ll->charge[i]);
    }
  }
}

double cluster_variation(loop *loop) {
  double variation = 0;

  std::vector<int> distance = {0, 0, 0};
  link1 link(x_size, y_size, z_size, t_size);

  cluster_variation_recurrent(loop, variation, distance, link);

  return variation;
}

void cluster_variation_recurrent(loop *loop, double &variation,
                                 std::vector<int> distance, const link1 &link) {
  for (int i = 0; i < loop->link.size(); i++) {
    int difference = 0;
    int mu = 0;
    do {
      if (loop->link[i]->coordinate[mu] == link.lattice_size[mu] - 1 &&
          loop->coordinate[mu] == 0)
        difference = -1;
      else if (loop->link[i]->coordinate[mu] == 0 &&
               loop->coordinate[mu] == link.lattice_size[mu] - 1)
        difference = 1;
      else
        difference = loop->link[i]->coordinate[mu] - loop->coordinate[mu];

      mu++;
    } while (difference == 0);
    mu--;
    if (mu != 3) {
      distance[mu] += difference;
    }

    variation += distance[0] * distance[0] + distance[1] * distance[1] +
                 distance[2] * distance[2];

    cluster_variation_recurrent(loop->link[i], variation, distance, link);

    if (mu != 3)
      distance[mu] -= difference;
  }
}

int site_number(loop *loop) {
  int link_number = 0;
  std::unordered_map<int, int> loop_sites;

  site_number_recurrent(loop, link_number, loop_sites);

  return link_number;
}

void site_number_recurrent(loop *loop, int &link_number,
                           std::unordered_map<int, int> &loop_sites) {
  for (int i = 0; i < loop->link.size(); i++) {
    site_number_recurrent(loop->link[i], link_number, loop_sites);
  }

  int key = 1000000 * loop->coordinate[0] + 10000 * loop->coordinate[1] +
            100 * loop->coordinate[2] + loop->coordinate[3];

  if (loop_sites.find(key) == loop_sites.end()) {
    loop_sites[key];
    link_number++;
  }
}

bool check_wrappings_sum(std::vector<std::vector<int>> &wrappings,
                         std::vector<int> &positions) {
  std::vector<int> result = {0, 0, 0, 0};
  for (int i = 0; i < positions.size(); i++) {
    std::transform(result.begin(), result.end(),
                   wrappings[positions[i]].begin(), result.begin(),
                   std::plus<int>());
  }
  return std::all_of(result.begin(), result.end(),
                     [](int i) { return i == 0; });
}

bool increment_position(std::vector<int> &positions, int max_position) {
  if (positions[positions.size() - 1] < max_position) {
    positions[positions.size() - 1]++;
    return true;
  } else {
    for (int i = positions.size() - 2; i > 0; i--) {
      if (positions[i + 1] - positions[i] > 2) {
        positions[i]++;
        for (int j = i + 1; j < positions.size(); j++) {
          positions[j] = positions[i] + j - i;
        }
        return true;
      }
    }
  }
  return false;
}

bool iterate_sum_order(std::vector<std::vector<int>> &wrappings,
                       std::vector<int> &positions) {
  do {
    if (check_wrappings_sum(wrappings, positions)) {
      return true;
    }
  } while (increment_position(positions, wrappings.size() - 1));
  return false;
}

std::vector<int> group_percolating(std::vector<std::vector<int>> &wrappings) {
  std::vector<int> percolating_positions;
  for (int sum_order = 2; sum_order <= wrappings.size(); sum_order++) {
    percolating_positions = std::vector<int>(sum_order);
    std::iota(percolating_positions.begin(), percolating_positions.end(), 0);
    if (iterate_sum_order(wrappings, percolating_positions)) {
      return percolating_positions;
    }
  }
  return std::vector<int>();
}

template int find_current(link1 &link, std::vector<double> &J);
template std::vector<loop *> find_paths(std::vector<loop *> &neighbours,
                                        std::vector<double> &J);
template void find_cluster(loop *ll, std::vector<double> &J);
template std::vector<loop *> calculate_clusters(std::vector<double> &J);
template std::vector<loop_new *>
calculate_clusters(std::vector<double> &J,
                   DataPatternLexicographical &data_pattern);

template int find_current(link1 &link, std::vector<int> &J);
template std::vector<loop *> find_paths(std::vector<loop *> &neighbours,
                                        std::vector<int> &J);
template void find_cluster(loop *ll, std::vector<int> &J);
template std::vector<loop *> calculate_clusters(std::vector<int> &J);
template std::vector<loop_new *>
calculate_clusters(std::vector<int> &J,
                   DataPatternLexicographical &data_pattern);