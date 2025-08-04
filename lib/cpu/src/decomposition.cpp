#include "../include/decomposition.h"
#include "../include/link.h"
#include "../include/matrix.h"
#include "../include/monopoles.h"

#include <ctime>
#include <fstream>
#include <math.h>
#include <omp.h>

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

std::vector<double> read_double_angles(std::string &file_name, int bites_skip) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  std::vector<double> angles(data_size);

  std::ifstream stream(file_name);
  stream.ignore(bites_skip);
  if (!stream.read((char *)&angles[0], (data_size) * sizeof(double)))
    std::cout << "read_double_angles error: " << file_name << std::endl;
  return angles;
}

void write_double_angles(std::string &file_name,
                         const std::vector<double> &angles,
                         DataPatternLexicographical &data_pattern) {
  int data_size = data_pattern.get_data_size();
  std::ofstream stream(file_name);
  if (!stream.write((char *)&angles[0], (data_size) * sizeof(double)))
    std::cout << "write_double_angles error: " << file_name << std::endl;
}

std::vector<std::vector<double>>
read_double_angles_su3(std::string &file_name,
                       DataPatternLexicographical &data_pattern) {
  int data_size = data_pattern.get_data_size();
  std::vector<std::vector<double>> angles(3, std::vector<double>(data_size));
  std::ifstream stream(file_name);
  for (int i = 0; i < 3; i++) {
    if (!stream.read((char *)&angles[i][0], (data_size) * sizeof(double)))
      std::cout << "read_double_angles_su3 error: " << file_name << std::endl;
  }
  return angles;
}

void write_double_angles_su3(std::string &file_name,
                             const std::vector<std::vector<double>> &angles,
                             DataPatternLexicographical &data_pattern) {
  int data_size = data_pattern.get_data_size();
  std::cout << data_size << std::endl;
  std::ofstream stream(file_name);
  for (int i = 0; i < 3; i++) {
    if (!stream.write((char *)&angles[i][0], (data_size) * sizeof(double)))
      std::cout << "write_double_angles_su3 error: " << file_name << std::endl;
  }
}

void write_double_su2(std::string &file_name,
                      const std::vector<su2> &conf_su2) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  std::ofstream stream(file_name);
  if (!stream.write((char *)&conf_su2[0], data_size * 4 * sizeof(double)))
    std::cout << "write_double_su2 error: " << file_name << std::endl;
  stream.close();
}

void write_double_su3(std::string &file_name,
                      const std::vector<su3> &conf_su3) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  std::ofstream stream(file_name);
  if (!stream.write((char *)&conf_su3[0], data_size * 18 * sizeof(double)))
    std::cout << "write_double_su3 error: " << file_name << std::endl;
  stream.close();
}

std::vector<std::vector<double>>
get_angles_su3(const std::vector<su3> &conf_su3) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  std::vector<std::vector<double>> angles =
      std::vector<std::vector<double>>(3, std::vector<double>(data_size));

  for (int i = 0; i < data_size; i++) {
    for (int c = 0; c < 3; c++) {
      angles[c][i] = atan2(conf_su3[i].matrix(c, c).imag(),
                           conf_su3[i].matrix(c, c).real());
    }
  }

  return angles;
}

std::vector<su3_abelian> get_abelian(std::vector<su3> &conf) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  std::vector<su3_abelian> abelian(data_size);
  for (int i = 0; i < data_size; i++) {
    std::vector<double> angles(3);
    double sum = 0;
    for (int c = 0; c < 3; c++) {
      angles[c] =
          atan2(conf[i].matrix(c, c).imag(), conf[i].matrix(c, c).real());
      sum += angles[c];
    }
    while (sum >= M_PI) {
      sum -= 2 * M_PI;
    }
    while (sum < -M_PI) {
      sum += 2 * M_PI;
    }
    for (int c = 0; c < 3; c++) {
      abelian[i].matrix[c] = std::complex<double>(cos(angles[c] - sum / 3),
                                                  sin(angles[c] - sum / 3));
    }
  }
  return abelian;
}

std::vector<su3_angles> su3_to_su3_angles(std::vector<su3> &conf) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  std::vector<su3_angles> su3_angles(data_size);
  for (int i = 0; i < data_size; i++) {
    std::vector<double> angles(3);
    double sum = 0;
    for (int c = 0; c < 3; c++) {
      angles[c] =
          atan2(conf[i].matrix(c, c).imag(), conf[i].matrix(c, c).real());
      sum += angles[c];
    }
    while (sum >= M_PI) {
      sum -= 2 * M_PI;
    }
    while (sum < -M_PI) {
      sum += 2 * M_PI;
    }
    for (int c = 0; c < 3; c++) {
      su3_angles[i].matrix[c] = angles[c] - sum / 3;
    }
  }
  return su3_angles;
}

std::vector<abelian> get_abelian(std::vector<su2> &conf) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  std::vector<abelian> conf_abelian(data_size);

  for (int i = 0; i < data_size; i++) {
    conf_abelian[i] = abelian(1, atan2(conf[i].a3, conf[i].a0));
  }

  return conf_abelian;
}

std::vector<abelian>
get_abelian(const Data::LatticeData<DataPatternLexicographical, su2> &conf) {
  DataPatternLexicographical data_pattern(conf.lat_dim);
  int data_size = data_pattern.get_data_size();
  std::vector<abelian> conf_abelian(data_size);
  for (int i = 0; i < data_size; i++) {
    conf_abelian[i] = abelian(1, atan2(conf[i].a3, conf[i].a0));
  }
  return conf_abelian;
}

std::vector<su3> get_offdiagonal(std::vector<su3> &conf) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  std::vector<su3> conf_offdiagonal(data_size);
  su3 A;

  for (int i = 0; i < data_size; i++) {
    for (int c = 0; c < 3; c++) {
      A.matrix(c, c) =
          conf[i].matrix(c, c) / std::sqrt(std::norm(conf[i].matrix(c, c)));
    }
    conf_offdiagonal[i] = conf[i] ^ A;
  }

  return conf_offdiagonal;
}

std::vector<su2> get_offdiagonal(std::vector<su2> &conf) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  std::vector<su2> conf_offdiagonal(data_size);
  su2 A;
  double module;

  for (int i = 0; i < data_size; i++) {
    module = sqrt(conf[i].a0 * conf[i].a0 + conf[i].a3 * conf[i].a3);
    A = su2(conf[i].a0 / module, 0, 0, conf[i].a3 / module);
    conf_offdiagonal[i] = conf[i] ^ A;
  }

  return conf_offdiagonal;
}

std::vector<double>
merge_angles(std::vector<std::vector<double>> &angles,
             DataPatternLexicographical &data_pattern_conf) {
  int data_size = data_pattern_conf.get_lattice_size();
  std::vector<double> angles_merged;
  angles_merged.reserve(4 * data_size);
  for (int i = 0; i < data_size; i++) {
    for (int mu = 0; mu < 4; mu++) {
      angles_merged.push_back(angles[mu][i]);
    }
  }
  return angles_merged;
}

std::vector<double> merge_angles(std::vector<std::vector<double>> &angles) {
  int data_size = x_size * y_size * z_size * t_size;
  std::vector<double> angles_merged;
  angles_merged.reserve(4 * data_size);
  for (int i = 0; i < data_size; i++) {
    for (int mu = 0; mu < 4; mu++) {
      angles_merged.push_back(angles[mu][i]);
    }
  }
  return angles_merged;
}

std::vector<su2> get_monopoless(const std::vector<su2> &conf_su2,
                                const std::vector<double> &angles_monopole) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  std::vector<su2> conf_monopoless(data_size);

  su2 A;

  for (int i = 0; i < conf_su2.size(); i++) {
    A = su2(cos(angles_monopole[i]), 0, 0, sin(angles_monopole[i]));

    conf_monopoless[i] = conf_su2[i] * A;
  }

  return conf_monopoless;
}

void get_monopoless_optimized(std::vector<su2> &conf_su2,
                              const std::vector<double> &angles_monopole) {
  su2 A;
  for (int i = 0; i < conf_su2.size(); i++) {
    A = su2(cos(angles_monopole[i]), 0, 0, sin(angles_monopole[i]));
    conf_su2[i] = conf_su2[i] * A;
  }
}

void get_monopoless_optimized_su3(
    std::vector<su3> &conf_su3,
    const std::vector<std::vector<double>> &angles_monopole) {
  su3 A;
  double module;
  for (int i = 0; i < conf_su3.size(); i++) {
    A = su3();
    for (int j = 0; j < 3; j++) {
      A.matrix(j, j) = std::complex<double>(cos(angles_monopole[j][i]),
                                            sin(angles_monopole[j][i]));
    }
    conf_su3[i] = conf_su3[i] * A;
  }
}

std::vector<su2> get_initial_su2(std::vector<su2> &conf_monopoless,
                                 std::vector<double> &angles_monopole) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  std::vector<su2> conf_su2(data_size);

  su2 A;

  for (int i = 0; i < conf_su2.size(); i++) {
    A = su2(cos(angles_monopole[i]), 0, 0, sin(angles_monopole[i]));

    conf_su2[i] = conf_monopoless[i] ^ A;
  }

  return conf_su2;
}

std::vector<double> read_gauge_Landau(std::string &file_name, int bites_skip) {
  int data_size = x_size * y_size * z_size * t_size;

  std::vector<double> gauge(data_size);

  std::ifstream stream(file_name);
  stream.ignore(bites_skip);
  if (!stream.read((char *)&gauge[0], (data_size) * sizeof(double)))
    std::cout << "read_gauge_Landau error: " << file_name << std::endl;
  return gauge;
}

void apply_gauge_Landau(std::vector<su2> &conf_su2,
                        std::vector<double> &gauge) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  link1 link(x_size, y_size, z_size, t_size);

  su2 A;

  SPACE_ITER_START

  A = su2(cos(gauge[link.place / 4]), 0, 0, sin(gauge[link.place / 4]));

  for (int mu = 0; mu < 4; mu++) {

    conf_su2[link.place + mu] = A * conf_su2[link.place + mu];
  }

  for (int mu = 0; mu < 4; mu++) {

    link.move(mu, -1);

    conf_su2[link.place + mu] = conf_su2[link.place + mu] ^ A;

    link.move(mu, 1);
  }

  SPACE_ITER_END
}

std::vector<double>
read_inverse_laplacian(std::string &file_path,
                       DataPatternLexicographical &data_pattern) {
  std::vector<double> laplace(data_pattern.get_lattice_size());
  std::ifstream stream(file_path);
  if (!stream.read((char *)&laplace[0],
                   data_pattern.get_lattice_size() * sizeof(double)))
    std::cout << "read_inverse_laplacian error: " << file_path << std::endl;
  return laplace;
}

void write_inverse_laplacian(const std::vector<double> &inverse_laplacian,
                             std::string &file_path,
                             DataPatternLexicographical &data_pattern) {
  int data_size = data_pattern.get_lattice_size();
  std::ofstream stream(file_path);
  if (!stream.write((char *)&inverse_laplacian[0], data_size * sizeof(double)))
    std::cout << "write_inverse_laplacian error: " << file_path << std::endl;
  stream.close();
}

double get_monopole_angle(std::vector<std::vector<int>> &monopole_plaket,
                          link1 &link_tmp, std::vector<double> &laplace,
                          int mu) {
  link1 link(x_size, y_size, z_size, t_size);

  std::vector<int> laplace_size = {x_size / 2 + 1, y_size / 2 + 1,
                                   z_size / 2 + 1, t_size / 2 + 1};

  link1 link_laplace(laplace_size[0], laplace_size[1], laplace_size[2],
                     laplace_size[3]);

  double monopole_angle = 0;
  int angle_tmp;

  std::vector<int> laplace_coordinate(4);

  SPACE_ITER_START

  for (int nu = 0; nu < 4; nu++) {
    laplace_coordinate[nu] = abs(link_tmp.coordinate[nu] - link.coordinate[nu]);

    if (laplace_coordinate[nu] > link.lattice_size[nu] / 2)
      laplace_coordinate[nu] = link.lattice_size[nu] - laplace_coordinate[nu];
  }
  link_laplace.go_update(laplace_coordinate[0], laplace_coordinate[1],
                         laplace_coordinate[2], laplace_coordinate[3]);

  for (int nu = 0; nu < 4; nu++) {
    if (nu != mu) {
      int factor;
      int a, b;
      if (mu < nu) {
        a = mu;
        b = nu;
        factor = 1;
      } else {
        a = nu;
        b = mu;
        factor = -1;
      }
      int index = 0;
      for (int i = 0; i <= a; i++) {
        for (int j = i + 1; j < 4; j++) {
          if (i == a && j == b)
            break;
          index++;
        }
      }

      angle_tmp = monopole_plaket[index][link.place / 4];

      link.move(nu, -1);

      monopole_angle += factor * laplace[link_laplace.place / 4] *
                        (angle_tmp - monopole_plaket[index][link.place / 4]);

      link.move(nu, 1);
    }
  }

  SPACE_ITER_END

  return -2 * M_PI * monopole_angle;
}

std::vector<double> make_monopole_angles2(std::vector<double> &angles,
                                          std::vector<double> &laplace) {

  link1 link(x_size, y_size, z_size, t_size);
  std::vector<std::vector<int>> monopole_plaket =
      calculate_monopole_plaket_singular(angles);

  std::vector<double> monopole_angles(4 * x_size * y_size * z_size * t_size);

  SPACE_ITER_START

  for (int mu = 0; mu < 4; mu++) {

    monopole_angles[link.place + mu] =
        get_monopole_angle(monopole_plaket, link, laplace, mu);
  }

  SPACE_ITER_END

  return monopole_angles;
}

void dirac_plaket_difference_nonzero(
    const std::vector<std::vector<int>> &dirac_plaket,
    std::vector<std::vector<int>> &dirac_difference,
    std::vector<std::vector<int>> &dirac_coordinate) {

  link1 link(x_size, y_size, z_size, t_size);

  int angle_tmp;
  int diff;

  SPACE_ITER_START

  int factor;
  int a, b;
  int index;

  for (int mu = 0; mu < 4; mu++) {
    diff = 0;
    for (int nu = 0; nu < 4; nu++) {
      if (nu != mu) {
        if (mu < nu) {
          a = mu;
          b = nu;
          factor = 1;
        } else {
          a = nu;
          b = mu;
          factor = -1;
        }
        int index = 0;
        for (int i = 0; i <= a; i++) {
          for (int j = i + 1; j < 4; j++) {
            if (i == a && j == b)
              break;
            index++;
          }
        }

        angle_tmp = dirac_plaket[index][link.place / 4];

        link.move(nu, -1);

        diff += factor * (angle_tmp - dirac_plaket[index][link.place / 4]);

        link.move(nu, 1);
      }
    }

    if (diff != 0) {
      dirac_difference[mu].push_back(diff);
      for (int nu = 0; nu < 4; nu++) {
        dirac_coordinate[mu].push_back(link.coordinate[nu]);
      }
    }
  }

  SPACE_ITER_END
}

void dirac_plaket_difference_nonzero(
    const std::vector<std::vector<int>> &dirac_plaket,
    std::vector<std::vector<int>> &dirac_difference,
    std::vector<std::vector<int>> &dirac_coordinate,
    DataPatternLexicographical &data_pattern) {
  int angle_tmp;
  int diff;

  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};

          int factor;
          int a, b;
          int index;

          for (int mu = 0; mu < 4; mu++) {
            diff = 0;
            for (int nu = 0; nu < 4; nu++) {
              if (nu != mu) {
                if (mu < nu) {
                  a = mu;
                  b = nu;
                  factor = 1;
                } else {
                  a = nu;
                  b = mu;
                  factor = -1;
                }
                int index = 0;
                for (int i = 0; i <= a; i++) {
                  for (int j = i + 1; j < 4; j++) {
                    if (i == a && j == b)
                      break;
                    index++;
                  }
                }

                angle_tmp = dirac_plaket[index][data_pattern.get_index_site()];
                data_pattern.move_backward(1, nu);
                diff += factor *
                        (angle_tmp -
                         dirac_plaket[index][data_pattern.get_index_site()]);
                data_pattern.move_forward(1, nu);
              }
            }

            if (diff != 0) {
              dirac_difference[mu].push_back(diff);
              for (int nu = 0; nu < 4; nu++) {
                dirac_coordinate[mu].push_back(data_pattern.lat_coord[nu]);
              }
            }
          }
        }
      }
    }
  }
}

void decomposition_step(
    const std::vector<std::vector<int>> &monopole_difference,
    const std::vector<std::vector<int>> &monopole_coordinate,
    const std::vector<double> &laplace, std::vector<double> &angles_decomposed,
    int i, int mu) {
  std::vector<int> laplace_size = {x_size / 2 + 1, y_size / 2 + 1,
                                   z_size / 2 + 1, t_size / 2 + 1};
  std::vector<int> laplace_coordinate(4);

  int laplace_place;
  int angles_place;

  angles_place = 0;
  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {

          laplace_coordinate[0] = abs(x - monopole_coordinate[mu][i * 4]);
          if (laplace_coordinate[0] > x_size / 2)
            laplace_coordinate[0] = x_size - laplace_coordinate[0];

          laplace_coordinate[1] = abs(y - monopole_coordinate[mu][i * 4 + 1]);
          if (laplace_coordinate[1] > y_size / 2)
            laplace_coordinate[1] = y_size - laplace_coordinate[1];

          laplace_coordinate[2] = abs(z - monopole_coordinate[mu][i * 4 + 2]);
          if (laplace_coordinate[2] > z_size / 2)
            laplace_coordinate[2] = z_size - laplace_coordinate[2];

          laplace_coordinate[3] = abs(t - monopole_coordinate[mu][i * 4 + 3]);
          if (laplace_coordinate[3] > t_size / 2)
            laplace_coordinate[3] = t_size - laplace_coordinate[3];

          laplace_place =
              (laplace_coordinate[3]) * laplace_size[0] * laplace_size[1] *
                  laplace_size[2] +
              (laplace_coordinate[2]) * laplace_size[0] * laplace_size[1] +
              (laplace_coordinate[1]) * laplace_size[0] + laplace_coordinate[0];

          angles_decomposed[angles_place * 4 + mu] +=
              laplace[laplace_place] * monopole_difference[mu][i];

          angles_place++;
        }
      }
    }
  }
}

std::vector<double>
make_monopole_angles1(const std::vector<std::vector<int>> &dirac_plakets,
                      const std::vector<double> &laplace) {
  link1 link(x_size, y_size, z_size, t_size);

  int data_size = 4 * x_size * y_size * z_size * t_size;

  std::vector<double> angles_decomposed(data_size);

  std::vector<std::vector<int>> dirac_difference(4, std::vector<int>());
  std::vector<std::vector<int>> dirac_coordinate(4, std::vector<int>());

  dirac_plaket_difference_nonzero(dirac_plakets, dirac_difference,
                                  dirac_coordinate);

  std::vector<int> laplace_size = {x_size / 2 + 1, y_size / 2 + 1,
                                   z_size / 2 + 1, t_size / 2 + 1};

  double dirac_angle = 0;
  int angle_tmp;

  std::vector<int> laplace_coordinate(4);
  int laplace_place;
  int angles_place;

  for (int mu = 0; mu < 4; mu++) {
    for (int i = 0; i < dirac_difference[mu].size(); i++) {

      decomposition_step(dirac_difference, dirac_coordinate, laplace,
                         angles_decomposed, i, mu);
    }
  }

  for (int i = 0; i < angles_decomposed.size(); i++) {
    angles_decomposed[i] = -2 * M_PI * angles_decomposed[i];
  }
  return angles_decomposed;
}

void decomposition_step3(
    const std::vector<std::vector<int>> &monopole_difference,
    const std::vector<std::vector<int>> &monopole_coordinate,
    const std::vector<double> &laplace,
    std::vector<std::vector<double>> &angles_decomposed, int i, int mu) {
  std::vector<int> laplace_size = {x_size / 2 + 1, y_size / 2 + 1,
                                   z_size / 2 + 1, t_size / 2 + 1};
  std::vector<int> laplace_coordinate(4);

  int laplace_place;
  int angles_place;

  angles_place = 0;
  laplace_place = 0;
  for (int t = 0; t < 1; t++) {
    laplace_coordinate[3] = abs(t - monopole_coordinate[mu][i * 4 + 3]);
    if (laplace_coordinate[3] > t_size / 2)
      laplace_coordinate[3] = t_size - laplace_coordinate[3];
    for (int z = 0; z < 1; z++) {
      laplace_coordinate[2] = abs(z - monopole_coordinate[mu][i * 4 + 2]);
      if (laplace_coordinate[2] > z_size / 2)
        laplace_coordinate[2] = z_size - laplace_coordinate[2];
      for (int y = 0; y < y_size; y++) {
        laplace_coordinate[1] = abs(y - monopole_coordinate[mu][i * 4 + 1]);
        if (laplace_coordinate[1] > y_size / 2)
          laplace_coordinate[1] = y_size - laplace_coordinate[1];
        for (int x = 0; x < x_size; x++) {

          laplace_coordinate[0] = abs(x - monopole_coordinate[mu][i * 4]);
          if (laplace_coordinate[0] > x_size / 2)
            laplace_coordinate[0] = x_size - laplace_coordinate[0];

          laplace_place =
              (laplace_coordinate[3]) * laplace_size[0] * laplace_size[1] *
                  laplace_size[2] +
              (laplace_coordinate[2]) * laplace_size[0] * laplace_size[1] +
              (laplace_coordinate[1]) * laplace_size[0] + laplace_coordinate[0];

          angles_decomposed[mu][angles_place] +=
              laplace[laplace_place] * monopole_difference[mu][i];

          angles_place++;
        }
      }
    }
  }
}

std::vector<std::vector<double>>
make_monopole_angles3(const std::vector<std::vector<int>> &dirac_plakets,
                      const std::vector<double> &laplace) {
  link1 link(x_size, y_size, z_size, t_size);

  int data_size = x_size * y_size * z_size * t_size;

  std::vector<std::vector<double>> angles_decomposed(
      4, std::vector<double>(data_size));

  std::vector<std::vector<int>> dirac_difference(4, std::vector<int>());
  std::vector<std::vector<int>> dirac_coordinate(4, std::vector<int>());

  dirac_plaket_difference_nonzero(dirac_plakets, dirac_difference,
                                  dirac_coordinate);

  std::vector<int> laplace_size = {x_size / 2 + 1, y_size / 2 + 1,
                                   z_size / 2 + 1, t_size / 2 + 1};

  double dirac_angle = 0;
  int angle_tmp;

  std::vector<int> laplace_coordinate(4);
  int laplace_place;
  int angles_place;

  for (int mu = 0; mu < 4; mu++) {
    for (int i = 0; i < dirac_difference[mu].size(); i++) {

      decomposition_step3(dirac_difference, dirac_coordinate, laplace,
                          angles_decomposed, i, mu);
    }
  }

  for (int mu = 0; mu < 4; mu++) {
    for (int i = 0; i < angles_decomposed[mu].size(); i++) {
      angles_decomposed[mu][i] = -2 * M_PI * angles_decomposed[mu][i];
    }
  }
  return angles_decomposed;
}

void decomposition_step(int monopole_difference,
                        const std::vector<int> &monopole_coordinate,
                        const std::vector<double> &laplace,
                        std::vector<double> &angles_decomposed) {
  std::vector<int> laplace_size = {x_size / 2 + 1, y_size / 2 + 1,
                                   z_size / 2 + 1, t_size / 2 + 1};
  std::vector<int> laplace_shift = {
      1, laplace_size[0], laplace_size[0] * laplace_size[1],
      laplace_size[0] * laplace_size[1] * laplace_size[2]};
  std::vector<int> laplace_coordinate(4);

  int laplace_place;
  int angles_place;

  angles_place = 0;

  for (int t = 0; t < t_size; t++) {

    laplace_coordinate[3] = abs(t - monopole_coordinate[3]);
    if (laplace_coordinate[3] > t_size / 2)
      laplace_coordinate[3] = t_size - laplace_coordinate[3];

    for (int z = 0; z < z_size; z++) {

      laplace_coordinate[2] = abs(z - monopole_coordinate[2]);
      if (laplace_coordinate[2] > z_size / 2)
        laplace_coordinate[2] = z_size - laplace_coordinate[2];

      for (int y = 0; y < y_size; y++) {

        laplace_coordinate[1] = abs(y - monopole_coordinate[1]);
        if (laplace_coordinate[1] > y_size / 2)
          laplace_coordinate[1] = y_size - laplace_coordinate[1];

        laplace_place = laplace_coordinate[3] * laplace_shift[3] +
                        laplace_coordinate[2] * laplace_shift[2] +
                        laplace_coordinate[1] * laplace_shift[1];

        if (monopole_coordinate[0] > x_size / 2) {

          for (int j = x_size - monopole_coordinate[0]; j <= x_size / 2; j++) {

            angles_decomposed[angles_place] +=
                laplace[laplace_place + j] * monopole_difference;

            angles_place++;
          }

          for (int j = x_size / 2 - 1; j > 0; j--) {

            angles_decomposed[angles_place] +=
                laplace[laplace_place + j] * monopole_difference;

            angles_place++;
          }

          for (int j = 0; j < x_size - monopole_coordinate[0]; j++) {

            angles_decomposed[angles_place] +=
                laplace[laplace_place + j] * monopole_difference;

            angles_place++;
          }
        } else {

          for (int j = monopole_coordinate[0]; j >= 0; j--) {

            angles_decomposed[angles_place] +=
                laplace[laplace_place + j] * monopole_difference;

            angles_place++;
          }

          for (int j = 1; j < x_size / 2; j++) {

            angles_decomposed[angles_place] +=
                laplace[laplace_place + j] * monopole_difference;

            angles_place++;
          }

          for (int j = x_size / 2; j > monopole_coordinate[0]; j--) {

            angles_decomposed[angles_place] +=
                laplace[laplace_place + j] * monopole_difference;

            angles_place++;
          }
        }
      }
    }
  }
}

std::vector<double>
make_monopole_angles(std::vector<std::vector<int>> &dirac_plakets,
                     const std::vector<double> &laplace) {

  link1 link(x_size, y_size, z_size, t_size);

  int data_size = x_size * y_size * z_size * t_size;

  std::vector<std::vector<int>> dirac_difference(4, std::vector<int>());
  std::vector<std::vector<int>> dirac_coordinate(4, std::vector<int>());

  dirac_plaket_difference_nonzero(dirac_plakets, dirac_difference,
                                  dirac_coordinate);

  for (int mu = 0; mu < dirac_plakets.size(); mu++) {
    dirac_plakets[mu].clear();
    dirac_plakets[mu].shrink_to_fit();
  }

  std::vector<std::vector<double>> angles_decomposed(
      4, std::vector<double>(data_size));

  std::vector<int> coordinate(4);

  for (int mu = 0; mu < 4; mu++) {

    for (int i = 0; i < dirac_difference[mu].size(); i++) {

      coordinate = {
          dirac_coordinate[mu][4 * i], dirac_coordinate[mu][4 * i + 1],
          dirac_coordinate[mu][4 * i + 2], dirac_coordinate[mu][4 * i + 3]};

      decomposition_step(dirac_difference[mu][i], coordinate, laplace,
                         angles_decomposed[mu]);
    }
  }

  for (int mu = 0; mu < 4; mu++) {
    for (int i = 0; i < angles_decomposed[mu].size(); i++) {
      angles_decomposed[mu][i] = -2 * M_PI * angles_decomposed[mu][i];
    }
  }
  return merge_angles(angles_decomposed);
}

void decomposition_step_parallel(int monopole_difference,
                                 const std::vector<int> &monopole_coordinate,
                                 const std::vector<double> &laplace,
                                 std::vector<double> &angles_decomposed) {
  std::vector<int> laplace_size = {x_size / 2 + 1, y_size / 2 + 1,
                                   z_size / 2 + 1, t_size / 2 + 1};
  std::vector<int> laplace_shift = {
      laplace_size[0], laplace_size[0] * laplace_size[1],
      laplace_size[0] * laplace_size[1] * laplace_size[2]};
  std::vector<int> angles_shift = {x_size, y_size * x_size,
                                   z_size * y_size * x_size};
  std::vector<int> laplace_coordinate(3);

  int laplace_place;
  int angles_place;

  angles_place = 0;

#pragma omp parallel for collapse(2) firstprivate(                             \
        angles_place, laplace_place, laplace_coordinate, laplace_shift,        \
            angles_shift, x_size, monopole_difference, monopole_coordinate)
  for (int t = 0; t < t_size; t++) {

    for (int z = 0; z < z_size; z++) {

      laplace_coordinate[2] = abs(t - monopole_coordinate[3]);
      if (laplace_coordinate[2] > t_size / 2)
        laplace_coordinate[2] = t_size - laplace_coordinate[2];

      laplace_coordinate[1] = abs(z - monopole_coordinate[2]);
      if (laplace_coordinate[1] > z_size / 2)
        laplace_coordinate[1] = z_size - laplace_coordinate[1];

      for (int y = 0; y < y_size; y++) {

        laplace_coordinate[0] = abs(y - monopole_coordinate[1]);
        if (laplace_coordinate[0] > y_size / 2)
          laplace_coordinate[0] = y_size - laplace_coordinate[0];

        laplace_place = laplace_coordinate[2] * laplace_shift[2] +
                        laplace_coordinate[1] * laplace_shift[1] +
                        laplace_coordinate[0] * laplace_shift[0];

        angles_place =
            angles_shift[2] * t + angles_shift[1] * z + angles_shift[0] * y;

        if (monopole_coordinate[0] > x_size / 2) {

          for (int j = x_size - monopole_coordinate[0]; j <= x_size / 2; j++) {

            angles_decomposed[angles_place] +=
                laplace[laplace_place + j] * monopole_difference;

            angles_place++;
          }

          for (int j = x_size / 2 - 1; j > 0; j--) {

            angles_decomposed[angles_place] +=
                laplace[laplace_place + j] * monopole_difference;

            angles_place++;
          }

          for (int j = 0; j < x_size - monopole_coordinate[0]; j++) {

            angles_decomposed[angles_place] +=
                laplace[laplace_place + j] * monopole_difference;

            angles_place++;
          }
        } else {

          for (int j = monopole_coordinate[0]; j >= 0; j--) {

            angles_decomposed[angles_place] +=
                laplace[laplace_place + j] * monopole_difference;

            angles_place++;
          }

          for (int j = 1; j < x_size / 2; j++) {

            angles_decomposed[angles_place] +=
                laplace[laplace_place + j] * monopole_difference;

            angles_place++;
          }

          for (int j = x_size / 2; j > monopole_coordinate[0]; j--) {

            angles_decomposed[angles_place] +=
                laplace[laplace_place + j] * monopole_difference;

            angles_place++;
          }
        }
      }
    }
  }
}

void decomposition_step_parallel3_simple_positive(
    const std::vector<int> &monopole_coordinate,
    const std::vector<double> &laplace,
    std::vector<double> &angles_decomposed) {
  std::vector<int> laplace_size = {x_size / 2 + 1, y_size / 2 + 1,
                                   z_size / 2 + 1, t_size / 2 + 1};
  std::vector<int> laplace_shift = {
      laplace_size[0], laplace_size[0] * laplace_size[1],
      laplace_size[0] * laplace_size[1] * laplace_size[2]};
  std::vector<int> angles_shift = {x_size, y_size * x_size,
                                   z_size * y_size * x_size};
  std::vector<int> laplace_coordinate(3);

  int laplace_place;
  int angles_place;

  angles_place = 0;

#pragma omp parallel for collapse(2)                                           \
    firstprivate(angles_place, laplace_place, laplace_coordinate,              \
                     laplace_shift, angles_shift, x_size, monopole_coordinate)
  for (int t = 0; t < t_size; t++) {

    for (int z = 0; z < z_size; z++) {

      laplace_coordinate[2] = abs(t - monopole_coordinate[3]);
      if (laplace_coordinate[2] > t_size / 2)
        laplace_coordinate[2] = t_size - laplace_coordinate[2];

      laplace_coordinate[1] = abs(z - monopole_coordinate[2]);
      if (laplace_coordinate[1] > z_size / 2)
        laplace_coordinate[1] = z_size - laplace_coordinate[1];

      for (int y = 0; y < y_size; y++) {

        laplace_coordinate[0] = abs(y - monopole_coordinate[1]);
        if (laplace_coordinate[0] > y_size / 2)
          laplace_coordinate[0] = y_size - laplace_coordinate[0];

        laplace_place = laplace_coordinate[2] * laplace_shift[2] +
                        laplace_coordinate[1] * laplace_shift[1] +
                        laplace_coordinate[0] * laplace_shift[0];

        angles_place =
            angles_shift[2] * t + angles_shift[1] * z + angles_shift[0] * y;

        if (monopole_coordinate[0] > x_size / 2) {

          for (int j = x_size - monopole_coordinate[0]; j <= x_size / 2; j++) {

            angles_decomposed[angles_place] += laplace[laplace_place + j];

            angles_place++;
          }

          for (int j = x_size / 2 - 1; j > 0; j--) {

            angles_decomposed[angles_place] += laplace[laplace_place + j];

            angles_place++;
          }

          for (int j = 0; j < x_size - monopole_coordinate[0]; j++) {

            angles_decomposed[angles_place] += laplace[laplace_place + j];

            angles_place++;
          }
        } else {

          for (int j = monopole_coordinate[0]; j >= 0; j--) {

            angles_decomposed[angles_place] += laplace[laplace_place + j];

            angles_place++;
          }

          for (int j = 1; j < x_size / 2; j++) {

            angles_decomposed[angles_place] += laplace[laplace_place + j];

            angles_place++;
          }

          for (int j = x_size / 2; j > monopole_coordinate[0]; j--) {

            angles_decomposed[angles_place] += laplace[laplace_place + j];

            angles_place++;
          }
        }
      }
    }
  }
}

void decomposition_step_parallel3_simple_negative(
    const std::vector<int> &monopole_coordinate,
    const std::vector<double> &laplace,
    std::vector<double> &angles_decomposed) {
  std::vector<int> laplace_size = {x_size / 2 + 1, y_size / 2 + 1,
                                   z_size / 2 + 1, t_size / 2 + 1};
  std::vector<int> laplace_shift = {
      laplace_size[0], laplace_size[0] * laplace_size[1],
      laplace_size[0] * laplace_size[1] * laplace_size[2]};
  std::vector<int> angles_shift = {x_size, y_size * x_size,
                                   z_size * y_size * x_size};
  std::vector<int> laplace_coordinate(3);

  int laplace_place;
  int angles_place;

  angles_place = 0;

#pragma omp parallel for collapse(2)                                           \
    firstprivate(angles_place, laplace_place, laplace_coordinate,              \
                     laplace_shift, angles_shift, x_size, monopole_coordinate)
  for (int t = 0; t < t_size; t++) {

    for (int z = 0; z < z_size; z++) {

      laplace_coordinate[2] = abs(t - monopole_coordinate[3]);
      if (laplace_coordinate[2] > t_size / 2)
        laplace_coordinate[2] = t_size - laplace_coordinate[2];

      laplace_coordinate[1] = abs(z - monopole_coordinate[2]);
      if (laplace_coordinate[1] > z_size / 2)
        laplace_coordinate[1] = z_size - laplace_coordinate[1];

      for (int y = 0; y < y_size; y++) {

        laplace_coordinate[0] = abs(y - monopole_coordinate[1]);
        if (laplace_coordinate[0] > y_size / 2)
          laplace_coordinate[0] = y_size - laplace_coordinate[0];

        laplace_place = laplace_coordinate[2] * laplace_shift[2] +
                        laplace_coordinate[1] * laplace_shift[1] +
                        laplace_coordinate[0] * laplace_shift[0];

        angles_place =
            angles_shift[2] * t + angles_shift[1] * z + angles_shift[0] * y;

        if (monopole_coordinate[0] > x_size / 2) {

          for (int j = x_size - monopole_coordinate[0]; j <= x_size / 2; j++) {

            angles_decomposed[angles_place] -= laplace[laplace_place + j];

            angles_place++;
          }

          for (int j = x_size / 2 - 1; j > 0; j--) {

            angles_decomposed[angles_place] -= laplace[laplace_place + j];

            angles_place++;
          }

          for (int j = 0; j < x_size - monopole_coordinate[0]; j++) {

            angles_decomposed[angles_place] -= laplace[laplace_place + j];

            angles_place++;
          }
        } else {

          for (int j = monopole_coordinate[0]; j >= 0; j--) {

            angles_decomposed[angles_place] -= laplace[laplace_place + j];

            angles_place++;
          }

          for (int j = 1; j < x_size / 2; j++) {

            angles_decomposed[angles_place] -= laplace[laplace_place + j];

            angles_place++;
          }

          for (int j = x_size / 2; j > monopole_coordinate[0]; j--) {

            angles_decomposed[angles_place] -= laplace[laplace_place + j];

            angles_place++;
          }
        }
      }
    }
  }
}

void decomposition_step_parallel1(int monopole_difference,
                                  const std::vector<int> &monopole_coordinate,
                                  const std::vector<double> &laplace,
                                  std::vector<double> &angles_decomposed) {
  std::vector<int> laplace_size = {x_size / 2 + 1, y_size / 2 + 1,
                                   z_size / 2 + 1, t_size / 2 + 1};
  std::vector<int> laplace_shift = {
      laplace_size[0], laplace_size[0] * laplace_size[1],
      laplace_size[0] * laplace_size[1] * laplace_size[2]};
  std::vector<int> angles_shift = {x_size, y_size * x_size,
                                   z_size * y_size * x_size};
  std::vector<int> laplace_coordinate(3);

  int laplace_place;
  int angles_place;

  angles_place = 0;

#pragma omp parallel for collapse(2) firstprivate(                             \
        angles_place, laplace_place, laplace_coordinate, laplace_shift,        \
            angles_shift, x_size, monopole_difference, monopole_coordinate)
  for (int t = 0; t < t_size; t++) {

    for (int z = 0; z < z_size; z++) {

      laplace_coordinate[2] = abs(t - monopole_coordinate[3]);
      if (laplace_coordinate[2] > t_size / 2)
        laplace_coordinate[2] = t_size - laplace_coordinate[2];

      laplace_coordinate[1] = abs(z - monopole_coordinate[2]);
      if (laplace_coordinate[1] > z_size / 2)
        laplace_coordinate[1] = z_size - laplace_coordinate[1];

      for (int y = 0; y < y_size; y++) {

        laplace_coordinate[0] = abs(y - monopole_coordinate[1]);
        if (laplace_coordinate[0] > y_size / 2)
          laplace_coordinate[0] = y_size - laplace_coordinate[0];

        laplace_place = laplace_coordinate[2] * laplace_shift[2] +
                        laplace_coordinate[1] * laplace_shift[1] +
                        laplace_coordinate[0] * laplace_shift[0];

        angles_place =
            angles_shift[2] * t + angles_shift[1] * z + angles_shift[0] * y;

        for (int j = x_size - monopole_coordinate[0]; j <= x_size / 2; j++) {

          angles_decomposed[angles_place] +=
              laplace[laplace_place + j] * monopole_difference;

          angles_place++;
        }

        for (int j = x_size / 2 - 1; j > 0; j--) {

          angles_decomposed[angles_place] +=
              laplace[laplace_place + j] * monopole_difference;

          angles_place++;
        }

        for (int j = 0; j < x_size - monopole_coordinate[0]; j++) {

          angles_decomposed[angles_place] +=
              laplace[laplace_place + j] * monopole_difference;

          angles_place++;
        }
      }
    }
  }
}

void decomposition_step_parallel2(int monopole_difference,
                                  const std::vector<int> &monopole_coordinate,
                                  const std::vector<double> &laplace,
                                  std::vector<double> &angles_decomposed) {
  std::vector<int> laplace_size = {x_size / 2 + 1, y_size / 2 + 1,
                                   z_size / 2 + 1, t_size / 2 + 1};
  std::vector<int> laplace_shift = {
      laplace_size[0], laplace_size[0] * laplace_size[1],
      laplace_size[0] * laplace_size[1] * laplace_size[2]};
  std::vector<int> angles_shift = {x_size, y_size * x_size,
                                   z_size * y_size * x_size};
  std::vector<int> laplace_coordinate(3);

  int laplace_place;
  int angles_place;

  angles_place = 0;

#pragma omp parallel for collapse(2) firstprivate(                             \
        angles_place, laplace_place, laplace_coordinate, laplace_shift,        \
            angles_shift, x_size, monopole_difference, monopole_coordinate)
  for (int t = 0; t < t_size; t++) {

    for (int z = 0; z < z_size; z++) {

      laplace_coordinate[2] = abs(t - monopole_coordinate[3]);
      if (laplace_coordinate[2] > t_size / 2)
        laplace_coordinate[2] = t_size - laplace_coordinate[2];

      laplace_coordinate[1] = abs(z - monopole_coordinate[2]);
      if (laplace_coordinate[1] > z_size / 2)
        laplace_coordinate[1] = z_size - laplace_coordinate[1];

      for (int y = 0; y < y_size; y++) {

        laplace_coordinate[0] = abs(y - monopole_coordinate[1]);
        if (laplace_coordinate[0] > y_size / 2)
          laplace_coordinate[0] = y_size - laplace_coordinate[0];

        laplace_place = laplace_coordinate[2] * laplace_shift[2] +
                        laplace_coordinate[1] * laplace_shift[1] +
                        laplace_coordinate[0] * laplace_shift[0];

        angles_place =
            angles_shift[2] * t + angles_shift[1] * z + angles_shift[0] * y;

        for (int j = monopole_coordinate[0]; j >= 0; j--) {

          angles_decomposed[angles_place] +=
              laplace[laplace_place + j] * monopole_difference;

          angles_place++;
        }

        for (int j = 1; j < x_size / 2; j++) {

          angles_decomposed[angles_place] +=
              laplace[laplace_place + j] * monopole_difference;

          angles_place++;
        }

        for (int j = x_size / 2; j > monopole_coordinate[0]; j--) {

          angles_decomposed[angles_place] +=
              laplace[laplace_place + j] * monopole_difference;

          angles_place++;
        }
      }
    }
  }
}

std::vector<double>
make_monopole_angles_parallel(std::vector<std::vector<int>> &dirac_plakets,
                              const std::vector<double> &laplace) {

  link1 link(x_size, y_size, z_size, t_size);

  int data_size = x_size * y_size * z_size * t_size;

  std::vector<std::vector<int>> dirac_difference(4, std::vector<int>());
  std::vector<std::vector<int>> dirac_coordinate(4, std::vector<int>());

  dirac_plaket_difference_nonzero(dirac_plakets, dirac_difference,
                                  dirac_coordinate);

  for (int mu = 0; mu < dirac_plakets.size(); mu++) {
    dirac_plakets[mu].clear();
    dirac_plakets[mu].shrink_to_fit();
  }

  std::vector<std::vector<double>> angles_decomposed(
      4, std::vector<double>(data_size));

  std::vector<int> coordinate(4);

  for (int mu = 0; mu < 4; mu++) {

    for (int i = 0; i < dirac_difference[mu].size(); i++) {

      coordinate = {
          dirac_coordinate[mu][4 * i], dirac_coordinate[mu][4 * i + 1],
          dirac_coordinate[mu][4 * i + 2], dirac_coordinate[mu][4 * i + 3]};

      if (dirac_coordinate[mu][i * 4] > x_size / 2) {
        decomposition_step_parallel1(dirac_difference[mu][i], coordinate,
                                     laplace, angles_decomposed[mu]);
      } else {
        decomposition_step_parallel2(dirac_difference[mu][i], coordinate,
                                     laplace, angles_decomposed[mu]);
      }
    }
  }

  for (int mu = 0; mu < 4; mu++) {
    for (int i = 0; i < angles_decomposed[mu].size(); i++) {
      angles_decomposed[mu][i] = -2 * M_PI * angles_decomposed[mu][i];
    }
  }
  return merge_angles(angles_decomposed);
}

void decomposition_step_parallel1(
    int monopole_difference, const std::array<int, 4> &monopole_coordinate,
    const std::vector<double> &laplace, std::vector<double> &angles_decomposed,
    DataPatternLexicographical &data_pattern_conf,
    DataPatternLexicographical &data_pattern_laplacian) {
  std::vector<int> laplace_size = {
      data_pattern_laplacian.lat_dim[0], data_pattern_laplacian.lat_dim[1],
      data_pattern_laplacian.lat_dim[2], data_pattern_laplacian.lat_dim[3]};
  std::vector<int> laplace_shift = {
      laplace_size[0], laplace_size[0] * laplace_size[1],
      laplace_size[0] * laplace_size[1] * laplace_size[2]};
  std::vector<int> angles_shift = {
      data_pattern_conf.lat_dim[0],
      data_pattern_conf.lat_dim[1] * data_pattern_conf.lat_dim[0],
      data_pattern_conf.lat_dim[2] * data_pattern_conf.lat_dim[1] *
          data_pattern_conf.lat_dim[0]};
  std::vector<int> laplace_coordinate(3);

  int laplace_place;
  int angles_place;

  angles_place = 0;

#pragma omp parallel for collapse(2) firstprivate(                             \
        angles_place, laplace_place, laplace_coordinate, laplace_shift,        \
            angles_shift, data_pattern_conf, data_pattern_laplacian,           \
            monopole_difference, monopole_coordinate)
  for (int t = 0; t < data_pattern_conf.lat_dim[3]; t++) {

    for (int z = 0; z < data_pattern_conf.lat_dim[2]; z++) {

      laplace_coordinate[2] = abs(t - monopole_coordinate[3]);
      if (laplace_coordinate[2] > data_pattern_conf.lat_dim[3] / 2)
        laplace_coordinate[2] =
            data_pattern_conf.lat_dim[3] - laplace_coordinate[2];

      laplace_coordinate[1] = abs(z - monopole_coordinate[2]);
      if (laplace_coordinate[1] > data_pattern_conf.lat_dim[2] / 2)
        laplace_coordinate[1] =
            data_pattern_conf.lat_dim[2] - laplace_coordinate[1];

      for (int y = 0; y < data_pattern_conf.lat_dim[1]; y++) {

        laplace_coordinate[0] = abs(y - monopole_coordinate[1]);
        if (laplace_coordinate[0] > data_pattern_conf.lat_dim[1] / 2)
          laplace_coordinate[0] =
              data_pattern_conf.lat_dim[1] - laplace_coordinate[0];

        laplace_place = laplace_coordinate[2] * laplace_shift[2] +
                        laplace_coordinate[1] * laplace_shift[1] +
                        laplace_coordinate[0] * laplace_shift[0];

        angles_place =
            angles_shift[2] * t + angles_shift[1] * z + angles_shift[0] * y;

        for (int j = data_pattern_conf.lat_dim[0] - monopole_coordinate[0];
             j <= data_pattern_conf.lat_dim[0] / 2; j++) {

          angles_decomposed[angles_place] +=
              laplace[laplace_place + j] * monopole_difference;

          angles_place++;
        }

        for (int j = data_pattern_conf.lat_dim[0] / 2 - 1; j > 0; j--) {

          angles_decomposed[angles_place] +=
              laplace[laplace_place + j] * monopole_difference;

          angles_place++;
        }

        for (int j = 0;
             j < data_pattern_conf.lat_dim[0] - monopole_coordinate[0]; j++) {

          angles_decomposed[angles_place] +=
              laplace[laplace_place + j] * monopole_difference;

          angles_place++;
        }
      }
    }
  }
}

void decomposition_step_parallel2(
    int monopole_difference, const std::array<int, 4> &monopole_coordinate,
    const std::vector<double> &laplace, std::vector<double> &angles_decomposed,
    DataPatternLexicographical &data_pattern_conf,
    DataPatternLexicographical &data_pattern_laplacian) {
  std::vector<int> laplace_size = {
      data_pattern_laplacian.lat_dim[0], data_pattern_laplacian.lat_dim[1],
      data_pattern_laplacian.lat_dim[2], data_pattern_laplacian.lat_dim[3]};
  std::vector<int> laplace_shift = {
      laplace_size[0], laplace_size[0] * laplace_size[1],
      laplace_size[0] * laplace_size[1] * laplace_size[2]};
  std::vector<int> angles_shift = {
      data_pattern_conf.lat_dim[0],
      data_pattern_conf.lat_dim[1] * data_pattern_conf.lat_dim[0],
      data_pattern_conf.lat_dim[2] * data_pattern_conf.lat_dim[1] *
          data_pattern_conf.lat_dim[0]};
  std::vector<int> laplace_coordinate(3);

  int laplace_place;
  int angles_place;

  angles_place = 0;

#pragma omp parallel for collapse(2) firstprivate(                             \
        angles_place, laplace_place, laplace_coordinate, laplace_shift,        \
            angles_shift, data_pattern_conf, data_pattern_laplacian,           \
            monopole_difference, monopole_coordinate)
  for (int t = 0; t < data_pattern_conf.lat_dim[3]; t++) {

    for (int z = 0; z < data_pattern_conf.lat_dim[2]; z++) {

      laplace_coordinate[2] = abs(t - monopole_coordinate[3]);
      if (laplace_coordinate[2] > data_pattern_conf.lat_dim[3] / 2)
        laplace_coordinate[2] =
            data_pattern_conf.lat_dim[3] - laplace_coordinate[2];

      laplace_coordinate[1] = abs(z - monopole_coordinate[2]);
      if (laplace_coordinate[1] > data_pattern_conf.lat_dim[2] / 2)
        laplace_coordinate[1] =
            data_pattern_conf.lat_dim[2] - laplace_coordinate[1];

      for (int y = 0; y < data_pattern_conf.lat_dim[1]; y++) {

        laplace_coordinate[0] = abs(y - monopole_coordinate[1]);
        if (laplace_coordinate[0] > data_pattern_conf.lat_dim[1] / 2)
          laplace_coordinate[0] =
              data_pattern_conf.lat_dim[1] - laplace_coordinate[0];

        laplace_place = laplace_coordinate[2] * laplace_shift[2] +
                        laplace_coordinate[1] * laplace_shift[1] +
                        laplace_coordinate[0] * laplace_shift[0];

        angles_place =
            angles_shift[2] * t + angles_shift[1] * z + angles_shift[0] * y;

        for (int j = monopole_coordinate[0]; j >= 0; j--) {

          angles_decomposed[angles_place] +=
              laplace[laplace_place + j] * monopole_difference;

          angles_place++;
        }

        for (int j = 1; j < data_pattern_conf.lat_dim[0] / 2; j++) {

          angles_decomposed[angles_place] +=
              laplace[laplace_place + j] * monopole_difference;

          angles_place++;
        }

        for (int j = data_pattern_conf.lat_dim[0] / 2;
             j > monopole_coordinate[0]; j--) {

          angles_decomposed[angles_place] +=
              laplace[laplace_place + j] * monopole_difference;

          angles_place++;
        }
      }
    }
  }
}

std::vector<double>
make_monopole_angles(std::vector<std::vector<int>> &dirac_plakets,
                     const std::vector<double> &laplace,
                     DataPatternLexicographical &data_pattern_conf,
                     DataPatternLexicographical &data_pattern_laplacian) {
  std::vector<std::vector<int>> dirac_difference(4, std::vector<int>());
  std::vector<std::vector<int>> dirac_coordinate(4, std::vector<int>());
  dirac_plaket_difference_nonzero(dirac_plakets, dirac_difference,
                                  dirac_coordinate, data_pattern_conf);
  for (int mu = 0; mu < dirac_plakets.size(); mu++) {
    dirac_plakets[mu].clear();
    dirac_plakets[mu].shrink_to_fit();
  }
  std::vector<std::vector<double>> angles_decomposed(
      4, std::vector<double>(data_pattern_conf.get_lattice_size()));
  std::array<int, 4> coordinate;
  for (int mu = 0; mu < 4; mu++) {
    for (int i = 0; i < dirac_difference[mu].size(); i++) {
      coordinate = {
          dirac_coordinate[mu][4 * i], dirac_coordinate[mu][4 * i + 1],
          dirac_coordinate[mu][4 * i + 2], dirac_coordinate[mu][4 * i + 3]};
      if (dirac_coordinate[mu][i * 4] > data_pattern_conf.lat_dim[0] / 2) {
        decomposition_step_parallel1(dirac_difference[mu][i], coordinate,
                                     laplace, angles_decomposed[mu],
                                     data_pattern_conf, data_pattern_laplacian);
      } else {
        decomposition_step_parallel2(dirac_difference[mu][i], coordinate,
                                     laplace, angles_decomposed[mu],
                                     data_pattern_conf, data_pattern_laplacian);
      }
    }
  }
  for (int mu = 0; mu < 4; mu++) {
    for (int i = 0; i < angles_decomposed[mu].size(); i++) {
      angles_decomposed[mu][i] = -2 * M_PI * angles_decomposed[mu][i];
    }
  }
  return merge_angles(angles_decomposed, data_pattern_conf);
}

void calculate_inverse_laplacian(std::vector<double> &inverse_laplacian_real,
                                 std::vector<double> &inverse_laplacian_imag) {
  int data_size = x_size * y_size * z_size * t_size;
  std::vector<double> momentum(4);
  double inv_moment_sum;
  link1 link(x_size, y_size, z_size, t_size);
  link1 link_laplacian(x_size / 2 + 1, y_size / 2 + 1, z_size / 2 + 1,
                       t_size / 2 + 1);
  double scalar_mult;
  inverse_laplacian_real =
      std::vector<double>((x_size / 2 + 1) * (y_size / 2 + 1) *
                          (z_size / 2 + 1) * (t_size / 2 + 1));
  inverse_laplacian_imag =
      std::vector<double>((x_size / 2 + 1) * (y_size / 2 + 1) *
                          (z_size / 2 + 1) * (t_size / 2 + 1));

  for (int t1 = 0; t1 < t_size / 2 + 1; t1++) {
    for (int z1 = 0; z1 < z_size / 2 + 1; z1++) {
      for (int y1 = 0; y1 < y_size / 2 + 1; y1++) {
        for (int x1 = 0; x1 < x_size / 2 + 1; x1++) {
          link_laplacian.go_update(x1, y1, z1, t1);

          for (int t = 0; t < t_size; t++) {
            for (int z = 0; z < z_size; z++) {
              for (int y = 0; y < y_size; y++) {
                for (int x = 0; x < x_size; x++) {
                  if (x != 0 && y != 0 && z != 0 && t != 0) {
                    link.go_update(x, y, z, t);
                    for (int mu = 0; mu < 4; mu++) {
                      momentum[mu] = 2 * M_PI * link.coordinate[mu] /
                                     link.lattice_size[mu];
                    }
                    inv_moment_sum = 0;
                    for (int mu = 0; mu < 4; mu++) {
                      inv_moment_sum += 2 - 2 * cos(momentum[mu]);
                    }
                    inv_moment_sum = 1 / inv_moment_sum;
                    scalar_mult = 0;
                    for (int mu = 0; mu < 4; mu++) {
                      scalar_mult += link.coordinate[mu] * momentum[mu];
                    }
                    inverse_laplacian_real[link_laplacian.place] +=
                        cos(scalar_mult) * inv_moment_sum;
                    inverse_laplacian_imag[link_laplacian.place] +=
                        sin(scalar_mult) * inv_moment_sum;
                  }
                }
              }
            }
          }
          inverse_laplacian_real[link_laplacian.place] =
              inverse_laplacian_real[link_laplacian.place] / data_size;
          inverse_laplacian_imag[link_laplacian.place] =
              inverse_laplacian_imag[link_laplacian.place] / data_size;
        }
      }
    }
  }
}

std::vector<double> calculate_inverse_laplacian(
    DataPatternLexicographical &data_pattern_conf,
    DataPatternLexicographical &data_pattern_laplacian) {
  int data_size = data_pattern_conf.get_lattice_size();
  std::vector<double> momentum(4);
  std::vector<double> inverse_laplacian(
      data_pattern_laplacian.get_lattice_size());
  double inv_moment_sum;
  double scalar_mult;
  double sin_p;
  int index_laplacian;
  for (int t1 = 0; t1 < data_pattern_laplacian.lat_dim[3]; t1++) {
    for (int z1 = 0; z1 < data_pattern_laplacian.lat_dim[2]; z1++) {
      for (int y1 = 0; y1 < data_pattern_laplacian.lat_dim[1]; y1++) {
        for (int x1 = 0; x1 < data_pattern_laplacian.lat_dim[0]; x1++) {
          data_pattern_laplacian.lat_coord = {x1, y1, z1, t1};
          index_laplacian = data_pattern_laplacian.get_index_site();
          for (int t = 0; t < data_pattern_conf.lat_dim[3]; t++) {
            for (int z = 0; z < data_pattern_conf.lat_dim[2]; z++) {
              for (int y = 0; y < data_pattern_conf.lat_dim[1]; y++) {
                for (int x = 0; x < data_pattern_conf.lat_dim[0]; x++) {
                  if (!(x == 0 && y == 0 && z == 0 && t == 0)) {
                    data_pattern_conf.lat_coord = {x, y, z, t};
                    for (int mu = 0; mu < 4; mu++) {
                      momentum[mu] = 2 * M_PI *
                                     data_pattern_conf.lat_coord[mu] /
                                     data_pattern_conf.lat_dim[mu];
                    }
                    inv_moment_sum = 0;
                    for (int mu = 0; mu < 4; mu++) {
                      sin_p = sin(momentum[mu] / 2);
                      inv_moment_sum += 4 * sin_p * sin_p;
                    }
                    inv_moment_sum = 1 / inv_moment_sum;
                    scalar_mult = 0;
                    for (int mu = 0; mu < 4; mu++) {
                      scalar_mult +=
                          data_pattern_laplacian.lat_coord[mu] * momentum[mu];
                    }
                    inverse_laplacian[index_laplacian] +=
                        cos(scalar_mult) * inv_moment_sum;
                  }
                }
              }
            }
          }
          inverse_laplacian[index_laplacian] /= data_size;
        }
      }
    }
  }
  return inverse_laplacian;
}

std::vector<double>
calculate_inverse_laplacian(DataPatternLexicographical &data_pattern_conf,
                            DataPatternLexicographical &data_pattern_laplacian,
                            const std::vector<double> &inverse_momenta) {
  int data_size = data_pattern_conf.get_lattice_size();
  std::vector<double> inverse_laplacian(
      data_pattern_laplacian.get_lattice_size());
  double scalar_mult;
  double sin_p;
  int index_laplacian;
  double inverse_laplacian_tmp;
  for (int t1 = 0; t1 < data_pattern_laplacian.lat_dim[3]; t1++) {
    for (int z1 = 0; z1 < data_pattern_laplacian.lat_dim[2]; z1++) {
      for (int y1 = 0; y1 < data_pattern_laplacian.lat_dim[1]; y1++) {
        for (int x1 = 0; x1 < data_pattern_laplacian.lat_dim[0]; x1++) {
          data_pattern_laplacian.lat_coord = {x1, y1, z1, t1};
          index_laplacian = data_pattern_laplacian.get_index_site();
          inverse_laplacian_tmp = 0;
#pragma omp parallel for collapse(4) firstprivate(                             \
        data_pattern_laplacian, data_pattern_conf) private(scalar_mult, sin_p) \
    reduction(+ : inverse_laplacian_tmp)
          for (int t = 0; t < data_pattern_conf.lat_dim[3]; t++) {
            for (int z = 0; z < data_pattern_conf.lat_dim[2]; z++) {
              for (int y = 0; y < data_pattern_conf.lat_dim[1]; y++) {
                for (int x = 0; x < data_pattern_conf.lat_dim[0]; x++) {
                  if (!(x == 0 && y == 0 && z == 0 && t == 0)) {
                    data_pattern_conf.lat_coord = {x, y, z, t};
                    scalar_mult = 0;
                    for (int mu = 0; mu < 4; mu++) {
                      scalar_mult += data_pattern_laplacian.lat_coord[mu] * 2 *
                                     M_PI * data_pattern_conf.lat_coord[mu] /
                                     data_pattern_conf.lat_dim[mu];
                    }
                    inverse_laplacian_tmp +=
                        cos(scalar_mult) *
                        inverse_momenta[data_pattern_conf.get_index_site()];
                  }
                }
              }
            }
          }
          inverse_laplacian[index_laplacian] =
              inverse_laplacian_tmp / data_size;
        }
      }
    }
  }
  return inverse_laplacian;
}

std::vector<double> calculate_inverse_laplacian_spatially_symmetric(
    DataPatternLexicographical &data_pattern_conf,
    DataPatternLexicographical &data_pattern_laplacian,
    const std::vector<double> &inverse_momenta) {
  int data_size = data_pattern_conf.get_lattice_size();
  std::vector<double> inverse_laplacian(
      data_pattern_laplacian.get_lattice_size());
  double scalar_mult;
  double sin_p;
  double inverse_laplacian_tmp;
  std::array<int, 3> momentum_spatial;
  for (int t1 = 0; t1 < data_pattern_laplacian.lat_dim[3]; t1++) {
    for (int z1 = 0; z1 < data_pattern_laplacian.lat_dim[2]; z1++) {
      for (int y1 = z1; y1 < data_pattern_laplacian.lat_dim[1]; y1++) {
        for (int x1 = y1; x1 < data_pattern_laplacian.lat_dim[0]; x1++) {
          data_pattern_laplacian.lat_coord = {x1, y1, z1, t1};
          inverse_laplacian_tmp = 0;
#pragma omp parallel for collapse(4) firstprivate(                             \
        data_pattern_laplacian, data_pattern_conf) private(scalar_mult, sin_p) \
    reduction(+ : inverse_laplacian_tmp)
          for (int t = 0; t < data_pattern_conf.lat_dim[3]; t++) {
            for (int z = 0; z < data_pattern_conf.lat_dim[2]; z++) {
              for (int y = 0; y < data_pattern_conf.lat_dim[1]; y++) {
                for (int x = 0; x < data_pattern_conf.lat_dim[0]; x++) {
                  if (!(x == 0 && y == 0 && z == 0 && t == 0)) {
                    data_pattern_conf.lat_coord = {x, y, z, t};
                    scalar_mult = 0;
                    for (int mu = 0; mu < 4; mu++) {
                      scalar_mult += data_pattern_laplacian.lat_coord[mu] * 2 *
                                     M_PI * data_pattern_conf.lat_coord[mu] /
                                     data_pattern_conf.lat_dim[mu];
                    }
                    inverse_laplacian_tmp +=
                        cos(scalar_mult) *
                        inverse_momenta[data_pattern_conf.get_index_site()];
                  }
                }
              }
            }
          }
          inverse_laplacian[data_pattern_laplacian.get_index_site()] =
              inverse_laplacian_tmp / data_size;
          momentum_spatial = {z1, y1, x1};
          while (std::next_permutation(momentum_spatial.begin(),
                                       momentum_spatial.end())) {
            for (int i = 0; i < 3; i++) {
              data_pattern_laplacian.lat_coord[i] = momentum_spatial[2 - i];
            }
            inverse_laplacian[data_pattern_laplacian.get_index_site()] =
                inverse_laplacian_tmp / data_size;
          }
        }
      }
    }
  }
  return inverse_laplacian;
}

std::vector<double>
make_inverse_momenta(DataPatternLexicographical &data_pattern_conf) {
  std::vector<double> inverse_momenta(data_pattern_conf.get_lattice_size());
  double inv_moment_sum;
  double sin_p;
  std::array<double, 4> momentum;
  for (int t = 0; t < data_pattern_conf.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern_conf.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern_conf.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern_conf.lat_dim[0]; x++) {
          if (!(x == 0 && y == 0 && z == 0 && t == 0)) {
            data_pattern_conf.lat_coord = {x, y, z, t};
            inv_moment_sum = 0;
            for (int mu = 0; mu < 4; mu++) {
              sin_p = sin(M_PI * data_pattern_conf.lat_coord[mu] /
                          data_pattern_conf.lat_dim[mu]);
              inv_moment_sum += 4 * sin_p * sin_p;
            }
            inverse_momenta[data_pattern_conf.get_index_site()] =
                1 / inv_moment_sum;
          }
        }
      }
    }
  }
  return inverse_momenta;
}