#include "../include/decomposition.h"
#include "../include/monopoles.h"

#include <ctime>
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
    std::cout << "read_double<abelian> error: " << file_name << std::endl;
  return angles;
}

void write_double_angles(std::string &file_name, std::vector<double> &angles) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  std::ofstream stream(file_name);
  if (!stream.write((char *)&angles[0], (data_size) * sizeof(double)))
    std::cout << "write_double<abelian> error: " << file_name << std::endl;
}

void write_double_su2(std::string &file_name, std::vector<su2> &conf_su2) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  std::ofstream stream(file_name);
  if (!stream.write((char *)&conf_su2[0], data_size * 4 * sizeof(double)))
    std::cout << "write_double_su2 error: " << file_name << std::endl;
  stream.close();
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

std::vector<su2> get_monopoless(std::vector<su2> &conf_su2,
                                std::vector<double> &angles_monopole) {
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
                              std::vector<double> &angles_monopole) {
  su2 A;

  for (int i = 0; i < conf_su2.size(); i++) {
    A = su2(cos(angles_monopole[i]), 0, 0, sin(angles_monopole[i]));

    conf_su2[i] = conf_su2[i] * A;
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

std::vector<double> read_inverse_laplacian(std::string &file_path) {
  int data_size =
      (x_size / 2 + 1) * (y_size / 2 + 1) * (z_size / 2 + 1) * (t_size / 2 + 1);
  std::vector<double> laplace(data_size);
  std::vector<double> v(data_size);
  std::ifstream stream(file_path);
  stream.ignore(4);
  if (!stream.read((char *)&v[0], data_size * sizeof(double)))
    std::cout << "read_inverse_laplacian error: " << file_path << std::endl;

  link1 link_laplace(x_size / 2 + 1, y_size / 2 + 1, z_size / 2 + 1,
                     t_size / 2 + 1);

  int index = 0;

  for (int x = 0; x < x_size / 2 + 1; x++) {
    for (int y = 0; y < y_size / 2 + 1; y++) {
      for (int z = 0; z < z_size / 2 + 1; z++) {
        for (int t = 0; t < t_size / 2 + 1; t++) {
          link_laplace.go_update(x, y, z, t);
          laplace[link_laplace.place / 4] = v[index];
          index++;
        }
      }
    }
  }

  return laplace;
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

void monopole_plaket_difference_nonzero(
    std::vector<std::vector<int>> &monopole_plaket,
    std::vector<std::vector<int>> &monopole_difference,
    std::vector<std::vector<int>> &monopole_coordinate) {

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

        angle_tmp = monopole_plaket[index][link.place / 4];

        link.move(nu, -1);

        diff += factor * (angle_tmp - monopole_plaket[index][link.place / 4]);

        link.move(nu, 1);
      }
    }

    if (diff != 0) {
      monopole_difference[mu].push_back(diff);
      for (int nu = 0; nu < 4; nu++) {
        monopole_coordinate[mu].push_back(link.coordinate[nu]);
      }
    }
  }

  SPACE_ITER_END
}

void decomposition_step(std::vector<std::vector<int>> &monopole_difference,
                        std::vector<std::vector<int>> &monopole_coordinate,
                        std::vector<double> &laplace,
                        std::vector<double> &angles_decomposed, int i, int mu) {
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

std::vector<double> make_monopole_angles1(std::vector<double> &angles,
                                          std::vector<double> &laplace) {
  link1 link(x_size, y_size, z_size, t_size);

  int data_size = 4 * x_size * y_size * z_size * t_size;

  std::vector<double> angles_decomposed(data_size);

  std::vector<std::vector<int>> monopole_plaket =
      calculate_monopole_plaket_singular(angles);

  std::vector<std::vector<int>> monopole_difference(4, std::vector<int>());
  std::vector<std::vector<int>> monopole_coordinate(4, std::vector<int>());

  monopole_plaket_difference_nonzero(monopole_plaket, monopole_difference,
                                     monopole_coordinate);

  std::vector<int> laplace_size = {x_size / 2 + 1, y_size / 2 + 1,
                                   z_size / 2 + 1, t_size / 2 + 1};

  double monopole_angle = 0;
  int angle_tmp;

  std::vector<int> laplace_coordinate(4);
  int laplace_place;
  int angles_place;

  for (int mu = 0; mu < 4; mu++) {
    for (int i = 0; i < monopole_difference[mu].size(); i++) {

      decomposition_step(monopole_difference, monopole_coordinate, laplace,
                         angles_decomposed, i, mu);
    }
  }

  for (int i = 0; i < angles_decomposed.size(); i++) {
    angles_decomposed[i] = -2 * M_PI * angles_decomposed[i];
  }
  return angles_decomposed;
}

void decomposition_step3(std::vector<std::vector<int>> &monopole_difference,
                         std::vector<std::vector<int>> &monopole_coordinate,
                         std::vector<double> &laplace,
                         std::vector<std::vector<double>> &angles_decomposed,
                         int i, int mu) {
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
make_monopole_angles3(std::vector<double> &angles,
                      std::vector<double> &laplace) {
  link1 link(x_size, y_size, z_size, t_size);

  int data_size = x_size * y_size * z_size * t_size;

  std::vector<std::vector<double>> angles_decomposed(
      4, std::vector<double>(data_size));

  std::vector<std::vector<int>> monopole_plaket =
      calculate_monopole_plaket_singular(angles);

  std::vector<std::vector<int>> monopole_difference(4, std::vector<int>());
  std::vector<std::vector<int>> monopole_coordinate(4, std::vector<int>());

  monopole_plaket_difference_nonzero(monopole_plaket, monopole_difference,
                                     monopole_coordinate);

  std::vector<int> laplace_size = {x_size / 2 + 1, y_size / 2 + 1,
                                   z_size / 2 + 1, t_size / 2 + 1};

  double monopole_angle = 0;
  int angle_tmp;

  std::vector<int> laplace_coordinate(4);
  int laplace_place;
  int angles_place;

  for (int mu = 0; mu < 4; mu++) {
    for (int i = 0; i < monopole_difference[mu].size(); i++) {

      decomposition_step3(monopole_difference, monopole_coordinate, laplace,
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
                        std::vector<int> &monopole_coordinate,
                        std::vector<double> &laplace,
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

std::vector<double> make_monopole_angles(std::vector<double> &angles,
                                         std::vector<double> &laplace) {

  link1 link(x_size, y_size, z_size, t_size);

  int data_size = x_size * y_size * z_size * t_size;

  std::vector<std::vector<int>> monopole_plaket =
      calculate_monopole_plaket_singular(angles);

  std::vector<std::vector<int>> monopole_difference(4, std::vector<int>());
  std::vector<std::vector<int>> monopole_coordinate(4, std::vector<int>());

  monopole_plaket_difference_nonzero(monopole_plaket, monopole_difference,
                                     monopole_coordinate);

  for (int mu = 0; mu < monopole_plaket.size(); mu++) {
    monopole_plaket[mu].clear();
    monopole_plaket[mu].shrink_to_fit();
  }

  std::vector<std::vector<double>> angles_decomposed(
      4, std::vector<double>(data_size));

  std::vector<int> coordinate(4);

  for (int mu = 0; mu < 4; mu++) {

    for (int i = 0; i < monopole_difference[mu].size(); i++) {

      coordinate = {monopole_coordinate[mu][4 * i],
                    monopole_coordinate[mu][4 * i + 1],
                    monopole_coordinate[mu][4 * i + 2],
                    monopole_coordinate[mu][4 * i + 3]};

      decomposition_step(monopole_difference[mu][i], coordinate, laplace,
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
                                 std::vector<int> &monopole_coordinate,
                                 std::vector<double> &laplace,
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
    angles_place, laplace_place, laplace_coordinate, laplace_shift,            \
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
    std::vector<int> &monopole_coordinate, std::vector<double> &laplace,
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
    std::vector<int> &monopole_coordinate, std::vector<double> &laplace,
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
                                  std::vector<int> &monopole_coordinate,
                                  std::vector<double> &laplace,
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
    angles_place, laplace_place, laplace_coordinate, laplace_shift,            \
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
                                  std::vector<int> &monopole_coordinate,
                                  std::vector<double> &laplace,
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
    angles_place, laplace_place, laplace_coordinate, laplace_shift,            \
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
make_monopole_angles_parallel(std::vector<double> &angles,
                              std::vector<double> &laplace) {

  link1 link(x_size, y_size, z_size, t_size);

  int data_size = x_size * y_size * z_size * t_size;

  std::vector<std::vector<int>> monopole_plaket =
      calculate_monopole_plaket_singular(angles);

  std::vector<std::vector<int>> monopole_difference(4, std::vector<int>());
  std::vector<std::vector<int>> monopole_coordinate(4, std::vector<int>());

  monopole_plaket_difference_nonzero(monopole_plaket, monopole_difference,
                                     monopole_coordinate);

  for (int mu = 0; mu < monopole_plaket.size(); mu++) {
    monopole_plaket[mu].clear();
    monopole_plaket[mu].shrink_to_fit();
  }

  std::vector<std::vector<double>> angles_decomposed(
      4, std::vector<double>(data_size));

  std::vector<int> coordinate(4);

  for (int mu = 0; mu < 4; mu++) {

    for (int i = 0; i < monopole_difference[mu].size(); i++) {

      coordinate = {monopole_coordinate[mu][4 * i],
                    monopole_coordinate[mu][4 * i + 1],
                    monopole_coordinate[mu][4 * i + 2],
                    monopole_coordinate[mu][4 * i + 3]};

      if (monopole_coordinate[mu][i * 4] > x_size / 2) {

        decomposition_step_parallel1(monopole_difference[mu][i], coordinate,
                                     laplace, angles_decomposed[mu]);
      } else {
        decomposition_step_parallel2(monopole_difference[mu][i], coordinate,
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