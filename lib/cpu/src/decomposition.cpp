#include "../include/decomposition.h"
#include "../include/monopoles.h"

#include <math.h>

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

std::vector<double> read_inverse_laplacian(std::string &file_path) {
  int data_size =
      (x_size / 2 + 1) * (y_size / 2 + 1) * (z_size / 2 + 1) * (t_size / 2 + 1);
  std::vector<double> laplace(data_size);
  std::vector<double> v(data_size);
  std::ifstream stream(file_path);
  stream.ignore(4);
  if (!stream.read((char *)&v[0], data_size * sizeof(double)))
    std::cout << "read_float_convert_abelian<abelian> error: " << file_path
              << std::endl;

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

std::vector<double> make_monopole_angles(std::vector<double> &angles,
                                         std::vector<double> &laplace) {

  link1 link(x_size, y_size, z_size, t_size);
  std::vector<std::vector<int>> monopole_plaket =
      calculate_monopole_plaket_singular(angles);

  std::vector<double> monopole_angles(4 * x_size * y_size * z_size * t_size);

  SPACE_ITER_START

  for (int mu = 0; mu < 4; mu++) {

    monopole_angles[link.place + mu] =
        get_monopole_angle(monopole_plaket, link, laplace, mu);
    // std::cout << monopole_angles[link.place + mu] << std::endl;
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
  for (int t = 0; t < 1; t++) {
    for (int z = 0; z < 1; z++) {
      for (int y = 0; y < 1; y++) {
        for (int x = 0; x < 1; x++) {

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

  for (int mu = 0; mu < 4; mu++) {
    std::cout << "monopole_difference size " << monopole_difference[mu].size()
              << std::endl;
  }

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

      // angles_place = 0;
      // for (int t = 0; t < 1; t++) {
      //   for (int z = 0; z < 1; z++) {
      //     for (int y = 0; y < 1; y++) {
      //       for (int x = 0; x < x_size; x++) {

      //         laplace_coordinate[0] = abs(x - monopole_coordinate[mu][i *
      //         4]); if (laplace_coordinate[0] > x_size / 2)
      //           laplace_coordinate[0] = x_size - laplace_coordinate[0];

      //         laplace_coordinate[1] =
      //             abs(y - monopole_coordinate[mu][i * 4 + 1]);
      //         if (laplace_coordinate[1] > y_size / 2)
      //           laplace_coordinate[1] = y_size - laplace_coordinate[1];

      //         laplace_coordinate[2] =
      //             abs(z - monopole_coordinate[mu][i * 4 + 2]);
      //         if (laplace_coordinate[2] > z_size / 2)
      //           laplace_coordinate[2] = z_size - laplace_coordinate[2];

      //         laplace_coordinate[3] =
      //             abs(t - monopole_coordinate[mu][i * 4 + 3]);
      //         if (laplace_coordinate[3] > t_size / 2)
      //           laplace_coordinate[3] = t_size - laplace_coordinate[3];

      //         laplace_place =
      //             (laplace_coordinate[3]) * laplace_size[0] * laplace_size[1]
      //             *
      //                 laplace_size[2] +
      //             (laplace_coordinate[2]) * laplace_size[0] * laplace_size[1]
      //             + (laplace_coordinate[1]) * laplace_size[0] +
      //             laplace_coordinate[0];

      //         angles_decomposed[angles_place * 4 + mu] +=
      //             laplace[laplace_place] * monopole_difference[mu][i];

      //         angles_place++;
      //       }
      //     }
      //   }
      // }
    }
  }
  return angles_decomposed;
}