#include "../include/monopoles.h"
#include <cmath>
#include <link.h>

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
      (t)*x_size *y_size *z_size * 4 + (z)*x_size *y_size * 4 +                \
      (y)*x_size * 4 + (x)*4

// read configuration of angles
std::vector<FLOAT> read_angles_float_fortran(std::string &file_path) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  std::vector<FLOAT> angles;
  angles.reserve(data_size1);
  std::ifstream stream(file_path);
  std::vector<float> v(data_size1 + 2);
  if (!stream.read((char *)&v[0], (data_size1 + 2) * sizeof(float)))
    std::cout << "read_angles_float_fortran error: " << file_path << std::endl;
  for (int i = 0; i < data_size1; i++) {
    angles.push_back((FLOAT)v[i + 1]);
  }
  stream.close();
  return angles;
}

std::vector<FLOAT> read_angles_double_fortran(std::string &file_path) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  std::vector<FLOAT> angles;
  angles.reserve(data_size1);
  std::ifstream stream(file_path);
  std::vector<double> v(data_size1);
  stream.ignore(4);
  if (!stream.read((char *)&v[0], (data_size1) * sizeof(double)))
    std::cout << "read_angles_double_fortran error: " << file_path << std::endl;
  for (int i = 0; i < data_size1; i++) {
    angles.push_back((FLOAT)v[i]);
  }
  stream.close();
  return angles;
}

// read configuration of su2 matrices and extract abelian degrees of freedom
std::vector<FLOAT> read_float_fortran_convet_abelian(std::string &file_path) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  std::vector<FLOAT> angles;
  angles.reserve(data_size1 * 4);
  std::ifstream stream(file_path);
  std::vector<float> v(data_size1 * 4 + 1);
  if (!stream.read((char *)&v[0], (data_size1 * 4 + 1) * sizeof(float)))
    std::cout << "read_float_convert_abelian<abelian> error: " << file_path
              << std::endl;
  for (int i = 0; i < data_size1; i++) {
    angles.push_back((FLOAT)atan2(v[i * 4 + 4], v[i * 4 + 1]));
  }
  stream.close();
  return angles;
}

std::vector<FLOAT> read_double_fortran_convet_abelian(std::string &file_path) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  std::vector<FLOAT> angles;
  angles.reserve(data_size1 * 4);
  std::ifstream stream(file_path);
  std::vector<double> v(data_size1 * 4);
  stream.ignore(4);
  if (!stream.read((char *)&v[0], (data_size1 * 4) * sizeof(double)))
    std::cout << "read_float_convert_abelian<abelian> error: " << file_path
              << std::endl;
  for (int i = 0; i < data_size1; i++) {
    angles.push_back((FLOAT)atan2(v[i * 4 + 3], v[i * 4]));
  }
  stream.close();
  return angles;
}

std::vector<FLOAT> read_double_qc2dstag_convet_abelian(std::string &file_path) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  std::vector<FLOAT> angles;
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
    angles.push_back((FLOAT)atan2(v[PLACE_QC2DSTAG + 3], v[PLACE_QC2DSTAG]));
  }
  SPACE_ITER_END_SIMPLE
  stream.close();
  return angles;
}

std::vector<FLOAT> read_inverse_laplacian(std::string &file_path) {
  int data_size =
      (x_size / 2 + 1) * (y_size / 2 + 1) * (z_size / 2 + 1) * (t_size / 2 + 1);
  std::vector<FLOAT> laplace(data_size);
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

// calculate monopole_plaketes on the lattice
std::vector<std::vector<FLOAT>>
calculate_monopole_plaket(std::vector<FLOAT> &angles) {
  int data_size = x_size * y_size * z_size * t_size;
  std::vector<std::vector<FLOAT>> plakets(6, std::vector<FLOAT>(data_size));
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

double get_monopole_angle(std::vector<std::vector<FLOAT>> &monopole_plaket,
                          link1 &link_tmp, std::vector<FLOAT> &laplace,
                          int mu) {
  link1 link(x_size, y_size, z_size, t_size);

  link1 link_laplace(x_size / 2 + 1, y_size / 2 + 1, z_size / 2 + 1,
                     t_size / 2 + 1);

  double monopole_angle = 0;
  double angle_tmp;

  SPACE_ITER_START

  link_laplace.go_update(0, 0, 0, 0);
  link_laplace.move(0, link_tmp.coordinate[0] - x);
  link_laplace.move(0, link_tmp.coordinate[1] - y);
  link_laplace.move(0, link_tmp.coordinate[2] - z);
  link_laplace.move(0, link_tmp.coordinate[3] - t);

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
    }
  }

  SPACE_ITER_END

  return monopole_angle;
}

std::vector<FLOAT> make_monopole_angles(std::vector<FLOAT> &angles,
                                        std::vector<FLOAT> &laplace) {

  link1 link(x_size, y_size, z_size, t_size);
  std::vector<std::vector<FLOAT>> monopole_plaket =
      calculate_monopole_plaket(angles);

  std::vector<FLOAT> monopole_angles(4 * x_size * y_size * z_size * t_size);

  SPACE_ITER_START

  for (int mu = 0; mu < 4; mu++) {

    monopole_angles[link.place + mu] =
        get_monopole_angle(monopole_plaket, link, laplace, mu);
  }

  SPACE_ITER_END

  return monopole_angles;
}

std::vector<FLOAT> calculate_current(std::vector<FLOAT> &angles) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  link1 link(x_size, y_size, z_size, t_size);

  std::vector<FLOAT> J(data_size);

  std::vector<std::vector<FLOAT>> monopole_plaket =
      calculate_monopole_plaket(angles);

  SPACE_ITER_START

  link.get_current(monopole_plaket, &J[link.place], angles);

  SPACE_ITER_END

  return J;
}

int find_current(link1 &link, std::vector<FLOAT> &J) {
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

// find all directions with current and add loops with that neighbours
std::vector<loop *> find_paths(std::vector<loop *> &neighbours,
                               std::vector<FLOAT> &J) {
  // std::vector for new sites with current
  std::vector<loop *> neighbours_new;
  loop *loop_tmp;
  link1 link(x_size, y_size, z_size, t_size);
  int J_tmp;

  // for all previously found sites of cluster find new sites with current
  for (int i = 0; i < neighbours.size(); i++) {
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
        neighbours.push_back(loop_tmp);
        link.move(mu, -1);
      }
    }
    for (int mu = 0; mu < 4; mu++) {
      link.move(mu, -1);
      if (J[link.place + mu] < -0.3) {
        J_tmp = std::lround(J[link.place + mu]);
        J[link.place + mu] = 0.;
        loop_tmp = new loop(link);
        // loop_tmp->link.push_back(neighbours[i]);
        neighbours[i]->link.push_back(loop_tmp);
        neighbours[i]->charge.push_back(J_tmp);
        neighbours.push_back(loop_tmp);
      }
      link.move(mu, 1);
    }
  }

  return neighbours_new;
}

// find cluster which has site ll
void find_cluster(loop *ll, std::vector<FLOAT> &J) {
  std::vector<loop *> neighbours = {ll};

  // while find_path finds new sites with current
  do {
    // find new sites with current and add them to loops
    neighbours = find_paths(neighbours, J);
  } while (neighbours.size() > 0);
}

// find all clusters on a lattice not using recurrence
std::vector<loop *> calculate_clusters(std::vector<FLOAT> &J) {
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

// functions for obtaining information about clusters for testing

void print_currents(loop *ll) {
  std::cout << ll->coordinate[0] << " " << ll->coordinate[1] << " "
            << ll->coordinate[2] << " " << ll->coordinate[3] << std::endl;
  for (int i = 0; i < ll->link.size(); i++) {
    print_currents(ll->link[i]);
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
  link1 link(x_size, y_size, z_size, t_size);

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