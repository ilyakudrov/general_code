#include "../include/monopoles.h"

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
    angles.push_back(
        (FLOAT)atan2(v[PLACE_QC2DSTAG + 3], v[PLACE_QC2DSTAG + 1]));
  }
  SPACE_ITER_END_SIMPLE
  stream.close();
  return angles;
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

// int loop::get_dir(int i) {
//   loop *tmp = link[i - 1];
//   if (tmp == NULL)
//     return 0;

//   int delta = (tmp->node.coordinate[0]) - node.coordinate[0];
//   if (delta == (1 - x_size))
//     return 1;
//   if (delta == (x_size - 1))
//     return -1;
//   if (delta != 0)
//     return 1 * delta;

//   delta = (tmp->node.coordinate[1]) - node.coordinate[1];
//   if (delta == (1 - y_size))
//     return 2;
//   if (delta == (y_size - 1))
//     return -2;
//   if (delta != 0)
//     return delta * 2;

//   delta = (tmp->node.coordinate[2]) - node.coordinate[2];
//   if (delta == (1 - z_size))
//     return 3;
//   if (delta == (z_size - 1))
//     return -3;
//   if (delta != 0)
//     return delta * 3;

//   delta = (tmp->node.coordinate[3]) - node.coordinate[3];
//   if (delta == (1 - t_size))
//     return 4;
//   if (delta == (t_size - 1))
//     return -4;
//   if (delta != 0)
//     return delta * 4;
//   return 0;
// }

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

  // for all previously found sites of cluster find new sites with current
  for (int i = 0; i < neighbours.size(); i++) {
    // go to site
    link.go_update(neighbours[i]->coordinate[0], neighbours[i]->coordinate[1],
                   neighbours[i]->coordinate[2], neighbours[i]->coordinate[3]);
    // check all directions for a current
    // if a current is found, add it to cluster and to std::vector of new sites
    for (int mu = 0; mu < 4; mu++) {
      if (J[link.place + mu] > 0.3) {
        J[link.place + mu] = 0.;
        link.move(mu, 1);
        loop_tmp = new loop(link);
        neighbours[i]->link.push_back(loop_tmp);
        neighbours_new.push_back(loop_tmp);
        link.move(mu, -1);
      } else if (J[link.place + mu] < -0.3) {
        J[link.place + mu] = 0.;
        link.move(mu, 1);
        loop_tmp = new loop(link);
        loop_tmp->link.push_back(neighbours[i]);
        neighbours_new.push_back(loop_tmp);
        link.move(mu, -1);
      }
    }
    for (int mu = 0; mu < 4; mu++) {
      link.move(mu, -1);
      if (J[link.place + mu] > 0.3) {
        J[link.place + mu] = 0.;
        loop_tmp = new loop(link);
        loop_tmp->link.push_back(neighbours[i]);
        neighbours_new.push_back(loop_tmp);
      }
      link.move(mu, 1);
      if (J[link.place + mu] < -0.3) {
        J[link.place + mu] = 0.;
        loop_tmp = new loop(link);
        neighbours[i]->link.push_back(loop_tmp);
        neighbours_new.push_back(loop_tmp);
      }
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

// monopole observables

void cluster_length(loop *ll, int &length) {
  for (int i = 0; i < ll->link.size(); i++) {
    length++;
    cluster_length(ll->link[i], length);
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

void length_mu(loop *ll, std::vector<int> &lengths_mu) {
  int dir = 0;
  for (int i = 0; i < ll->link.size(); i++) {
    length_mu(ll->link[i], lengths_mu);
    for (int mu = 0; mu < 4; mu++) {
      lengths_mu[mu] += ll->link[i]->coordinate[mu] - ll->coordinate[mu];
    }
  }
}