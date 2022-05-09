#include "../include/link.h"
#include "../include/data.h"
#include <cmath>

#define Pi 3.141592653589793238462643383279502884
#define data_size                                                              \
  4 * lattice_size[0] * lattice_size[1] * lattice_size[2] * lattice_size[3]
#define PLACE3_LINK_DIR                                                        \
  (coordinate[3]) * 3 * lattice_size[0] * lattice_size[1] * lattice_size[2] +  \
      (coordinate[2]) * 3 * lattice_size[0] * lattice_size[1] +                \
      (coordinate[1]) * 3 * lattice_size[0] + (coordinate[0]) * 3 + direction

#define PLACE1_LINK_NODIR                                                      \
  (coordinate[3]) * lattice_size[0] * lattice_size[1] * lattice_size[2] +      \
      (coordinate[2]) * lattice_size[0] * lattice_size[1] +                    \
      (coordinate[1]) * lattice_size[0] + (coordinate[0])

link1::link1(int lattice_size_x, int lattice_size_y, int lattice_size_z,
             int lattice_size_t) {
  lattice_size[0] = lattice_size_x;
  lattice_size[1] = lattice_size_y;
  lattice_size[2] = lattice_size_z;
  lattice_size[3] = lattice_size_t;
  multiplier[0] = 4;
  multiplier[1] = lattice_size[0] * 4;
  multiplier[2] = lattice_size[0] * lattice_size[1] * 4;
  multiplier[3] = lattice_size[0] * lattice_size[1] * lattice_size[2] * 4;
  direction = 0;
  for (int i = 0; i < 4; i++) {
    coordinate[i] = 0;
    coordinate_old[i] = 0;
  }
  place = 0;
}

link1::link1(const link1 &link) {
  lattice_size[0] = link.lattice_size[0];
  lattice_size[1] = link.lattice_size[1];
  lattice_size[2] = link.lattice_size[2];
  lattice_size[3] = link.lattice_size[3];
  direction = link.direction;
  multiplier[0] = 4;
  multiplier[1] = lattice_size[0] * 4;
  multiplier[2] = lattice_size[0] * lattice_size[1] * 4;
  multiplier[3] = lattice_size[0] * lattice_size[1] * lattice_size[2] * 4;
  for (int i = 0; i < 4; i++) {
    coordinate[i] = link.coordinate[i];
    coordinate_old[i] = link.coordinate_old[i];
  }
  place = link.place;
}

link1::link1() {
  lattice_size[0] = x_size;
  lattice_size[1] = y_size;
  lattice_size[2] = z_size;
  lattice_size[3] = t_size;
  direction = 0;
  multiplier[0] = 4;
  multiplier[1] = lattice_size[0] * 4;
  multiplier[2] = lattice_size[0] * lattice_size[1] * 4;
  multiplier[3] = lattice_size[0] * lattice_size[1] * lattice_size[2] * 4;
  for (int i = 0; i < 4; i++) {
    coordinate[i] = 0;
    coordinate_old[i] = 0;
  }
  place = 0;
}

std::ostream &operator<<(std::ostream &os, const link1 &link) {
  os << "x: " << link.coordinate[0] << " y: " << link.coordinate[1]
     << " z: " << link.coordinate[2] << " t: " << link.coordinate[3]
     << " dir: " << link.direction << std::endl;
  return os;
}

void link1::move(int dir, int step) {
  coordinate[dir] += (lattice_size[dir] + step);
  coordinate[dir] = coordinate[dir] % lattice_size[dir];
  place += (coordinate[dir] - coordinate_old[dir]) * multiplier[dir];
  coordinate_old[dir] = coordinate[dir];
}

void link1::update(int dir) {
  place += (coordinate[dir] - coordinate_old[dir]) * multiplier[dir];
  coordinate_old[dir] = coordinate[dir];
}

void link1::go(int x, int y, int z, int t) {
  coordinate[0] = x;
  coordinate[1] = y;
  coordinate[2] = z;
  coordinate[3] = t;
}
void link1::move_dir(int dir) { direction = dir; }

void link1::go_update(int x, int y, int z, int t) {
  coordinate[0] = x;
  coordinate[1] = y;
  coordinate[2] = z;
  coordinate[3] = t;

  update(0);
  update(1);
  update(2);
  update(3);
}

template <class T> T link1::plaket_mu(const std::vector<T> &array, int mu) {
  int dir = direction;
  T A = array[place + direction];
  move(dir, 1);
  move_dir(mu);
  A = A * array[place + mu];
  move_dir(dir);
  move(dir, -1);
  move(mu, 1);
  A = A ^ array[place + dir];
  move_dir(mu);
  move(mu, -1);
  A = A ^ array[place + mu];
  move_dir(dir);
  return A;
}

template <class T>
T link1::plaket_schwinger_average(const std::vector<T> &array, int mu) {
  int dir = direction;
  T result;
  T A = array[place + direction];
  move(dir, 1);
  move_dir(mu);
  A = A * array[place + direction];
  move_dir(dir);
  move(dir, -1);
  move(mu, 1);
  A = A ^ array[place + direction];
  move_dir(mu);
  move(mu, -1);
  A = A ^ array[place + direction];

  result = A;

  A = array[place + direction];
  move(mu, 1);
  move(dir, -1);
  move_dir(dir);
  A = A ^ array[place + direction];
  move(mu, -1);
  move_dir(mu);
  A = A ^ array[place + direction];
  move_dir(dir);
  A = A * array[place + direction];

  result = result + A;

  A = array[place + direction].conj();
  move(mu, -1);
  move_dir(mu);
  A = A ^ array[place + direction];
  move_dir(dir);
  A = A * array[place + direction];
  move(dir, 1);
  move(mu, 1);
  move_dir(mu);
  A = A * array[place + direction];

  result = result + A;

  A = array[place + direction].conj();
  move_dir(dir);
  A = A * array[place + direction];
  move(dir, 1);
  move_dir(mu);
  A = A * array[place + direction];
  move(mu, 1);
  move(dir, -1);
  move_dir(dir);
  A = A ^ array[place + direction];

  result = result + A;

  return result * (1. / 4);
}

// TODO: elaborate
template <class T>
T link1::schwinger_line(const std::vector<T> &array, int d, int dir, int x) {
  int dir1 = direction;
  T A;
  for (int i = 0; i < d; i++) {
    A = A * array[place + direction];
    move(dir1, 1);
  }
  move_dir(dir);
  for (int i = 0; i < x; i++) {
    A = A * array[place + direction];
    move(dir, 1);
  }
  move_dir(dir1);
  move(dir, -x);
  move(dir1, -d);
  return A;
}

template <class T> T link1::polyakov_loop(const std::vector<T> &array) {
  T A;
  for (int i = 0; i < lattice_size[3]; i++) {
    A = A * array[place + direction];
    move(direction, 1);
  }
  return A;
}

template <class T>
T link1::wilson_loop(const std::vector<T> &array, int r, int t) {
  int dir = direction;
  T A;
  for (int i = 0; i < r; i++) {
    A = A * array[place + direction];
    move(dir, 1);
  }
  move_dir(3);
  for (int i = 0; i < t; i++) {
    A = A * array[place + direction];
    move(3, 1);
  }
  move_dir(dir);
  for (int i = 0; i < r; i++) {
    move(dir, -1);
    A = A ^ array[place + direction];
  }
  move_dir(3);
  for (int i = 0; i < t; i++) {
    move(3, -1);
    A = A ^ array[place + direction];
  }
  move_dir(dir);
  return A;
}

double link1::wilson_loop_abelian(const std::vector<double> &array, int r,
                                  int t) {
  int dir = direction;
  double A = 0;
  for (int i = 0; i < r; i++) {
    A = A + array[place + dir];
    move(dir, 1);
  }
  for (int i = 0; i < t; i++) {
    A = A + array[place + 3];
    move(3, 1);
  }
  for (int i = 0; i < r; i++) {
    move(dir, -1);
    A = A - array[place + dir];
  }
  for (int i = 0; i < t; i++) {
    move(3, -1);
    A = A - array[place + 3];
  }
  move_dir(dir);
  return A;
}

template <class T>
T link1::wilson_loop_schwinger(const std::vector<T> &array, int r, int t) {
  int dir = direction;
  T A;
  move_dir(3);
  for (int i = 0; i < t / 2; i++) {
    A = A * array[place + direction];
    move(3, 1);
  }
  move_dir(dir);
  for (int i = 0; i < r; i++) {
    A = A * array[place + direction];
    move(dir, 1);
  }
  move_dir(3);
  move(3, -1);
  for (int i = 0; i < t - 1; i++) {
    A = A ^ array[place + direction];
    move(3, -1);
  }
  move_dir(dir);
  move(dir, -1);
  for (int i = 0; i < r - 1; i++) {
    A = A ^ array[place + direction];
    move(dir, -1);
  }
  move_dir(3);
  for (int i = 0; i < t / 2; i++) {
    A = A * array[place + direction];
    move(3, 1);
  }
  move_dir(dir);
  return A;
}

template <class T>
T link1::wilson_line(const std::vector<T> &array, int length) {
  int dir = direction;
  T A;
  for (int i = 0; i < length; i++) {
    A = A * array[place + direction];
    move(dir, 1);
  }
  return A;
}

template <class T>
T link1::wilson_line_single(const std::vector<T> &array, int length) {
  T A;
  for (int i = 0; i < length; i++) {
    A = A * array[place / 4];
    move(direction, 1);
  }
  return A;
}

// calculates offaxis spatial line of wilson loop according to pattern
// pattern defines directions of links along the line
// it consists of positive and negative integers which correspond to spatial
// directions x - 1, y - 2, z - 3
template <class T>
T link1::wilson_line_offaxis(const std::vector<T> &array,
                             const std::vector<int> &pattern) {
  T A;

  // iterate through pattern
  for (int i = 0; i < pattern.size(); i++) {
    move_dir(abs(pattern[i]) - 1);

    // if positive direction
    if (pattern[i] > 0) {
      A = A * array[place + direction];
      move(abs(pattern[i]) - 1, 1);
    }

    // if negative direction
    else if (pattern[i] < 0) {
      move(abs(pattern[i]) - 1, -1);
      A = A ^ array[place + direction];
    } else {
      std::cout << "wilson_line_offaxis pattern direction error " << pattern[i]
                << std::endl;
    }
  }
  return A;
}

// TODO: elaborate directions
template <class T>
double link1::field1(const std::vector<std::vector<T>> &schwinger_line,
                     const std::vector<T> &plaket,
                     const std::vector<T> &polyakov_loop, int d, int D, int dir,
                     int x) {
  int dir1 = direction;
  T C = schwinger_line[dir - 1][PLACE3_LINK_DIR];
  move(dir1, d);
  move(dir, x);
  T B = plaket[PLACE3_LINK_DIR];
  move(dir, -x);
  move(dir1, -d);
  move_dir(-3);
  T A = polyakov_loop[PLACE1_LINK_NODIR];
  A = A.conj() * C * B;
  A = A * C.conj();
  move(dir1, D);
  move_dir(3);
  B = polyakov_loop[PLACE1_LINK_NODIR];
  move(dir1, -D);
  move_dir(dir1);
  return A.tr() * B.tr();
}

// TODO: elaborate directions
template <class T>
double link1::field2(const std::vector<T> &plaket,
                     const std::vector<T> &polyakov_loop, int d, int D, int dir,
                     int x) {
  int dir1 = direction;
  move_dir(-4);
  T A = polyakov_loop[PLACE1_LINK_NODIR];
  move(dir1, d);
  move(dir, x);
  move_dir(dir1);
  T B = plaket[PLACE3_LINK_DIR];
  move(dir, -x);
  move(dir1, D - d);
  move_dir(4);
  T C = polyakov_loop[PLACE1_LINK_NODIR];
  move(dir1, -D);
  move_dir(dir1);
  return B.tr() * C.tr() * A.conj().tr();
}

// TODO: elaborate directions
template <class T>
double link1::field3(const std::vector<T> &polyakov_loop, int D, int x) {
  int dir1 = direction;
  move_dir(-4);
  T A = polyakov_loop[PLACE1_LINK_NODIR];
  move(dir1, D);
  move_dir(4);
  T B = polyakov_loop[PLACE1_LINK_NODIR];
  move(dir1, -D);
  move_dir(dir1);
  return A.conj().tr() * B.tr();
}

// monopoles

// calclate monopole plaket from angles in link.direction-mu plane
double link1::monopole_plaket_mu(std::vector<double> &angles, int mu) {
  double angle = angles[place + direction];
  move(direction, 1);
  angle += angles[place + mu];
  move(direction, -1);
  move(mu, 1);
  angle -= angles[place + direction];
  move(mu, -1);
  angle -= angles[place + mu];
  while ((angle > Pi) || (angle < -Pi)) {
    if (angle > Pi)
      angle = angle - 2 * Pi;
    if (angle < -Pi)
      angle = angle + 2 * Pi;
  }
  return angle;
}

void link1::get_current(std::vector<std::vector<double>> &monopole_plaket,
                        double *J, std::vector<double> &angles) {
  double j0, j1, j2, j3;
  link1 linkx(*this);
  linkx.move(0, 1);
  link1 linky(*this);
  linky.move(1, 1);
  link1 linkz(*this);
  linkz.move(2, 1);
  link1 linkt(*this);
  linkt.move(3, 1);

  move(3, 1);
  linkx.move(3, 1);
  linky.move(3, 1);
  linkz.move(3, 1);
  j3 = monopole_plaket[3][linkx.place / 4] - monopole_plaket[3][place / 4] -
       (monopole_plaket[1][linky.place / 4] - monopole_plaket[1][place / 4]) +
       monopole_plaket[0][linkz.place / 4] - monopole_plaket[0][place / 4];
  move(3, -1);
  linkx.move(3, -1);
  linky.move(3, -1);
  linkz.move(3, -1);
  move(0, 1);
  linkt.move(0, 1);
  linky.move(0, 1);
  linkz.move(0, 1);
  j0 = -(monopole_plaket[3][linkt.place / 4] - monopole_plaket[3][place / 4]) -
       (monopole_plaket[5][linky.place / 4] - monopole_plaket[5][place / 4]) +
       monopole_plaket[4][linkz.place / 4] - monopole_plaket[4][place / 4];
  move(0, -1);
  linkt.move(0, -1);
  linky.move(0, -1);
  linkz.move(0, -1);
  move(1, 1);
  linkt.move(1, 1);
  linkx.move(1, 1);
  linkz.move(1, 1);
  j1 = monopole_plaket[1][linkt.place / 4] - monopole_plaket[1][place / 4] +
       monopole_plaket[5][linkx.place / 4] - monopole_plaket[5][place / 4] -
       (monopole_plaket[2][linkz.place / 4] - monopole_plaket[2][place / 4]);
  move(1, -1);
  linkt.move(1, -1);
  linkx.move(1, -1);
  linkz.move(1, -1);
  move(2, 1);
  linkt.move(2, 1);
  linkx.move(2, 1);
  linky.move(2, 1);
  j2 = -(monopole_plaket[0][linkt.place / 4] - monopole_plaket[0][place / 4]) -
       (monopole_plaket[4][linkx.place / 4] - monopole_plaket[4][place / 4]) +
       monopole_plaket[2][linky.place / 4] - monopole_plaket[2][place / 4];
  move(2, -1);
  linkt.move(2, -1);
  linkx.move(2, -1);
  linky.move(2, -1);

  J[0] = j0 / 2 / Pi;
  J[1] = j1 / 2 / Pi;
  J[2] = j2 / 2 / Pi;
  J[3] = j3 / 2 / Pi;
}

// specializations

// su2
template su2 link1::plaket_mu(const std::vector<su2> &array, int mu);
template su2 link1::plaket_schwinger_average(const std::vector<su2> &array,
                                             int mu);
template su2 link1::schwinger_line(const std::vector<su2> &array, int d,
                                   int dir, int x);
template su2 link1::polyakov_loop(const std::vector<su2> &array);
template su2 link1::wilson_loop(const std::vector<su2> &array, int r, int t);
template su2 link1::wilson_loop_schwinger(const std::vector<su2> &array, int r,
                                          int t);
template su2 link1::wilson_line(const std::vector<su2> &array, int length);
template su2 link1::wilson_line_single(const std::vector<su2> &array,
                                       int length);
template su2 link1::wilson_line_offaxis(const std::vector<su2> &array,
                                        const std::vector<int> &pattern);
template double
link1::field1(const std::vector<std::vector<su2>> &schwinger_line,
              const std::vector<su2> &plaket,
              const std::vector<su2> &polyakov_loop, int d, int D, int dir,
              int x);
template double link1::field2(const std::vector<su2> &plaket,
                              const std::vector<su2> &polyakov_loop, int d,
                              int D, int dir, int x);
template double link1::field3(const std::vector<su2> &polyakov_loop, int D,
                              int x);

// abelian
template abelian link1::plaket_mu(const std::vector<abelian> &array, int mu);
template abelian
link1::plaket_schwinger_average(const std::vector<abelian> &array, int mu);
template abelian link1::schwinger_line(const std::vector<abelian> &array, int d,
                                       int dir, int x);
template abelian link1::polyakov_loop(const std::vector<abelian> &array);
template abelian link1::wilson_loop(const std::vector<abelian> &array, int r,
                                    int t);
template abelian link1::wilson_loop_schwinger(const std::vector<abelian> &array,
                                              int r, int t);
template abelian link1::wilson_line(const std::vector<abelian> &array,
                                    int length);
template abelian link1::wilson_line_single(const std::vector<abelian> &array,
                                           int length);
template abelian link1::wilson_line_offaxis(const std::vector<abelian> &array,
                                            const std::vector<int> &pattern);
template double
link1::field1(const std::vector<std::vector<abelian>> &schwinger_line,
              const std::vector<abelian> &plaket,
              const std::vector<abelian> &polyakov_loop, int d, int D, int dir,
              int x);
template double link1::field2(const std::vector<abelian> &plaket,
                              const std::vector<abelian> &polyakov_loop, int d,
                              int D, int dir, int x);
template double link1::field3(const std::vector<abelian> &polyakov_loop, int D,
                              int x);

// su3
template su3 link1::plaket_mu(const std::vector<su3> &array, int mu);
template su3 link1::plaket_schwinger_average(const std::vector<su3> &array,
                                             int mu);
template su3 link1::schwinger_line(const std::vector<su3> &array, int d,
                                   int dir, int x);
template su3 link1::polyakov_loop(const std::vector<su3> &array);
template su3 link1::wilson_loop(const std::vector<su3> &array, int r, int t);
template su3 link1::wilson_loop_schwinger(const std::vector<su3> &array, int r,
                                          int t);
template su3 link1::wilson_line(const std::vector<su3> &array, int length);
template su3 link1::wilson_line_single(const std::vector<su3> &array,
                                       int length);
template su3 link1::wilson_line_offaxis(const std::vector<su3> &array,
                                        const std::vector<int> &pattern);
template double
link1::field1(const std::vector<std::vector<su3>> &schwinger_line,
              const std::vector<su3> &plaket,
              const std::vector<su3> &polyakov_loop, int d, int D, int dir,
              int x);
template double link1::field2(const std::vector<su3> &plaket,
                              const std::vector<su3> &polyakov_loop, int d,
                              int D, int dir, int x);
template double link1::field3(const std::vector<su3> &polyakov_loop, int D,
                              int x);