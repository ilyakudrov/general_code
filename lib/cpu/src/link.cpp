#include "../include/link.h"
#include "../include/data.h"
#include <cmath>

#define Pi 3.141592653589793238462643383279502884
#define data_size                                                              \
  4 * lattice_size[0] * lattice_size[1] * lattice_size[2] * lattice_size[3]
#define PLACE3_LINK_DIR                                                        \
  (coordinate[3]) * 3 * lattice_size[0] * lattice_size[1] * lattice_size[2] +  \
      (coordinate[2]) * 3 * lattice_size[0] * lattice_size[1] +                \
      (coordinate[1]) * 3 * lattice_size[0] + (coordinate[0]) * 3 +            \
      abs(direction) - 1
#define PLACE4_LINK_DIR                                                        \
  (coordinate[3]) * 4 * lattice_size[0] * lattice_size[1] * lattice_size[2] +  \
      (coordinate[2]) * 4 * lattice_size[0] * lattice_size[1] +                \
      (coordinate[1]) * 4 * lattice_size[0] + (coordinate[0]) * 4 +            \
      abs(direction) - 1
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
  direction = 1;
  for (int i = 0; i < 4; i++) {
    coordinate[i] = 0;
  }
}

ostream &operator<<(ostream &os, const link1 &link) {
  os << "x: " << link.coordinate[0] << " y: " << link.coordinate[1]
     << " z: " << link.coordinate[2] << " t: " << link.coordinate[3]
     << " dir: " << link.direction << endl;
  return os;
}

void link1::move(int dir, int step) {
  coordinate[abs(dir) - 1] +=
      (step * (dir / abs(dir)) + lattice_size[abs(dir) - 1]);
  coordinate[abs(dir) - 1] =
      coordinate[abs(dir) - 1] % lattice_size[abs(dir) - 1];
}
void link1::go(int x, int y, int z, int t) {
  coordinate[0] = x;
  coordinate[1] = y;
  coordinate[2] = z;
  coordinate[3] = t;
}
void link1::move_dir(int dir) { direction = dir; }
int link1::get_place() { return PLACE4_LINK_DIR; }
int link1::get_place1() { return PLACE1_LINK_NODIR; }

template <class T> T link1::get_matrix(const vector<T> &vec) {
  T A;
  if (direction / abs(direction) == 1) {
    return vec[PLACE4_LINK_DIR];
  }
  if (direction / abs(direction) == -1) {
    move(direction, 1);
    A = vec[PLACE4_LINK_DIR];
    move(direction, -1);
    return A.conj();
  } else {
    cout << "error in get_matrix" << endl;
    return A;
  }
}

template <class T> T link1::plaket_mu(const vector<T> &array, int mu) {
  int dir = direction;
  T A = get_matrix(array);
  move(dir, 1);
  move_dir(-mu); // the same temporal direction as polyakov loop, connected to
  // schwinger line, has
  A = A * get_matrix(array);
  move(-mu, 1);
  move_dir(-dir);
  A = A * get_matrix(array);
  move(-dir, 1);
  move_dir(mu);
  A = A * get_matrix(array);
  move(mu, 1);
  move_dir(dir);
  return A;
}

template <class T>
T link1::schwinger_line(const vector<T> &array, int d, int dir, int x) {
  int dir1 = direction;
  T A;
  for (int i = 0; i < d; i++) {
    A = A * get_matrix(array);
    move(dir1, 1);
  }
  move_dir(dir);
  for (int i = 0; i < x; i++) {
    A = A * get_matrix(array);
    move(dir, 1);
  }
  move_dir(dir1);
  move(-dir, x);
  move(-dir1, d);
  return A;
}

template <class T> T link1::polyakov_loop(const vector<T> &array) {
  T A;
  for (int i = 0; i < lattice_size[3]; i++) {
    A = A * get_matrix(array);
    move(direction, 1);
  }
  return A;
}

template <class T> T link1::wilson_loop(const vector<T> &array, int r, int t) {
  int dir = direction;
  T A;
  for (int i = 0; i < r; i++) {
    A = A * get_matrix(array);
    move(dir, 1);
  }
  move_dir(4);
  for (int i = 0; i < t; i++) {
    A = A * get_matrix(array);
    move(4, 1);
  }
  move_dir(-dir);
  for (int i = 0; i < r; i++) {
    A = A * get_matrix(array);
    move(-dir, 1);
  }
  move_dir(-4);
  for (int i = 0; i < t; i++) {
    A = A * get_matrix(array);
    move(-4, 1);
  }
  move_dir(dir);
  return A;
}

template <class T> T link1::wilson_line(const vector<T> &array, int length) {
  int dir = direction;
  T A;
  for (int i = 0; i < length; i++) {
    A = A * get_matrix(array);
    move(dir, 1);
  }
  return A;
}

template <class T>
FLOAT link1::field1(const vector<vector<T>> &schwinger_line,
                    const vector<T> &plaket, const vector<T> &polyakov_loop,
                    int d, int D, int dir, int x) {
  int dir1 = direction;
  T C = schwinger_line[dir - 1][PLACE3_LINK_DIR];
  move(dir1, d);
  move(dir, x);
  T B = plaket[PLACE3_LINK_DIR];
  move(-dir, x);
  move(-dir1, d);
  move_dir(-4);
  T A = polyakov_loop[PLACE1_LINK_NODIR];
  A = A.conj() * C * B;
  A = A * C.conj();
  move(dir1, D);
  move_dir(4);
  B = polyakov_loop[PLACE1_LINK_NODIR];
  move(-dir1, D);
  move_dir(dir1);
  return A.tr() * B.tr();
}

template <class T>
FLOAT link1::field2(const vector<T> &plaket, const vector<T> &polyakov_loop,
                    int d, int D, int dir, int x) {
  int dir1 = direction;
  move_dir(-4);
  T A = polyakov_loop[PLACE1_LINK_NODIR];
  move(dir1, d);
  move(dir, x);
  move_dir(dir1);
  T B = plaket[PLACE3_LINK_DIR];
  move(-dir, x);
  move(dir1, D - d);
  move_dir(4);
  T C = polyakov_loop[PLACE1_LINK_NODIR];
  move(-dir1, D);
  move_dir(dir1);
  return B.tr() * C.tr() * A.conj().tr();
}

template <class T>
FLOAT link1::field3(const vector<T> &polyakov_loop, int D, int x) {
  int dir1 = direction;
  move_dir(-4);
  T A = polyakov_loop[PLACE1_LINK_NODIR];
  move(dir1, D);
  move_dir(4);
  T B = polyakov_loop[PLACE1_LINK_NODIR];
  move(-dir1, D);
  move_dir(dir1);
  return A.conj().tr() * B.tr();
}

// monopoles

/*template <class T>
FLOAT link1::monopole_plaket(data<T> &conf, int i, int j) {
  int dir = direction;
  move_dir(i);
  float angle = get_angle_abelian(conf.array);
  move(direction, 1);
  move_dir(j);
  angle += get_angle_abelian(conf.array);
  move(direction, 1);
  move_dir(-i);
  angle += get_angle_abelian(conf.array);
  move(direction, 1);
  move_dir(-j);
  angle += get_angle_abelian(conf.array);
  move(direction, 1);
  move_dir(dir);
  while ((angle > Pi) || (angle < -Pi)) {
    if (angle > Pi)
      angle = angle - 2 * Pi;
    if (angle < -Pi)
      angle = angle + 2 * Pi;
  }
  return angle;
}

template <class T> FLOAT link1::get_current(data<T> &conf) {
  FLOAT jj = 0.;
  link1 Lx(coordinate[0], coordinate[1], coordinate[2], coordinate[3],
           lattice_size[0], lattice_size[1], lattice_size[2], lattice_size[3]);
  Lx.move(1, 1);
  link1 Ly(coordinate[0], coordinate[1], coordinate[2], coordinate[3],
           lattice_size[0], lattice_size[1], lattice_size[2], lattice_size[3]);
  Ly.move(2, 1);
  link1 Lz(coordinate[0], coordinate[1], coordinate[2], coordinate[3],
           lattice_size[0], lattice_size[1], lattice_size[2], lattice_size[3]);
  Lz.move(3, 1);
  link1 Lt(coordinate[0], coordinate[1], coordinate[2], coordinate[3],
           lattice_size[0], lattice_size[1], lattice_size[2], lattice_size[3]);
  Lt.move(4, 1);
  if (direction == 4)
    jj = Lx.monopole_plaket(conf, 2, 3) - monopole_plaket(conf, 2, 3) -
         (Ly.monopole_plaket(conf, 1, 3) - monopole_plaket(conf, 1, 3)) +
         Lz.monopole_plaket(conf, 1, 2) - monopole_plaket(conf, 1, 2);
  if (direction == 1)
    jj = -(Lt.monopole_plaket(conf, 2, 3) - monopole_plaket(conf, 2, 3)) +
         (Ly.monopole_plaket(conf, 4, 3) - monopole_plaket(conf, 4, 3)) -
         (Lz.monopole_plaket(conf, 4, 2) - monopole_plaket(conf, 4, 2));
  if (direction == 2)
    jj = Lt.monopole_plaket(conf, 1, 3) - monopole_plaket(conf, 1, 3) -
         (Lx.monopole_plaket(conf, 4, 3) - monopole_plaket(conf, 4, 3)) +
         Lz.monopole_plaket(conf, 4, 1) - monopole_plaket(conf, 4, 1);
  if (direction == 3)
    jj = -(Lt.monopole_plaket(conf, 1, 2) - monopole_plaket(conf, 1, 2)) +
         (Lx.monopole_plaket(conf, 4, 2) - monopole_plaket(conf, 4, 2)) -
         (Ly.monopole_plaket(conf, 4, 1) - monopole_plaket(conf, 4, 1));
  return jj / 2. / Pi;
}

template <class T> int link1::current_test(FLOAT *J) {
  for (int i = 1; i <= 4; i++) {
    move_dir(i);
    if ((J[get_place()] > 0.3) || (J[get_place()] < -0.3)) {
      return i;
    }
  }

  for (int i = 1; i <= 4; i++) {
    move_dir(i);
    move(i, -1);
    if ((J[get_place()] > 0.3) || (J[get_place()] < -0.3)) {
      return -i;
    }
    move(i, 1);
  }
  return 0;
}*/

// specializations

// su2
template su2 link1::get_matrix(const vector<su2> &vec);
template su2 link1::plaket_mu(const vector<su2> &array, int mu);
template su2 link1::schwinger_line(const vector<su2> &array, int d, int dir,
                                   int x);
template su2 link1::polyakov_loop(const vector<su2> &array);
template su2 link1::wilson_loop(const vector<su2> &array, int r, int t);
template su2 link1::wilson_line(const vector<su2> &array, int length);
template FLOAT link1::field1(const vector<vector<su2>> &schwinger_line,
                             const vector<su2> &plaket,
                             const vector<su2> &polyakov_loop, int d, int D,
                             int dir, int x);
template FLOAT link1::field2(const vector<su2> &plaket,
                             const vector<su2> &polyakov_loop, int d, int D,
                             int dir, int x);
template FLOAT link1::field3(const vector<su2> &polyakov_loop, int D, int x);

// abelian
template abelian link1::get_matrix(const vector<abelian> &vec);
template abelian link1::plaket_mu(const vector<abelian> &array, int mu);
template abelian link1::schwinger_line(const vector<abelian> &array, int d,
                                       int dir, int x);
template abelian link1::polyakov_loop(const vector<abelian> &array);
template abelian link1::wilson_loop(const vector<abelian> &array, int r, int t);
template abelian link1::wilson_line(const vector<abelian> &array, int length);
template FLOAT link1::field1(const vector<vector<abelian>> &schwinger_line,
                             const vector<abelian> &plaket,
                             const vector<abelian> &polyakov_loop, int d, int D,
                             int dir, int x);
template FLOAT link1::field2(const vector<abelian> &plaket,
                             const vector<abelian> &polyakov_loop, int d, int D,
                             int dir, int x);
template FLOAT link1::field3(const vector<abelian> &polyakov_loop, int D,
                             int x);