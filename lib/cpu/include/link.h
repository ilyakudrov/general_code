#pragma once

#include "math.h"
#include "matrix.h"
#include "stdlib.h"

#include <iostream>
#include <vector>

extern int x_size;
extern int y_size;
extern int z_size;
extern int t_size;

// Class link1 goes on a periodic 4D lattice and calculates some obseervables
// can take matrices from vector of matrices in usual order (see data.h)
// numeration of directions is as follows
// x - 0, y - 1, z - 2, t - 3 (in arrays numeratioon starts with 0)
class link1 {
public:
  // lattice sizes
  int lattice_size[4];
  // coordinates of the site
  // coordinate[i] = 0..lattise_size[i]-1
  int coordinate[4];
  // holds previous values of coordinates
  int coordinate_old[4];
  // current direction of the link
  // only positive values make sense
  int direction;
  // holds a place of matrix in array at direction 0 and current lattice site
  int place;
  // holds distance to move in array after moving on the lattice in each
  // direction
  int multiplier[4];

  link1(int lattice_size_x, int lattice_size_y, int lattice_size_z,
        int lattice_size_t);

  link1(const link1 &link);

  link1();

  // moves link in direction dir on step steps and updates place and
  // coordinates_old
  // dir is positive, step may be negative
  void move(int dir, int step);

  // update coordinates after some changes
  void update(int dir);

  // goes on site with following coordinates
  void go(int x, int y, int z, int t);
  void go(std::vector<int> &coordinate);

  // changes direction to dir
  void move_dir(int dir);

  // go and update
  void go_update(int x, int y, int z, int t);
  void go_update(std::vector<int> &coordinate);

  const spin *get_spin(const std::vector<spin> &vec);

  const spin *get_consecutive_spin(const std::vector<spin> &vec, int mu);

  template <class T> const T *get_matrix(const std::vector<T> &array);

  // calculates plaket matrix in current direction and mu plane
  // end of the name means orientation of starting/end point of plaket
  template <class T> T plaket_left_down(const std::vector<T> &array, int mu);
  template <class T> T plaket_left_up(const std::vector<T> &array, int mu);
  template <class T> T plaket_right_down(const std::vector<T> &array, int mu);
  template <class T> T plaket_right_up(const std::vector<T> &array, int mu);

  // calculates polyakov loop in currect direction
  template <class T> T polyakov_loop(const std::vector<T> &array);
  template <class T> T polyakov_loop(const std::vector<std::vector<T>> &array);

  // calculate wilson loop of r*t size in current direction and 4-th direction
  // plane
  template <class T> T wilson_loop(const std::vector<T> &array, int r, int t);

  double wilson_loop_abelian(const std::vector<double> &array, int r, int t);

  template <class T>
  T wilson_loop_schwinger(const std::vector<T> &array, int r, int t);

  template <class T>
  T wilson_loop_schwinger_opposite(const std::vector<T> &array, int r, int t);

  // used in wilson_loop
  // calculates straight lines for wilson loop in advance
  template <class T> T wilson_line(const std::vector<T> &array, int length);

  template <class T>
  T wilson_line_single(const std::vector<T> &array, int length);

  template <class T>
  T wilson_line_offaxis(const std::vector<T> &array,
                        const std::vector<int> &pattern);

  // monopoles
  double monopole_plaket_mu(std::vector<double> &angles, int mu);
  int monopole_plaket_singular_mu(std::vector<double> &angles, int mu);
  int current_test(double *J);
  void get_current(std::vector<std::vector<double>> &monopole_plaket,
                   double *J);

  void get_current_singular(std::vector<std::vector<int>> &monopole_plaket,
                            int *J);
};

std::ostream &operator<<(std::ostream &os, const link1 &link);