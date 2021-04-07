#pragma once

using namespace std;

extern int x_size;
extern int y_size;
extern int z_size;
extern int t_size;

// size of floating point numbers used in code
#ifdef DOUBLE
#define FLOAT double
#else
#define FLOAT float
#endif

#include "data.h"
#include "math.h"
#include "matrix.h"
#include "stdlib.h"
#include <iostream>
#include <unordered_map>
#include <vector>

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
  int multiplier[4] = {4, x_size * 4, x_size *y_size * 4,
                       x_size *y_size *z_size * 4};

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

  // changes direction to dir
  void move_dir(int dir);

  // gets matrix on current link, only gets matrices in positive directions
  // in order to get conjugated matrix use .conj() or ^ operator in matrix
  // returns pointer for better performance
  template <class T> const T *get_matrix(const vector<T> &array);

  // calculates plaket matrix in current direction and mu plane
  // only positive directions make sense
  template <class T> T plaket_mu(const vector<T> &array, int mu);

  // calculate schwinger line matrix in current direction and dir plane of
  template <class T>
  T schwinger_line(const vector<T> &array, int d, int dir, int x);

  // calculates polyakov loop in currect direction
  template <class T> T polyakov_loop(const vector<T> &array);

  // calculate wilson loop of r*t size in current direction and 4-th direction
  // plane
  template <class T> T wilson_loop(const vector<T> &array, int r, int t);

  // used in wilson_loop
  // calculates straight lines for wilson loop in advance
  template <class T> T wilson_line(const vector<T> &array, int length);

  template <class T>
  T wilson_line_offaxis(const vector<T> &array, const vector<int> &pattern);

  // d - distance between "left" source and plaket
  // D - distance between sources
  template <class T>
  FLOAT field1(
      const vector<vector<T>> &schwinger_line, const vector<T> &plaket,
      const vector<T> &polyakov_loop, int d, int D, int dir,
      int x); // Link is attached to the "left" source, dir points to the plaket
  // second numerator
  template <class T>
  FLOAT
  field2(const vector<T> &plaket, const vector<T> &polyakov_loop, int d, int D,
         int dir,
         int x); // attached to the "left" source, dir points to the plaket
  // denominator
  template <class T>
  FLOAT field3(const vector<T> &polyakov_loop, int D,
               int x); // attached to the "left" source and points to another

  // monopoles
  FLOAT monopole_plaket_mu(vector<FLOAT> &angles, int mu);
  int current_test(FLOAT *J);
  void get_current(vector<vector<FLOAT>> &monopole_plaket, FLOAT *J);
};

ostream &operator<<(ostream &os, const link1 &link);