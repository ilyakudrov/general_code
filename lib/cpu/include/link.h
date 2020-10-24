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
// x - 1, y - 2, z - 3, t - 4 (in arrays numeratioon starts with 0)
template <class T> class link1 {
public:
  // lattice sizes
  int lattice_size[4];
  // coordinates of the site
  // coordinate[i] = 0..lattise_size[i]-1
  int coordinate[4];
  // current direction of the link
  // negative values means opposite direction
  int direction;

  link1(int lattice_size_x, int lattice_size_y, int lattice_size_z,
        int lattice_size_t);

  // prints link information
  void print_link();

  // moves link in direction dir on step steps
  void move(int dir, int step);

  // goes on site with following coordinates
  void go(int x, int y, int z, int t);

  // changes direction to dir
  void move_dir(int dir); // changes the direction

  // gets place in different arrays (see link.cpp)
  int get_place();
  int get_place1();

  // gets matrix on carrent link, gets conjugate if direction is opposite
  T get_matrix(const vector<T> &vec);

  // calculates plaket matrix in current direction and mu plane
  T plaket_mu(const vector<T> &array, int mu);

  // calculates polyakov loop in currect direction
  T polyakov_loop(const vector<T> &array);

  // calculate wilson loop of r*t size in current direction and 4-th direction
  // plane
  T wilson_loop(const vector<T> &array, int r, int t);

  // used in wilson_loop
  T wilson_line(const vector<T> &array, int length);

  // d - distance between "left" source and plaket
  // D - distance between sources
  FLOAT field1(
      const vector<vector<T>> &schwinger_line, const vector<T> &plaket,
      const vector<T> &polyakov_loop, int d, int D, int dir,
      int x); // Link is attached to the "left" source, dir points to the plaket
  // second numerator
  FLOAT
  field2(const vector<T> &plaket, const vector<T> &polyakov_loop, int d, int D,
         int dir,
         int x); // attached to the "left" source, dir points to the plaket
  // denominator
  FLOAT field3(const vector<T> &polyakov_loop, int D,
               int x); // attached to the "left" source and points to another
  T staples_first(const vector<T> &vec, int eta);
  T staples_second(const vector<vector<T>> &smearing_first,
                   unordered_map<int, int> &indexes, int rho, int mu, int nu);
  T staples_second_refresh(const vector<T> &vec, int eta, int nu,
                           FLOAT alpha3); // staples for refreshing
                                          // algorythm(refresh link every step)
  T staples_third(const vector<vector<T>> &smearing_second,
                  unordered_map<int, int> indexes, int nu, int mu);
  T staples_third_refresh(const vector<T> &vec, int eta, FLOAT alpha2,
                          FLOAT alpha3);
  vector<T> smearing_first(const vector<T> &array, FLOAT alpha3, int mu, int nu,
                           int rho);
  vector<vector<T>> smearing_first_full(const vector<T> &array, FLOAT alpha3);
  vector<T> smearing_second(const vector<T> &array,
                            vector<vector<T>> &smearing_first, FLOAT alpha2,
                            int mu, int nu);
  vector<vector<T>> smearing_second_full(const vector<T> &array,
                                         vector<vector<T>> &smearing_first,
                                         FLOAT alpha2);
  vector<T> smearing_HYP(const vector<T> &array,
                         vector<vector<T>> &smearing_second, FLOAT alpha1);
  vector<T> smearing_APE(const vector<T> &array, FLOAT alpha_APE);
  T smearing_first_refresh(const vector<T> &vec, int nu, int rho,
                           FLOAT alpha3); // refresh link every step
  T smearing_second_refresh(const vector<T> &vec, int nu, FLOAT alpha2,
                            FLOAT alpha3); // refresh link every step
  vector<T> smearing_HYP_refresh(data<T> &conf, FLOAT alpha1, FLOAT alpha2,
                                 FLOAT alpha3); // refresh link every step
  vector<T> smearing_APE_refresh(data<T> &conf,
                                 FLOAT alpha_APE); // refresh link every step
  vector<T> smearing_stout(data<T> &conf, FLOAT rho);
  T stout_factor(data<T> &conf, FLOAT rho);
  T stout_omega(data<T> &conf, FLOAT rho);
  void gauge_transform(data<T> &conf);
  // monopoles
  /*FLOAT monopole_plaket(data<T> &conf, int i, int j); // monopole_plaket
  FLOAT get_current(data<T> &conf);
  int current_test(FLOAT *J);*/
};

template <class T> ostream &operator<<(ostream &os, const link1<T> &link);
void make_map_first(unordered_map<int, int> &indexes);
void make_map_second(unordered_map<int, int> &indexes);