#pragma once

using namespace std;

extern int x_size;
extern int y_size;
extern int z_size;
extern int t_size;

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
#include <vector>

template <class T> class link1 {
public:
  int lattice_size[4];
  int coordinate[4]; // 0..lattise_size[i]-1
  int direction;

  link1(int lattice_size_x, int lattice_size_y, int lattice_size_z,
        int lattice_size_t);

  void print_link();
  void move(int dir, int step);
  void go(int x, int y, int z, int t);
  void move_dir(int dir); // changes the direction
  int get_place();
  int get_place1();
  T get_matrix(const vector<T> &vec); // works with negative
                                      // directions(takes inverse matrix)
  FLOAT border_sign(int mu);
  FLOAT get_angle_abelian(const vector<T> &vec);
  T schwinger_line(
      const vector<T> &array, int d, int dir,
      int x); // link is attached to "left" source and directed to the plaket
  T plaket_mu(const vector<T> &array,
              int mu); // mu is the second direction
  T plaket_average4(const vector<T> &array, int mu);
  T polyakov_loop(const vector<T> &array); // attached to where it loops and
                                           // directed to the temporal direction
  T wilson_loop(const vector<T> &array, int r, int t);
  T wilson_line(const vector<T> &array, int length);
  // first numeratoк
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
  T staples_second(const vector<vector<T>> &smearing_first, int eta, int nu);
  T staples_second_refresh(const vector<T> &vec, int eta, int nu,
                           FLOAT alpha3); // staples for refreshing
                                          // algorythm(refresh link every step)
  T staples_third(const vector<vector<T>> &smearing_second, int eta);
  T staples_third_refresh(const vector<T> &vec, int eta, FLOAT alpha2,
                          FLOAT alpha3);
  vector<T> smearing_first(const vector<T> &array, FLOAT alpha3, int nu,
                           int rho);
  vector<vector<T>> smearing_first_full(const vector<T> &array, FLOAT alpha3);
  vector<T> smearing_second(const vector<T> &array,
                            vector<vector<T>> &smearing_first, FLOAT alpha2,
                            int nu);
  vector<vector<T>> smearing_second_full(const vector<T> &array,
                                         vector<vector<T>> &smearing_first,
                                         FLOAT alpha2);
  vector<T> smearing_HYP(const vector<T> &array,
                         vector<vector<T>> &smearing_second, FLOAT alpha1);
  inline int position_first(int a, int b);
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