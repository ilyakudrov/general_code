#pragma once

#ifdef DOUBLE
#define FLOAT double
#else
#define FLOAT float
#endif

#include "data.h"
#include "link.h"
#include "matrix.h"
#include "result.h"

using namespace std;

// Wilson_plaket_schwinger_correlator
template <class T> vector<T> calculate_plaket_time(const vector<T> array);
template <class T> vector<T> calculate_plaket_space(const vector<T> &array);
template <class T> vector<T> calculate_polyakov_loop(const vector<T> &array);
template <class T>
vector<T> calculate_wilson_loop(const vector<T> &array, int r, int time);
template <class T>
vector<vector<T>> calculate_schwinger_line(const vector<T> &array, int d,
                                           int x_trans);
template <class T>
void wilson_plaket_schwinger_electric(const vector<vector<T>> &schwinger_line,
                                      const vector<T> &plaket,
                                      const vector<T> &polyakov_loop,
                                      vector<vector<result>> &field1,
                                      vector<vector<result>> &field2,
                                      vector<result> &field3, int d, int D,
                                      int x_trans);

// Wilson_plaket_correlator
template <class T>
vector<FLOAT> calculate_plaket_time_tr(const vector<T> &array);
template <class T>
vector<FLOAT> calculate_plaket_space_tr(const vector<T> &array);
FLOAT plaket4_time(const vector<FLOAT> &plaket_tr, link1 &link);
FLOAT plaket4_space(const vector<FLOAT> &plaket_tr, link1 &link, int nu);
template <class T>
vector<FLOAT> calculate_wilson_loop_tr(const vector<T> &array, int r, int time);

vector<FLOAT>
wilson_plaket_correlator_electric(const vector<FLOAT> &wilson_loop_tr,
                                  const vector<FLOAT> &plaket_tr, int r,
                                  int time, int x_trans, int d_min, int d_max);
vector<FLOAT> wilson_plaket_correlator_electric_x(
    const vector<FLOAT> &wilson_loop_tr, const vector<FLOAT> &plaket_tr, int r,
    int time, int x_trans_min, int x_trans_max, int d);
vector<FLOAT>
wilson_plaket_correlator_magnetic(const vector<FLOAT> &wilson_loop_tr,
                                  const vector<FLOAT> &plaket_tr, int r,
                                  int time, int x_trans, int d_min, int d_max);
vector<FLOAT> wilson_plaket_correlator_magnetic_x(
    const vector<FLOAT> &wilson_loop_tr, const vector<FLOAT> &plaket_tr, int r,
    int time, int x_trans_min, int x_trans_max, int d);