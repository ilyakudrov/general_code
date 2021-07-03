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

// Wilson_plaket_schwinger_correlator
template <class T>
std::vector<T> calculate_plaket_time(const std::vector<T> array);
template <class T>
std::vector<T> calculate_plaket_space(const std::vector<T> &array);
template <class T>
std::vector<T> calculate_polyakov_loop(const std::vector<T> &array);
template <class T>
std::vector<T> calculate_wilson_loop(const std::vector<T> &array, int r,
                                     int time);
template <class T>
std::vector<T> calculate_plaket_schwinger_time(const std::vector<T> &array);
template <class T>
std::vector<std::vector<T>>
calculate_plaket_schwinger_space(const std::vector<T> &array);
template <class T>
std::vector<std::vector<T>>
calculate_schwinger_lines_short(const std::vector<T> &array, int d);
template <class T>
std::vector<std::vector<T>>
calculate_schwinger_line(const std::vector<T> &array, int d, int x_trans);
template <class T>
std::vector<FLOAT> wilson_plaket_schwinger_electric(
    const std::vector<T> &array,
    const std::vector<std::vector<std::vector<T>>> &schwinger_lines,
    const std::vector<T> &plaket, int d_min, int d_max, int t, int r);

// Wilson_plaket_correlator
template <class T>
std::vector<FLOAT> calculate_plaket_time_tr(const std::vector<T> &array);
template <class T>
std::vector<FLOAT> calculate_plaket_space_tr(const std::vector<T> &array);
FLOAT plaket4_time(const std::vector<FLOAT> &plaket_tr, link1 &link);
FLOAT plaket4_space(const std::vector<FLOAT> &plaket_tr, link1 &link, int nu);
template <class T>
std::vector<FLOAT> calculate_wilson_loop_tr(const std::vector<T> &array, int r,
                                            int time);

std::vector<FLOAT>
wilson_plaket_correlator_electric(const std::vector<FLOAT> &wilson_loop_tr,
                                  const std::vector<FLOAT> &plaket_tr, int r,
                                  int time, int x_trans, int d_min, int d_max);
std::vector<FLOAT>
wilson_plaket_correlator_electric_x(const std::vector<FLOAT> &wilson_loop_tr,
                                    const std::vector<FLOAT> &plaket_tr, int r,
                                    int time, int x_trans_min, int x_trans_max,
                                    int d);
std::vector<FLOAT>
wilson_plaket_correlator_magnetic(const std::vector<FLOAT> &wilson_loop_tr,
                                  const std::vector<FLOAT> &plaket_tr, int r,
                                  int time, int x_trans, int d_min, int d_max);
std::vector<FLOAT>
wilson_plaket_correlator_magnetic_x(const std::vector<FLOAT> &wilson_loop_tr,
                                    const std::vector<FLOAT> &plaket_tr, int r,
                                    int time, int x_trans_min, int x_trans_max,
                                    int d);