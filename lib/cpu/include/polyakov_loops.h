#pragma once

#include "matrix.h"

#include <complex>
#include <map>
#include <vector>

template <class T> double polyakov_loop(const std::vector<T> &array);
template <class T>
double polyakov_loop_parallel(const std::vector<std::vector<T>> &array);

template <class T>
std::vector<double> calculate_polyakov_loops_tr(const std::vector<T> &array);

std::vector<std::complex<double>>
calculate_polyakov_loops_tr(const std::vector<su3> &array);
std::vector<std::complex<double>>
calculate_polyakov_loops_tr(const std::vector<su3_abelian> &array);

template <class T>
std::vector<double>
calculate_polyakov_loops_tr(const std::vector<std::vector<T>> &array);

std::vector<std::complex<double>>
calculate_polyakov_loops_tr(const std::vector<std::vector<su3>> &array);
std::vector<std::complex<double>>
calculate_polyakov_loops_tr(const std::vector<std::vector<su3_abelian>> &array);

template <class T>
std::vector<double> polyakov_loop_correlator(const std::vector<T> &conf,
                                             int D_max);
template <class T>
std::vector<double>
polyakov_loop_correlator(const std::vector<std::vector<T>> &conf, int D_max);

template <class T>
std::vector<T> calculate_polyakov_loops(const std::vector<T> &array);
template <class T>
std::vector<T>
calculate_polyakov_loops(const std::vector<std::vector<T>> &array);

template <class T>
std::map<int, double>
polyakov_loop_correlator_singlet(const std::vector<T> &conf, int D_min,
                                 int D_max);
template <class T>
std::vector<double>
polyakov_loop_correlator_singlet(const std::vector<std::vector<T>> &conf,
                                 int D_max);

std::map<double, double>
polyakov_average_directions(const std::vector<double> &correlators, int D_max);

template <class T>
std::vector<double> polyakov_loop_correlator_singlet(const std::vector<T> &conf,
                                                     int D_max);