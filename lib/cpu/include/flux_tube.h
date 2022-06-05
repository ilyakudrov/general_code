#pragma once

#include "data.h"
#include "link.h"
#include "matrix.h"

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
std::vector<T> calculate_schwinger_lines_short(const std::vector<T> &array,
                                               int d);
template <class T>
std::vector<std::vector<T>>
calculate_schwinger_line(const std::vector<T> &array, int d, int x_trans);
template <class T>
std::map<int, double>
wilson_plaket_schwinger_electric(const std::vector<T> &array,
                                 const std::vector<T> &plaket, int d_min,
                                 int d_max, int t, int r);

// Wilson_plaket_correlator
template <class T>
std::vector<double> calculate_plaket_time_tr(const std::vector<T> &array);
template <class T>
std::vector<double> calculate_plaket_space_tr(const std::vector<T> &array);
double plaket4_time(const std::vector<double> &plaket_tr, link1 &link);
double plaket4_space(const std::vector<double> &plaket_tr, link1 &link, int nu);
template <class T>
std::vector<double> calculate_wilson_loop_tr(const std::vector<T> &array, int r,
                                             int time);

std::map<int, double>
wilson_plaket_correlator_electric(const std::vector<double> &wilson_loop_tr,
                                  const std::vector<double> &plaket_tr, int r,
                                  int time, int x_trans, int d_min, int d_max);
std::map<int, double>
wilson_plaket_correlator_electric_x(const std::vector<double> &wilson_loop_tr,
                                    const std::vector<double> &plaket_tr, int r,
                                    int time, int x_trans_max, int d);
std::map<int, double>
wilson_plaket_correlator_magnetic(const std::vector<double> &wilson_loop_tr,
                                  const std::vector<double> &plaket_tr, int r,
                                  int time, int x_trans, int d_min, int d_max);
std::map<int, double>
wilson_plaket_correlator_magnetic_x(const std::vector<double> &wilson_loop_tr,
                                    const std::vector<double> &plaket_tr, int r,
                                    int time, int x_trans_max, int d);

template <class T>
std::vector<double> wilson_plane_tr(std::vector<T> &wilson_lines_mu,
                                    std::vector<T> &wilson_lines_nu,
                                    int size_mu1, int size_mu2, int size_nu1,
                                    int size_nu2, int length_mu, int length_nu);

template <class T>
std::vector<double> plaket_plane_tr(std::vector<T> &conf_mu,
                                    std::vector<T> &conf_nu, int size_mu1,
                                    int size_mu2, int size_nu1, int size_nu2);

void plaket_plane_aver(std::vector<double> &plaket_aver_tr,
                       std::vector<double> &plaket_tr, int size_mu1,
                       int size_mu2, int size_nu1, int size_nu2);

template <class T>
std::vector<double> plaket_aver_tr_time(std::vector<std::vector<T>> conf);

template <class T>
std::vector<double> plaket_aver_tr_space(std::vector<std::vector<T>> conf);

void wilson_plaket_correlator_plane_longitudinal(
    std::vector<double> &correlator, const std::vector<double> &wilson_loop_tr,
    const std::vector<double> &plaket_tr, int size_mu1, int size_mu2,
    int x_trans, int d_min, int d_max);

void wilson_plaket_correlator_plane_transversal(
    std::vector<double> &correlator, const std::vector<double> &wilson_loop_tr,
    const std::vector<double> &plaket_tr, int size_mu1, int size_mu2,
    int size_nu1, int size_nu2, int d, int x_trans_min, int x_trans_max);

template <class T>
std::map<std::tuple<int, int, int>, double>
wilson_plaket_correlator(std::vector<std::vector<T>> &conf, int T_min,
                         int T_max, int R_min, int R_max, int main_coordinate,
                         int transverse_coordinate, std::string field,
                         std::string direction);