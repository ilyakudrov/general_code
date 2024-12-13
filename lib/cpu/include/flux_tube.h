#pragma once

#include <map>

#include "link.h"

// Wilson_plaket_schwinger_correlator
template <class T>
std::vector<T> calculate_plaket_time_left_down(const std::vector<T> &array);
template <class T>
std::vector<T> calculate_plaket_time_left_up(const std::vector<T> &array);
template <class T>
std::vector<T> calculate_plaket_time_right_down(const std::vector<T> &array);
template <class T>
std::vector<T> calculate_plaket_time_right_up(const std::vector<T> &array);
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
std::vector<T> calculate_wilson_loops_schwinger(const std::vector<T> &array,
                                                int r, int time);
template <class T>
std::vector<T>
calculate_wilson_loops_schwinger_opposite(const std::vector<T> &array, int r,
                                          int time);
template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_longitudinal_l(
    const std::vector<T> &array,
    const std::vector<std::vector<T>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_ouside);
template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_longitudinal_tr(
    const std::vector<T> &array,
    const std::vector<std::vector<T>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_ouside);
template <class T>
double schwinger_electric_long_tr_even_test(const std::vector<T> &conf, int t,
                                            int r);
template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_transversal_l(
    const std::vector<T> &array,
    const std::vector<std::vector<T>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_max);
template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_electric_transversal_tr(
    const std::vector<T> &array,
    const std::vector<std::vector<T>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_max);
template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_longitudinal_l(
    const std::vector<T> &array,
    const std::vector<std::vector<T>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_ouside);
template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_longitudinal_tr(
    const std::vector<T> &array,
    const std::vector<std::vector<T>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_ouside);
template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_transversal_l(
    const std::vector<T> &array,
    const std::vector<std::vector<T>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_max);
template <class T>
std::map<std::tuple<int, int, int>, double>
flux_schwinger_magnetic_transversal_tr(
    const std::vector<T> &array,
    const std::vector<std::vector<T>> &schwinger_lines_short, int T_min,
    int T_max, int R_min, int R_max, int d_max);
template <class T>
void flux_schwinger_all(
    const std::vector<T> &array,
    const std::vector<std::vector<T>> &schwinger_lines_short,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_electric_long_l,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_electric_long_tr,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_electric_trans_l,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_electric_trans_tr,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_magnetic_long_l,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_magnetic_long_tr,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_magnetic_trans_l,
    std::map<std::tuple<int, int, int>, double>
        &flux_tube_schwinger_magnetic_trans_tr,
    int T_min, int T_max, int R_min, int R_max, int d_ouside, int d_max);

// Wilson_plaket_correlator
template <class T>
std::vector<double> calculate_plaket_time_trace_l(const std::vector<T> &array);
template <class T>
std::vector<double> calculate_plaket_time_trace_tr(const std::vector<T> &array);
template <class T>
std::vector<double> calculate_plaket_space_trace_l(const std::vector<T> array);
template <class T>
std::vector<double>
calculate_plaket_space_trace_tr(const std::vector<T> &array);
template <class T>
std::vector<double> calculate_plaket_space_tr(const std::vector<T> &array);
double plaket4_time(const std::vector<double> &plaket_tr, link1 &link);
double plaket4_space(const std::vector<double> &plaket_tr, link1 &link, int nu);

template <class T>
std::vector<double> calculate_wilson_loop_tr(const std::vector<T> &array, int r,
                                             int time);
template <class T>
std::map<std::tuple<int, int, int>, double>
calculate_wilson_plaket_correlator_electric_longitudinal(
    const std::vector<double> &plaket_tr, const std::vector<T> &conf_wilson,
    int T_min, int T_max, int R_min, int R_max, int d_ouside);
template <class T>
std::map<std::tuple<int, int, int>, double>
calculate_wilson_plaket_correlator_electric_transversal(
    const std::vector<double> &plaket_tr, const std::vector<T> &conf_wilson,
    int T_min, int T_max, int R_min, int R_max, int d_ouside);

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
std::vector<double> wilson_plane_tr(const std::vector<T> &wilson_lines_mu,
                                    const std::vector<T> &wilson_lines_nu,
                                    int size_mu1, int size_mu2, int size_nu1,
                                    int size_nu2, int length_mu, int length_nu);

template <class T>
std::vector<double> plaket_plane_tr(const std::vector<T> &conf_mu,
                                    const std::vector<T> &conf_nu, int size_mu1,
                                    int size_mu2, int size_nu1, int size_nu2);

void plaket_plane_aver(std::vector<double> &plaket_aver_tr,
                       const std::vector<double> &plaket_tr, int size_mu1,
                       int size_mu2, int size_nu1, int size_nu2);

template <class T>
std::vector<double>
plaket_aver_tr_time(const std::vector<std::vector<T>> &conf);

template <class T>
std::vector<double>
plaket_aver_tr_space(const std::vector<std::vector<T>> &conf);

void wilson_plaket_correlator_plane_longitudinal(
    std::vector<double> &correlator, const std::vector<double> &wilson_loop_tr,
    const std::vector<double> &plaket_tr, int size_mu1, int size_mu2,
    int x_trans, int d_min, int d_max);

void wilson_plaket_correlator_plane_transversal(
    std::vector<double> &correlator, const std::vector<double> &wilson_loop_tr,
    const std::vector<double> &plaket_tr, int size_mu1, int size_mu2,
    int size_nu1, int size_nu2, int d, int x_trans_min, int x_trans_max);

template <class T>
std::vector<double> plaket_tr_single_dir_test(const std::vector<T> &conf,
                                              int mu);

template <class T>
std::map<std::tuple<int, int, int>, double>
wilson_plaket_correlator(const std::vector<double> &plaket_tr,
                         const std::vector<std::vector<T>> &conf_wilson,
                         int T_min, int T_max, int R_min, int R_max,
                         int main_coordinate, int transverse_coordinate,
                         std::string direction);