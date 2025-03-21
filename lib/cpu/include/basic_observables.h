#pragma once

#include "matrix.h"

#include <algorithm>
#include <map>

template <class T> double plaket(const std::vector<T> &conf);
template <class T> double plaket_time(const std::vector<T> &conf);
template <class T> double plaket_space(const std::vector<T> &conf);
template <class T>
std::vector<T> wilson_lines(const std::vector<T> &array, int mu, int length);
template <class T>
std::vector<T> wilson_line_increase(const std::vector<T> &array,
                                    const std::vector<T> &lines, int mu,
                                    int length);

// contains offaxis wilson loop information
struct wilson_result {
  double wilson_loop;
  int time_size;
  double space_size;
  int statistics_size;
};

void wilson_offaxis_reduce(std::vector<wilson_result> &wilson_offaxis_result);

template <class T>
std::vector<wilson_result>
wilson_offaxis(const std::vector<T> &array,
               const std::vector<std::vector<int>> directions, double r_min,
               double r_max, int time_min, int time_max);

template <class T>
std::map<std::tuple<int, double>, double>
wilson_offaxis_result(const std::vector<T> &array, double r_min, double r_max,
                      int time_min, int time_max);

template <class T>
std::vector<wilson_result>
wilson_offaxis_adjoint(const std::vector<T> &array,
                       const std::vector<std::vector<int>> directions,
                       double r_min, double r_max, int time_min, int time_max);

template <class T>
std::map<std::tuple<int, double>, double>
wilson_offaxis_adjoint_result(const std::vector<T> &array, double r_min,
                              double r_max, int time_min, int time_max);

std::vector<std::vector<int>>
generate_permutations(const std::vector<int> &direction);

std::vector<std::vector<int>>
generate_reflections(const std::vector<int> &direction);

template <class T>
double calculate_wilson_loop_offaxis(const std::vector<T> &time_lines, int time,
                                     const std::vector<T> &space_lines,
                                     const std::vector<int> &direction);

template <class T>
double calculate_wilson_loop_offaxis_adjoint(const std::vector<T> &time_lines,
                                             int time,
                                             const std::vector<T> &space_lines,
                                             const std::vector<int> &direction);

template <class T>
std::vector<T> wilson_lines_offaxis(const std::vector<T> &array,
                                    const std::vector<int> pattern);

template <class T>
std::vector<T> wilson_lines_offaxis_increase(const std::vector<T> &array,
                                             const std::vector<T> &lines1,
                                             const std::vector<int> pattern,
                                             const std::vector<int> direction);

std::vector<std::vector<int>> generate_directions(int length_max);

bool is_direction_same(const std::vector<int> &direction1,
                       const std::vector<int> &direction2);

std::vector<int> permutate_pattern(const std::vector<int> &pattern,
                                   const std::vector<int> permutation);

std::vector<int> reflect_pattern(const std::vector<int> &pattern,
                                 const std::vector<int> reflection);

std::vector<int> make_offaxis_pattern(const std::vector<int> &line_direction);

// Polyakov loop
template <class T> double polyakov_loop(const std::vector<T> &array);
template <class T>
double polyakov_loop_parallel(const std::vector<std::vector<T>> &array);

// Polyakov_correlator
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

template <class T>
std::vector<std::vector<T>> separate_wilson(const std::vector<T> &conf);

template <class T>
std::vector<T> merge_wilson(const std::vector<std::vector<T>> &conf_separated);

template <class T>
std::vector<T> wilson_lines(const std::vector<T> &separated, int length,
                            int size1, int size2);

template <class T>
double wilson_plane(const std::vector<T> &wilson_lines_mu,
                    const std::vector<T> &wilson_lines_nu, int size_mu1,
                    int size_mu2, int size_nu1, int size_nu2, int length_mu,
                    int length_nu);

template <class T>
double wilson_loop_test_time(const std::vector<std::vector<T>> &wilson_lines,
                             int length_R, int length_T);
template <class T>
std::map<std::tuple<int, int>, double> wilson_loop(const std::vector<T> &conf,
                                                   int r_min, int r_max,
                                                   int time_min, int time_max);
template <class T>
std::map<std::tuple<int, int>, double>
wilson_loop_adjoint(const std::vector<T> &conf, int r_min, int r_max,
                    int time_min, int time_max);
template <class T>
std::map<std::tuple<int, int>, double>
wilson_gevp_indexed(const std::vector<T> &conf1, const std::vector<T> &conf2,
                    int r_min, int r_max, int time_min, int time_max);
template <class T>
std::map<std::tuple<int, int>, double>
wilson_gevp_adjoint_indexed(const std::vector<T> &conf1,
                            const std::vector<T> &conf2, int r_min, int r_max,
                            int time_min, int time_max);
template <class T>
std::map<std::tuple<int, int>, double>
wilson_gevp_parallel(const std::vector<std::vector<T>> &conf1,
                     const std::vector<std::vector<T>> &conf2, int r_min,
                     int r_max, int time_min, int time_max);

template <class T>
std::map<std::tuple<int, int>, double>
wilson_gevp_adjoint_parallel(const std::vector<std::vector<T>> &conf1,
                             const std::vector<std::vector<T>> &conf2,
                             int r_min, int r_max, int time_min, int time_max);

template <class T>
double wilson_adjoint_plane(const std::vector<T> &wilson_lines_mu,
                            const std::vector<T> &wilson_lines_nu, int size_mu1,
                            int size_mu2, int size_nu1, int size_nu2,
                            int length_mu, int length_nu);

template <class T>
std::map<std::tuple<int, int>, double>
wilson_adjoint_parallel(const std::vector<std::vector<T>> &conf, int r_min,
                        int r_max, int time_min, int time_max);

template <class T>
std::map<std::tuple<int, int, int>, double>
wilson_spatial_3d_parallel(const std::vector<std::vector<T>> &conf, int r_min,
                           int r_max, int time_min, int time_max, double alpha,
                           int smearing_start, int smearing_end,
                           int smearing_step);

template <class T>
std::map<std::tuple<int, int, int>, double>
wilson_spatial_3d_indexed(const std::vector<T> &conf, int r_min, int r_max,
                          int time_min, int time_max, double alpha,
                          int smearing_start, int smearing_end,
                          int smearing_step);
