#pragma once

#include "data.h"
#include "link.h"
#include "matrix.h"
#include <map>
#include <unordered_map>

// Plaket
template <class T> double plaket_time(const std::vector<T> &array);
template <class T> double plaket_space(const std::vector<T> &array);
template <class T> double plaket(const std::vector<T> &array);

// Wilson loop
template <class T>
std::vector<double> wilson(const std::vector<T> &array, int r_min, int r_max,
                           int time_min, int time_max);

template <class T>
std::vector<T> wilson_lines_single(const std::vector<T> &array, int length);
template <class T>
double wilson_loop_single_size(std::vector<T> lines1, std::vector<T> lines2,
                               int mu, int nu, int r1, int r2);

template <class T>
std::map<std::tuple<int, int>, double>
wilson_spatial(const std::vector<T> &array,
               std::map<std::tuple<int, int>, std::vector<T>> smeared,
               int time_min, int time_max, int r_min, int r_max);

std::vector<double> read_abelian_fortran(std::string path_abelian);

std::vector<double> read_abelian_fortran_double(std::string path_abelian);

double wilson_abelian(const std::vector<double> &array, int r, int t);

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
std::vector<wilson_result>
wilson_offaxis_adjoint(const std::vector<T> &array,
                       const std::vector<std::vector<int>> directions,
                       double r_min, double r_max, int time_min, int time_max);

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

// Polyakov_correlator
template <class T>
std::vector<double> calculate_polyakov_loops_tr(const std::vector<T> &array);

template <class T>
std::map<int, double> polyakov_loop_correlator(const std::vector<T> &conf,
                                               int D_min, int D_max);

template <class T>
std::vector<T> calculate_polyakov_loops(const std::vector<T> &array);

template <class T>
std::map<int, double>
polyakov_loop_correlator_singlet(const std::vector<T> &conf, int D_min,
                                 int D_max);

std::map<double, double>
polyakov_average_directions(std::vector<double> correlators, int D_max);

template <class T>
std::vector<double> polyakov_loop_correlator_singlet(const std::vector<T> &conf,
                                                     int D_max);

template <class T>
double plaket_plane(std::vector<T> &conf_mu, std::vector<T> &conf_nu,
                    int size_mu1, int size_mu2, int size_nu1, int size_nu2);

template <class T>
double plaket_time_parallel(std::vector<std::vector<T>> conf);

template <class T>
double plaket_space_parallel(std::vector<std::vector<T>> conf);

template <class T> double plaket_parallel(std::vector<std::vector<T>> conf);

template <class T>
std::vector<std::vector<T>> separate_wilson(std::vector<T> &conf);

template <class T>
std::vector<T> merge_wilson(std::vector<std::vector<T>> &conf_separated);

template <class T>
std::vector<T> wilson_lines(std::vector<T> separated, int length, int size1,
                            int size2);

template <class T>
double wilson_plane(std::vector<T> &wilson_lines_mu,
                    std::vector<T> &wilson_lines_nu, int size_mu1, int size_mu2,
                    int size_nu1, int size_nu2, int length_mu, int length_nu);

template <class T>
double wilson_loop_test_time(std::vector<std::vector<T>> &wilson_lines,
                             int length_R, int length_T);

template <class T>
std::map<std::tuple<int, int>, double>
wilson_parallel(std::vector<std::vector<T>> conf, int r_min, int r_max,
                int time_min, int time_max);
