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

// Plaket
template <class T> FLOAT plaket_time(const std::vector<T> &array);
template <class T> FLOAT plaket_space(const std::vector<T> &array);
template <class T> FLOAT plaket(const std::vector<T> &array);

// Wilson loop
template <class T>
std::vector<FLOAT> wilson(const std::vector<T> &array, int r_min, int r_max,
                          int time_min, int time_max);
template <class T>
std::vector<T> wilson_lines(const std::vector<T> &array, int mu, int length);
template <class T>
std::vector<T> wilson_line_increase(const std::vector<T> &array,
                                    const std::vector<T> &lines, int mu,
                                    int length);

// contains offaxis wilson loop information
struct wilson_result {
  FLOAT wilson_loop;
  int time_size;
  FLOAT space_size;
  int statistics_size;
};

void wilson_offaxis_reduce(std::vector<wilson_result> &wilson_offaxis_result);

template <class T>
std::vector<wilson_result>
wilson_offaxis(const std::vector<T> &array,
               const std::vector<std::vector<int>> directions, FLOAT r_min,
               FLOAT r_max, int time_min, int time_max);

std::vector<std::vector<int>>
generate_permutations(const std::vector<int> &direction);

std::vector<std::vector<int>>
generate_reflections(const std::vector<int> &direction);

template <class T>
FLOAT calculate_wilson_loop_offaxis(const std::vector<T> &time_lines, int time,
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
template <class T> FLOAT polyakov(const std::vector<T> &array);

// Polyakov_correlator
template <class T>
FLOAT polyakov_loop_corelator(const std::vector<T> &array, int D);

FLOAT MAG_functional_su2(const std::vector<su2> &array);