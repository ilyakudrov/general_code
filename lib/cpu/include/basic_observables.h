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

// Plaket
template <class T> FLOAT plaket_time(const vector<T> &array);
template <class T> FLOAT plaket_space(const vector<T> &array);
template <class T> FLOAT plaket(const vector<T> &array);

// Wilson loop
template <class T>
vector<FLOAT> wilson(const vector<T> &array, int r_min, int r_max, int time_min,
                     int time_max);
template <class T>
vector<T> wilson_lines(const vector<T> &array, int mu, int length);
template <class T>
vector<T> wilson_line_increase(const vector<T> &array, const vector<T> &lines,
                               int mu, int length);

// contains offaxis wilson loop information
struct wilson_result {
  FLOAT wilson_loop;
  int time_size;
  FLOAT space_size;
  int statistics_size;
};

void wilson_offaxis_reduce(vector<wilson_result> &wilson_offaxis_result);

template <class T>
vector<wilson_result>
wilson_offaxis(const vector<T> &array, const vector<vector<int>> directions,
               FLOAT r_min, FLOAT r_max, int time_min, int time_max);

vector<vector<int>> generate_permutations(const vector<int> &direction);

vector<vector<int>> generate_reflections(const vector<int> &direction);

template <class T>
FLOAT calculate_wilson_loop_offaxis(const vector<T> &time_lines, int time,
                                    const vector<T> &space_lines,
                                    const vector<int> &direction);

template <class T>
vector<T> wilson_lines_offaxis(const vector<T> &array,
                               const vector<int> pattern);

template <class T>
vector<T> wilson_lines_offaxis_increase(const vector<T> &array,
                                        const vector<T> &lines1,
                                        const vector<int> pattern,
                                        const vector<int> direction);

vector<vector<int>> generate_directions(int length_max);

bool is_direction_same(const vector<int> &direction1,
                       const vector<int> &direction2);

vector<int> permutate_pattern(const vector<int> &pattern,
                              const vector<int> permutation);

vector<int> reflect_pattern(const vector<int> &pattern,
                            const vector<int> reflection);

vector<int> make_offaxis_pattern(const vector<int> &line_direction);

// Polyakov loop
template <class T> FLOAT polyakov(const vector<T> &array);

// Polyakov_correlator
template <class T> FLOAT polyakov_loop_corelator(const vector<T> &array, int D);

FLOAT MAG_functional_su2(const vector<su2> &array);