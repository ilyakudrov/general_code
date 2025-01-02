#include "../include/matrix.h"
#include <map>
#include <vector>

double plaket_indexed(const std::vector<su3> &conf);
template <class T>
std::map<std::tuple<int, int>, double>
wilson_parallel_indexed(const std::vector<T> &conf, int r_min, int r_max,
                        int time_min, int time_max);
template <class T>
std::map<std::tuple<int, int>, double>
wilson_parallel_indexed_single_rxt(const std::vector<T> &conf, int r_min,
                                   int r_max, int time_min, int time_max);