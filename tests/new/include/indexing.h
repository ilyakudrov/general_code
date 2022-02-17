#pragma once

#include <vector>

template <class T>
double plaket_plane(std::vector<T> &conf_mu, std::vector<T> &conf_nu, int size1,
                    int size2);

template <class T>
double plaket_time_test(std::vector<std::vector<T>> &separated);
template <class T>
double plaket_time_test1(std::vector<std::vector<T>> &separated);

template <class T>
std::vector<std::vector<T>> separate_3(std::vector<T> &conf, int mu);
