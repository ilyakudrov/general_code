#pragma once

#include <vector>

template <class T>
std::vector<std::vector<T>> separate_wilson(std::vector<T> &conf);

template <class T>
std::vector<T> wilson_lines_test1(std::vector<T> separated, int length,
                                  int step);

template <class T>
std::vector<T> wilson_lines_test2(std::vector<T> separated, int length,
                                  int step);

template <class T>
std::vector<T> wilson_lines_increase_test(std::vector<T> separated,
                                          std::vector<T> wilson_lines_old,
                                          int length, int step);

template <class T>
std::vector<std::vector<T>> separate_wilson_unchanged(std::vector<T> &conf);

template <class T>
std::vector<T> wilson_lines_test3(std::vector<T> separated, int length,
                                  int size1, int size2);

template <class T>
double wilson_plane_test1(std::vector<T> &wilson_lines_mu,
                          std::vector<T> &wilson_lines_nu, int size_mu1,
                          int size_mu2, int size_nu1, int size_nu2,
                          int length_mu, int length_nu, int thread_num);

template <class T>
double wilson_loop_test_time(std::vector<std::vector<T>> &wilson_lines,
                             int length_R, int length_T, int thread_num);