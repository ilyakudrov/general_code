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
double plaket_time_test2(std::vector<std::vector<T>> &separated);

template <class T>
double plaket_plane1(std::vector<T> &conf_mu, std::vector<T> &conf_nu,
                     int size1, int size2);

template <class T>
double plaket_plane2(std::vector<T> &conf_mu, std::vector<T> &conf_nu,
                     int size1, int size2);
template <class T>
double plaket_time_test4(std::vector<std::vector<T>> &separated);
template <class T>
double plaket_time_test3(std::vector<std::vector<T>> &separated);

template <class T> double plaket_time_test5(std::vector<T> &conf);

template <class T> double plaket_time_test6(std::vector<T> &conf);

template <class T>
double plaket_plane3(std::vector<T> &conf, int mu, int nu, int size_mu1,
                     int size_mu2, int size_nu1, int size_nu2);

template <class T> double plaket_time_test7(std::vector<T> &conf);

template <class T>
double plaket_plane4(std::vector<T> &conf_mu, std::vector<T> &conf_nu,
                     int size_mu1, int size_mu2, int size_nu1, int size_nu2);

template <class T>
double plaket_time_test8(std::vector<std::vector<T>> &separated);

template <class T>
double plaket_plane5(std::vector<T> &conf_mu, std::vector<T> &conf_nu,
                     int size_mu1, int size_mu2, int size_nu1, int size_nu2,
                     int step);
template <class T>
double plaket_time_test9(std::vector<std::vector<T>> &separated, int step);

template <class T>
double plaket_plane6(std::vector<T> &conf_mu, std::vector<T> &conf_nu,
                     int size_mu1, int size_mu2, int size_nu1, int size_nu2,
                     int thread_num);
template <class T>
double plaket_time_test10(std::vector<std::vector<T>> &separated,
                          int thread_num);

template <class T> std::vector<std::vector<T>> separate_3(std::vector<T> &conf);
