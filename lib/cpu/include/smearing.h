#pragma once

#include "link.h"
#include <iostream>
#include <unordered_map>
#include <vector>

extern int x_size;
extern int y_size;
extern int z_size;
extern int t_size;

template <class T>
T staples_first(const std::vector<T> &vec, link1 &link, int eta);
template <class T>
T staples_second(const std::vector<std::vector<T>> &smearing_first, link1 &link,
                 std::unordered_map<int, int> &indexes, int rho, int mu,
                 int nu);
template <class T>
T staples_second_refresh(const std::vector<T> &vec, link1 &link, int eta,
                         int nu,
                         double alpha3); // staples for refreshing
                                         // algorythm(refresh link every step)
template <class T>
T staples_third(const std::vector<std::vector<T>> &smearing_second, link1 &link,
                std::unordered_map<int, int> indexes, int nu, int mu);
template <class T>
T staples_third_refresh(const std::vector<T> &vec, link1 &link, int eta,
                        double alpha2, double alpha3);
template <class T>
std::vector<T> smearing_first(const std::vector<T> &array, double alpha3,
                              int mu, int nu, int rho);
template <class T>
std::vector<std::vector<T>> smearing_first_full(const std::vector<T> &array,
                                                double alpha3);
template <class T>
std::vector<T> smearing_second(const std::vector<T> &array,
                               std::vector<std::vector<T>> &smearing_first,
                               double alpha2, int mu, int nu);
template <class T>
std::vector<std::vector<T>>
smearing_second_full(const std::vector<T> &array,
                     std::vector<std::vector<T>> &smearing_first,
                     double alpha2);
template <class T>
std::vector<T> smearing_HYP(const std::vector<T> &array,
                            std::vector<std::vector<T>> &smearing_second,
                            double alpha1);
template <class T>
std::vector<T> smearing_APE(const std::vector<T> &array, double alpha_APE);
template <class T>
std::vector<T> smearing1_APE(const std::vector<T> &array, double alpha_APE);
template <class T>
std::map<std::tuple<int, int>, std::vector<T>>
smearing_APE_2d(const std::vector<T> &array, double alpha_APE);
template <class T>
void smearing_APE_2d_continue(
    std::map<std::tuple<int, int>, std::vector<T>> &smeared, double alpha_APE);
template <class T>
void smearing_APE_2d_continue_plane(
    std::map<std::tuple<int, int>, std::vector<T>> &smeared, int mu,
    double alpha_APE);
template <class T>
T staples_2d_continue(std::vector<T> &array1, std::vector<T> &array2,
                      link1 &link, int i, int j);
template <class T>
T smearing_first_refresh(const std::vector<T> &vec, link1 &link, int nu,
                         int rho,
                         double alpha3); // refresh link every step
template <class T>
T smearing_second_refresh(const std::vector<T> &vec, link1 &link, int nu,
                          double alpha2,
                          double alpha3); // refresh link every step
template <class T>
std::vector<T> smearing_HYP_refresh(data<T> &conf, double alpha1, double alpha2,
                                    double alpha3); // refresh link every step
template <class T>
std::vector<T>
smearing_APE_refresh(data<T> &conf,
                     double alpha_APE); // refresh link every step
std::vector<su2> smearing_stout(data<su2> &conf, double rho);
su2 stout_factor(data<su2> &conf, link1 &link, double rho);
su2 stout_omega(data<su2> &conf, link1 &link, double rho);

void make_map_first(std::unordered_map<int, int> &indexes);
void make_map_second(std::unordered_map<int, int> &indexes);