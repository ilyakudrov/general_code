#pragma once

#include "data.h"
#include "link.h"

#include <map>
#include <unordered_map>
#include <vector>

std::vector<su2> smearing_stout(Data::data<su2> &conf, double rho);
template <class T>
std::vector<std::vector<T>> separate_smearing(std::vector<T> &conf);
template <class T> void smearing_APE(std::vector<T> &conf, double alpha);
template <class T>
void smearing_APE_2d(std::vector<T> &conf1, std::vector<T> &conf2, int mu,
                     int nu, double alpha);
template <class T>
void smearing_APE_2d_test(std::vector<T> &conf1, std::vector<T> &conf2, int mu,
                          int nu, double alpha);
template <class T>
void smearing_HYP(std::vector<T> &array, double alpha1, double alpha2,
                  double alpha3);
