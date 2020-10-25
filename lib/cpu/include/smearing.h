#pragma once

#include "link.h"
#include <iostream>
#include <unordered_map>
#include <vector>

using namespace std;

extern int x_size;
extern int y_size;
extern int z_size;
extern int t_size;

template <class T> T staples_first(const vector<T> &vec, link1 &link, int eta);
template <class T>
T staples_second(const vector<vector<T>> &smearing_first, link1 &link,
                 unordered_map<int, int> &indexes, int rho, int mu, int nu);
template <class T>
T staples_second_refresh(const vector<T> &vec, link1 &link, int eta, int nu,
                         FLOAT alpha3); // staples for refreshing
                                        // algorythm(refresh link every step)
template <class T>
T staples_third(const vector<vector<T>> &smearing_second, link1 &link,
                unordered_map<int, int> indexes, int nu, int mu);
template <class T>
T staples_third_refresh(const vector<T> &vec, link1 &link, int eta,
                        FLOAT alpha2, FLOAT alpha3);
template <class T>
vector<T> smearing_first(const vector<T> &array, FLOAT alpha3, int mu, int nu,
                         int rho);
template <class T>
vector<vector<T>> smearing_first_full(const vector<T> &array, FLOAT alpha3);
template <class T>
vector<T> smearing_second(const vector<T> &array,
                          vector<vector<T>> &smearing_first, FLOAT alpha2,
                          int mu, int nu);
template <class T>
vector<vector<T>> smearing_second_full(const vector<T> &array,
                                       vector<vector<T>> &smearing_first,
                                       FLOAT alpha2);
template <class T>
vector<T> smearing_HYP(const vector<T> &array,
                       vector<vector<T>> &smearing_second, FLOAT alpha1);
template <class T>
vector<T> smearing_APE(const vector<T> &array, FLOAT alpha_APE);
template <class T>
T smearing_first_refresh(const vector<T> &vec, link1 &link, int nu, int rho,
                         FLOAT alpha3); // refresh link every step
template <class T>
T smearing_second_refresh(const vector<T> &vec, link1 &link, int nu,
                          FLOAT alpha2,
                          FLOAT alpha3); // refresh link every step
template <class T>
vector<T> smearing_HYP_refresh(data<T> &conf, FLOAT alpha1, FLOAT alpha2,
                               FLOAT alpha3); // refresh link every step
template <class T>
vector<T> smearing_APE_refresh(data<T> &conf,
                               FLOAT alpha_APE); // refresh link every step
vector<su2> smearing_stout(data<su2> &conf, FLOAT rho);
su2 stout_factor(data<su2> &conf, link1 &link, FLOAT rho);
su2 stout_omega(data<su2> &conf, link1 &link, FLOAT rho);

void make_map_first(unordered_map<int, int> &indexes);
void make_map_second(unordered_map<int, int> &indexes);