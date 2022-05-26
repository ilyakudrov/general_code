#pragma once

#include <vector>

template <class T>
std::vector<std::vector<T>> separate_smearing_unchanged(std::vector<T> &conf);

template <class T>
void smearing_plane_minor_test1(std::vector<T> &smeared,
                                std::vector<T> &conf_mu,
                                std::vector<T> &conf_nu, int size_mu1,
                                int size_mu2, int size_nu1, int size_nu2,
                                double alpha);

template <class T>
void smearing_plane_major_test1(std::vector<T> &smeared,
                                std::vector<T> &conf_mu,
                                std::vector<T> &conf_nu, int size_mu1,
                                int size_mu2, int size_nu1, int size_nu2,
                                double alpha);

template <class T>
void smearing_APE_test1(std::vector<std::vector<T>> &conf, double alpha);

template <class T>
void smearing_plane_HYP_minor_test1(std::vector<T> &smeared,
                                    std::vector<T> &conf_mu,
                                    std::vector<T> &conf_nu, int size_mu1,
                                    int size_mu2, int size_nu1, int size_nu2,
                                    double alpha, double divisor);

template <class T>
void smearing_plane_HYP_major_test1(std::vector<T> &smeared,
                                    std::vector<T> &conf_mu,
                                    std::vector<T> &conf_nu, int size_mu1,
                                    int size_mu2, int size_nu1, int size_nu2,
                                    double alpha, double divisor);

template <class T>
void smearing_HYP_test1(std::vector<std::vector<T>> &conf, double alpha1,
                        double alpha2, double alpha3);