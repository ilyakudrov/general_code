#pragma once

#include "../include/matrix.h"

#include <vector>

std::vector<double> convert_to_angles(std::vector<su2> conf_su2);

std::vector<abelian> convert_to_abelian(std::vector<su2> conf_su2);

std::vector<complex_t> convert_to_complex(std::vector<su2> conf_su2);

std::vector<double> generate_gauge_angles_uniform();

std::vector<abelian> generate_gauge_abelian_uniform();

std::vector<complex_t> generate_gauge_complex_uniform();

std::vector<double> generate_random_numbers(int vector_size);

void heat_bath_test1(abelian &g, abelian &K, double temperature,
                     double &random_numbers);

abelian contribution_site_test1(std::vector<abelian> &gauge_abelian,
                                std::vector<abelian> &conf_abelian, int x,
                                int y, int z, int t, int position,
                                std::vector<int> &shift);

void heat_bath_update_test1(std::vector<abelian> &gauge_abelian,
                            std::vector<abelian> &conf_abelian,
                            double temperature);

void normalize_abelian(std::vector<abelian> &abelian);

double Landau_functional_gauge_abelian(std::vector<abelian> &gauge_abelian,
                                       std::vector<abelian> &conf_abelian);

double Landau_functional_gauge_angles(std::vector<double> &gauge_angles,
                                      std::vector<abelian> &conf_abelian);

double Landau_functional_gauge(std::vector<double> &gauge_angles,
                               std::vector<double> &conf_angles);

double Landau_functional_gauge_complex(std::vector<complex_t> &gauge_complex,
                                       std::vector<complex_t> &conf_complex);

double Landau_functional(std::vector<su2> &conf);

double Landau_functional_abelian(std::vector<abelian> &conf_abelian);

void gauge_tranformation_abelian(std::vector<abelian> &gauge_abelian,
                                 std::vector<abelian> &conf_abelian);

void gauge_tranformation(std::vector<abelian> &gauge_abelian,
                         std::vector<su2> &conf_su2);

void overrelaxation_update_test1(std::vector<abelian> &gauge_abelian,
                                 std::vector<abelian> &conf_abelian);

std::tuple<double, double>
relaxation_update_test1(std::vector<abelian> &gauge_abelian,
                        std::vector<abelian> &conf_abelian);

void make_simulated_annealing_test1(std::vector<abelian> &gauge_abelian,
                                    std::vector<abelian> &conf_abelian,
                                    double T_init, double T_final,
                                    double T_step, int OR_steps,
                                    int thermalization_steps);

void make_maximization_approximate(std::vector<abelian> &gauge_abelian,
                                   std::vector<abelian> &conf_abelian,
                                   int OR_steps, int tolerance_digits);

void make_maximization_final(std::vector<abelian> &gauge_abelian,
                             std::vector<abelian> &conf_abelian, int OR_steps,
                             double tolerance_maximal,
                             double tolerance_average);

void heat_bath_test2(double &g, abelian &K, double temperature,
                     double &random_numbers);

abelian contribution_site_test2(std::vector<double> &gauge_abelian,
                                std::vector<abelian> &conf_abelian, int x,
                                int y, int z, int t, int position,
                                std::vector<int> &shift);

void heat_bath_update_test2(std::vector<double> &gauge_angles,
                            std::vector<abelian> &conf_abelian,
                            double temperature);

abelian contribution_site_test3(std::vector<double> &gauge_angles,
                                std::vector<double> &conf_angles, int x, int y,
                                int z, int t, int position,
                                std::vector<int> &shift);

void heat_bath_update_test3(std::vector<double> &gauge_angles,
                            std::vector<double> &conf_angles,
                            double temperature);

void heat_bath(complex_t &g, complex_t &K, double temperature,
               double &random_numbers);

complex_t contribution_site(std::vector<complex_t> &gauge_complex,
                            std::vector<complex_t> &conf_complex, int x, int y,
                            int z, int t, int position,
                            std::vector<int> &shift);

void heat_bath_update(std::vector<complex_t> &gauge_complex,
                      std::vector<complex_t> &conf_complex, double temperature);

void overrelaxation_update(std::vector<complex_t> &gauge_complex,
                           std::vector<complex_t> &conf_complex);

std::tuple<double, double>
relaxation_update(std::vector<complex_t> &gauge_complex,
                  std::vector<complex_t> &conf_complex);

void make_simulated_annealing(std::vector<complex_t> &gauge_complex,
                              std::vector<complex_t> &conf_complex,
                              double T_init, double T_final, double T_step,
                              int OR_steps, int thermalization_steps);

void make_maximization_final(std::vector<complex_t> &gauge_complex,
                             std::vector<complex_t> &conf_complex, int OR_steps,
                             double tolerance_maximal,
                             double tolerance_average);
