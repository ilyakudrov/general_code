#pragma once

#include "link.h"
#include "matrix.h"

#include <random>
#include <vector>

struct SA_data {
  double temperature;
  double functional;
};

double functional_Landau(std::vector<su2> &conf);

double functional_Landau(std::vector<abelian> &conf);

double functional_Landau(std::vector<abelian> &angles_conf,
                         std::vector<double> &angles_gauge);

std::vector<abelian> convert_su2_abelian(std::vector<su2> &conf_su2);

std::vector<double> generate_angles_gauge_random();

std::vector<double> generate_angles_gauge_trivial();

abelian angle_neighbour(std::vector<double> &angles_gauge,
                        std::vector<abelian> &angles_conf, link1 link);

void heat_bath_local(
    double *angles_gauge, abelian neighbour, double temperature,
    std::subtract_with_carry_engine<unsigned, 24, 10, 24> &random_generator);

void heat_bath_step(std::vector<double> &angles_gauge,
                    std::vector<abelian> &angles_conf, double temperature);

void relaxation_step(std::vector<double> &angles_gauge,
                     std::vector<abelian> &angles_conf);

void overrelaxation_step(std::vector<double> &angles_gauge,
                         std::vector<abelian> &angles_conf);

std::vector<double> calculate_binomial_coefficients(int w, int order);

void microcanonical_step(std::vector<double> &angles_gauge,
                         std::vector<abelian> &angles_conf, double w,
                         int order);

void thermalize(std::vector<double> &angles_gauge,
                std::vector<abelian> &angles_conf, double temperature,
                int therm_steps);

void simulated_annealing(std::vector<double> &angles_gauge,
                         std::vector<abelian> &angles_conf,
                         double temperature_start, double temperature_end,
                         double temperature_step);

std::vector<SA_data> simulated_annealing_test(std::vector<double> &angles_gauge,
                                              std::vector<abelian> &angles_conf,
                                              double temperature_start,
                                              double temperature_end,
                                              double temperature_step);
