#pragma once

#include "../include/matrix.h"

#include <array>
#include <complex>
#include <vector>

std::vector<std::array<double, 12>>
get_vector_potential(std::vector<su2> &conf);

std::vector<std::complex<double>>
get_furier_coefficients(std::array<double, 4> &momenta);

std::array<std::complex<double>, 144> calculate_gluon_propagator(
    std::vector<std::array<double, 12>> &vector_potential,
    std::vector<std::complex<double>> &furier_coefficients, double multiplier);

std::array<std::complex<double>, 12> calculate_gluon_propagator_diagonal(
    std::vector<std::array<double, 12>> &vector_potential,
    std::vector<std::complex<double>> &furier_coefficients, double multiplier);

std::vector<std::array<double, 3>>
get_vector_potential_longitudinal(std::vector<su2> &conf);

double calculate_gluon_propagator_longitudinal_zero_momentum(
    std::vector<std::array<double, 3>> &vector_potential, double multiplier);

std::complex<double> calculate_gluon_propagator_single(
    std::vector<std::array<double, 12>> &vector_potential,
    std::vector<std::complex<double>> &furier_coefficients, double multiplier,
    int m, int n);

std::vector<std::vector<std::array<double, 4>>> generate_momenta(int Ns,
                                                                 int Nt);