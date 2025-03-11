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
    std::vector<std::complex<double>> &furier_coefficients);
