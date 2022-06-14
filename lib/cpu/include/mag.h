#pragma once

#include "../include/matrix.h"

#include <vector>

std::vector<double> generate_random_numbers_sphere(int vector_size);

std::vector<spin> generate_spins_uniform();

std::vector<double> generate_random_numbers(int vector_size);

void heat_bath(spin &spins, spin &neighbour, double temperature,
               double *random_numbers);

spin contribution_site(std::vector<spin> &spins, std::vector<su2> &conf_su2,
                       int x, int y, int z, int t, int position,
                       std::vector<int> &shift);

void heat_bath_update_tets1(std::vector<spin> &spins,
                            std::vector<su2> &conf_su2, double temperature);

void normalize_spin(std::vector<spin> &spins);

double MAG_functional_su2_spin(std::vector<su2> &conf_su2,
                               std::vector<spin> &spins);

std::vector<su2> make_gauge(std::vector<spin> &spins);

std::vector<su2> gauge_tranformation(std::vector<su2> &conf_su2,
                                     std::vector<su2> &gauge);

std::vector<int> make_indices_qube(int qube_size);

spin contribution_site2(std::vector<spin> &spins, std::vector<su2> &conf_su2,
                        std::vector<int> &indices, int index);

void heat_bath_update_tets2(std::vector<spin> &spins,
                            std::vector<su2> &conf_su2,
                            std::vector<int> indices, double temperature);

void heat_bath_update_tets3(std::vector<spin> &spins,
                            std::vector<su2> &conf_su2, double temperature);

void make_indices_qube1(std::vector<int> &indices,
                        std::vector<char> &coordinates, int qube_size);

void heat_bath_update_tets4(std::vector<spin> &spins,
                            std::vector<su2> &conf_su2,
                            std::vector<int> &indices,
                            std::vector<char> &coordinates, double temperature);