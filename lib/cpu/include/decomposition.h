#pragma once

#include "../include/link.h"

#include <vector>

std::vector<double> read_double_angles(std::string &file_name, int bites_skip);

void write_double_angles(std::string &file_name, std::vector<double> &angles);

std::vector<double> merge_angles(std::vector<std::vector<double>> &angles);

std::vector<su2> get_monopoless(std::vector<su2> &conf_su2,
                                std::vector<double> &angles_monopole);

std::vector<su2> get_initial_su2(std::vector<su2> &conf_monopoless,
                                 std::vector<double> &angles_monopole);

std::vector<double> read_gauge_Landau(std::string &file_name, int bites_skip);

void apply_gauge_Landau(std::vector<su2> &conf_su2, std::vector<double> &gauge);

std::vector<double> read_inverse_laplacian(std::string &file_path);

double get_monopole_angle(std::vector<std::vector<int>> &monopole_plaket,
                          link1 &link_tmp, std::vector<double> &laplace,
                          int mu);

std::vector<double> make_monopole_angles2(std::vector<double> &angles,
                                          std::vector<double> &laplace);

void monopole_plaket_difference_nonzero(
    std::vector<std::vector<int>> &monopole_plaket,
    std::vector<std::vector<int>> &monopole_difference,
    std::vector<std::vector<int>> &monopole_coordinate);

std::vector<double> make_monopole_angles1(std::vector<double> &angles,
                                          std::vector<double> &laplace);

void decomposition_step2(std::vector<std::vector<int>> &monopole_difference,
                         std::vector<std::vector<int>> &monopole_coordinate,
                         std::vector<double> &laplace,
                         std::vector<double> &angles_decomposed, int i, int mu);

std::vector<double> make_monopole_angles2(std::vector<double> &angles,
                                          std::vector<double> &laplace);

void decomposition_step3(std::vector<std::vector<int>> &monopole_difference,
                         std::vector<std::vector<int>> &monopole_coordinate,
                         std::vector<double> &laplace,
                         std::vector<std::vector<double>> &angles_decomposed,
                         int i, int mu);

std::vector<std::vector<double>>
make_monopole_angles3(std::vector<double> &angles,
                      std::vector<double> &laplace);

void decomposition_step(std::vector<std::vector<int>> &monopole_difference,
                        std::vector<std::vector<int>> &monopole_coordinate,
                        std::vector<double> &laplace,
                        std::vector<std::vector<double>> &angles_decomposed,
                        int i, int mu);

std::vector<double> make_monopole_angles(std::vector<double> &angles,
                                         std::vector<double> &laplace);