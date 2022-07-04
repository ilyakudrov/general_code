#pragma once

#include "../include/link.h"

#include <vector>

std::vector<double> read_double_angles(std::string &file_name, int bites_skip);

std::vector<double> read_inverse_laplacian(std::string &file_path);

double get_monopole_angle(std::vector<std::vector<int>> &monopole_plaket,
                          link1 &link_tmp, std::vector<double> &laplace,
                          int mu);

std::vector<double> make_monopole_angles(std::vector<double> &angles,
                                         std::vector<double> &laplace);

void monopole_plaket_difference_nonzero(
    std::vector<std::vector<int>> &monopole_plaket,
    std::vector<std::vector<int>> &monopole_difference,
    std::vector<std::vector<int>> &monopole_coordinate);

std::vector<double> make_monopole_angles1(std::vector<double> &angles,
                                          std::vector<double> &laplace);