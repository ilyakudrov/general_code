#pragma once

#include "../include/link.h"

#include <vector>

std::vector<double> read_double_angles(std::string &file_name, int bites_skip);

void write_double_angles(std::string &file_name,
                         const std::vector<double> &angles);

std::vector<std::vector<double>> read_double_angles_su3(std::string &file_name);

void write_double_angles_su3(std::string &file_name,
                             const std::vector<std::vector<double>> &angles);

void write_double_su2(std::string &file_name, const std::vector<su2> &conf_su2);

void write_double_su3(std::string &file_name, const std::vector<su3> &conf_su3);

std::vector<std::vector<double>>
get_angles_su3(const std::vector<su3> &conf_su3);

std::vector<su3_abelian> get_abelian(std::vector<su3> &conf);

std::vector<abelian> get_abelian(std::vector<su2> &conf);

std::vector<su3> get_offdiagonal(std::vector<su3> &conf);

std::vector<su2> get_offdiagonal(std::vector<su2> &conf);

std::vector<double>
merge_angles(const std::vector<std::vector<double>> &angles);

std::vector<su2> get_monopoless(const std::vector<su2> &conf_su2,
                                const std::vector<double> &angles_monopole);

void get_monopoless_optimized(std::vector<su2> &conf_su2,
                              const std::vector<double> &angles_monopole);

void get_monopoless_optimized_su3(
    std::vector<su3> &conf_su3,
    const std::vector<std::vector<double>> &angles_monopole);

std::vector<su2> get_initial_su2(std::vector<su2> &conf_monopoless,
                                 std::vector<double> &angles_monopole);

std::vector<double> read_gauge_Landau(std::string &file_name, int bites_skip);

void apply_gauge_Landau(std::vector<su2> &conf_su2,
                        const std::vector<double> &gauge);

std::vector<double> read_inverse_laplacian(std::string &file_path);

double get_monopole_angle(std::vector<std::vector<int>> &monopole_plaket,
                          link1 &link_tmp, std::vector<double> &laplace,
                          int mu);

std::vector<double> make_monopole_angles2(std::vector<double> &angles,
                                          std::vector<double> &laplace);

void dirac_plaket_difference_nonzero(
    const std::vector<std::vector<int>> &dirac_plaket,
    std::vector<std::vector<int>> &dirac_difference,
    std::vector<std::vector<int>> &dirac_coordinate);

std::vector<double>
make_monopole_angles1(const std::vector<std::vector<int>> &dirac_plakets,
                      const std::vector<double> &laplace);

void decomposition_step2(const std::vector<std::vector<int>> &dirac_difference,
                         const std::vector<std::vector<int>> &dirac_coordinate,
                         const std::vector<double> &laplace,
                         std::vector<double> &angles_decomposed, int i, int mu);

void decomposition_step3(const std::vector<std::vector<int>> &dirac_difference,
                         const std::vector<std::vector<int>> &dirac_coordinate,
                         const std::vector<double> &laplace,
                         std::vector<std::vector<double>> &angles_decomposed,
                         int i, int mu);

std::vector<std::vector<double>>
make_monopole_angles3(const std::vector<std::vector<int>> &dirac_plakets,
                      const std::vector<double> &laplace);

void decomposition_step(const std::vector<std::vector<int>> &dirac_difference,
                        const std::vector<std::vector<int>> &dirac_coordinate,
                        const std::vector<double> &laplace,
                        std::vector<std::vector<double>> &angles_decomposed,
                        int i, int mu);

std::vector<double>
make_monopole_angles(std::vector<std::vector<int>> &dirac_plakets,
                     const std::vector<double> &laplace);

void decomposition_step_parallel(int dirac_difference,
                                 const std::vector<int> &dirac_coordinate,
                                 const std::vector<double> &laplace,
                                 std::vector<double> &angles_decomposed);

void decomposition_step_parallel3_simple_positive(
    const std::vector<int> &dirac_coordinate,
    const std::vector<double> &laplace, std::vector<double> &angles_decomposed);

void decomposition_step_parallel3_simple_negative(
    const std::vector<int> &dirac_coordinate,
    const std::vector<double> &laplace, std::vector<double> &angles_decomposed);

void decomposition_step_parallel1(
    const std::vector<std::vector<int>> &dirac_difference,
    const std::vector<std::vector<int>> &dirac_coordinate,
    const std::vector<double> &laplace,
    std::vector<std::vector<double>> &angles_decomposed, int i, int mu);

void decomposition_step_parallel2(
    const std::vector<std::vector<int>> &dirac_difference,
    const std::vector<std::vector<int>> &dirac_coordinate,
    const std::vector<double> &laplace,
    std::vector<std::vector<double>> &angles_decomposed, int i, int mu);

std::vector<double>
make_monopole_angles_parallel(std::vector<std::vector<int>> &dirac_plakets,
                              const std::vector<double> &laplace);

void calculate_inverse_laplacian(std::vector<double> &inverse_laplacian_real,
                                 std::vector<double> &inverse_laplacian_imag);