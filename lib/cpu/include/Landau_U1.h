#pragma once

#include "../include/data.h"
#include "../include/indexing.h"
#include "../include/matrix.h"

#include <map>
#include <vector>

std::vector<double> convert_to_angles(const std::vector<su2> &conf_su2);
std::vector<double> convert_to_angles(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf_su2);

std::vector<abelian> convert_to_abelian(const std::vector<su2> &conf_su2);

std::vector<std::complex<double>>
convert_to_complex(const std::vector<su2> &conf_su2);

std::vector<std::complex<double>> convert_to_complex(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf_su2);
std::vector<std::complex<double>> convert_to_complex(
    const Data::LatticeData<DataPatternLexicographical, abelian> &conf);

std::vector<double> convert_complex_to_angles(
    const std::vector<std::complex<double>> &conf_complex);
std::vector<abelian> convert_complex_to_abelian(
    const std::vector<std::complex<double>> &conf_complex,
    DataPatternLexicographical &data_pattern);

std::vector<double> generate_gauge_angles_uniform();

std::vector<abelian> generate_gauge_abelian_uniform();

std::vector<std::complex<double>> generate_gauge_complex_uniform();
std::vector<std::complex<double>>
generate_gauge_complex_uniform(DataPatternLexicographical &data_pattern);

std::vector<std::complex<double>> generate_gauge_complex_unity();

std::vector<double>
generate_gauge_angles_uniform(DataPatternLexicographical &data_pattern);

std::vector<abelian>
generate_gauge_abelian_uniform(DataPatternLexicographical &data_pattern);

std::vector<std::complex<double>>
generate_gauge_complex_uniform(DataPatternLexicographical &data_pattern);

std::vector<std::complex<double>>
generate_gauge_complex_unity(DataPatternLexicographical &data_pattern);

std::vector<double> generate_random_numbers(int vector_size);

void heat_bath_test1(abelian &g, const abelian &K, double temperature,
                     double &random_numbers);

abelian contribution_site_test1(const std::vector<abelian> &gauge_abelian,
                                const std::vector<abelian> &conf_abelian, int x,
                                int y, int z, int t, int position,
                                std::vector<int> &shift);

void heat_bath_update_test1(std::vector<abelian> &gauge_abelian,
                            const std::vector<abelian> &conf_abelian,
                            double temperature);

void normalize_abelian(std::vector<abelian> &abelian);

void normalize_complex(std::vector<std::complex<double>> &gauge_complex);

double
Landau_functional_gauge_abelian(const std::vector<abelian> &gauge_abelian,
                                const std::vector<abelian> &conf_abelian);

double Landau_functional_gauge_angles(const std::vector<double> &gauge_angles,
                                      const std::vector<abelian> &conf_abelian);

double Landau_functional_gauge(const std::vector<double> &gauge_angles,
                               const std::vector<double> &conf_angles);

double Landau_functional_conf_complex(
    const std::vector<std::complex<double>> &conf_complex,
    const std::vector<std::complex<double>> &gauge_complex,
    DataPatternLexicographical &data_pattern);

double Landau_functional(const std::vector<su2> &conf);

double Landau_functional_abelian(std::vector<abelian> &conf_abelian);

double Landau_functional_complex(
    const std::vector<std::complex<double>> &conf_complex);

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

void heat_bath(std::complex<double> &g, std::complex<double> &K,
               double temperature, double *random_numbers);

std::complex<double>
contribution_site(std::vector<std::complex<double>> &gauge_complex,
                  std::vector<std::complex<double>> &conf_complex, int x, int y,
                  int z, int t, int position, std::vector<int> &shift);

void heat_bath_update(std::vector<std::complex<double>> &gauge_complex,
                      std::vector<std::complex<double>> &conf_complex,
                      double temperature);

void overrelaxation_update(std::vector<std::complex<double>> &gauge_complex,
                           std::vector<std::complex<double>> &conf_complex);

std::tuple<double, double>
relaxation_update(std::vector<std::complex<double>> &gauge_complex,
                  std::vector<std::complex<double>> &conf_complex);

void make_simulated_annealing(std::vector<std::complex<double>> &gauge_complex,
                              std::vector<std::complex<double>> &conf_complex,
                              double T_init, double T_final, double T_step,
                              int OR_steps, int thermalization_steps);

void make_maximization_final(std::vector<std::complex<double>> &gauge_complex,
                             std::vector<std::complex<double>> &conf_complex,
                             int OR_steps, double tolerance_maximal,
                             double tolerance_average);

void apply_gauge_Landau_complex(
    std::vector<std::complex<double>> &gauge_complex,
    std::vector<std::complex<double>> &conf_complex);

void apply_gauge_Landau(std::vector<std::complex<double>> &gauge_complex,
                        std::vector<su2> &conf_su2);

void apply_gauge_Landau_complex(
    std::vector<std::complex<double>> &gauge_complex,
    std::vector<std::complex<double>> &conf_complex,
    DataPatternLexicographical &data_pattern);

void apply_gauge_Landau(
    std::vector<std::complex<double>> &gauge_complex,
    Data::LatticeData<DataPatternLexicographical, su2> &conf_su2);

std::map<double, double> simulated_annealing_thermalization_test(
    const std::vector<std::complex<double>> &conf_complex,
    std::vector<std::complex<double>> &gauge_complex,
    DataPatternLexicographical &data_pattern, double T_init, double T_final,
    double T_step, int OR_steps, int thermalization_steps,
    int local_thermalization_steps);

void make_simulated_annealing(
    const std::vector<std::complex<double>> &conf_complex,
    std::vector<std::complex<double>> &gauge_complex,
    DataPatternLexicographical &data_pattern, double T_init, double T_final,
    double T_step, int OR_steps, int thermalization_steps);

void make_maximization_final(
    const std::vector<std::complex<double>> &conf_complex,
    std::vector<std::complex<double>> &gauge_complex,
    DataPatternLexicographical &data_pattern, int OR_steps,
    double tolerance_maximal, double tolerance_average);
