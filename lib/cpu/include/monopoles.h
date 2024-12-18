#pragma once

#include "../include/link.h"
#include "../include/loop.h"

#include <unordered_map>

std::vector<double> read_angles_double_fortran(std::string &file_path);
std::vector<double> read_angles_double_fortran(std::string &file_path);
std::vector<double> read_double_fortran_convet_abelian(std::string &file_path);
std::vector<double> read_double_fortran_convet_abelian(std::string &file_path);
std::vector<std::vector<double>>
read_double_su3_convet_angles(std::string &file_path);
std::vector<double> read_double_qc2dstag_convet_abelian(std::string &file_path);
std::vector<double> convert_to_angles(std::vector<abelian> &conf_abelian);
std::vector<std::vector<double>>
convert_to_angles(std::vector<su3_abelian> &conf);
std::vector<std::vector<double>> convert_to_angles(std::vector<su3> &conf);

std::vector<std::vector<double>>
calculate_monopole_plaket(std::vector<double> &angles);

std::vector<std::vector<int>>
calculate_monopole_plaket_singular(std::vector<double> &angles);

template <class T> int find_current(link1 &link, std::vector<T> &J);

template <class T>
std::vector<loop *> find_paths(std::vector<loop *> &neighbours,
                               std::vector<T> &J);

template <class T> void find_cluster(loop *ll, std::vector<T> &J);

template <class T> std::vector<loop *> calculate_clusters(std::vector<T> &J);

std::vector<double> calculate_current(std::vector<double> &angles);

std::vector<int> calculate_current_singular(std::vector<double> &angles);

std::vector<double> calculate_current_monopole_plakets(
    std::vector<std::vector<double>> monopole_plakets);

std::vector<std::vector<std::vector<double>>>
make_monopole_plakets(std::vector<std::vector<double>> &angles);

std::vector<std::vector<std::vector<int>>>
make_monopole_plakets_singular(std::vector<std::vector<double>> &angles);

void make_plakets_both(
    std::vector<std::vector<double>> &angles,
    std::vector<std::vector<std::vector<double>>> &monopole_plakets,
    std::vector<std::vector<std::vector<int>>> &dirac_plakets);

// functions for obtaining information about clusters for testing

void print_currents(loop *ll);

void check_for_coordinate(loop *loop, int coordinate[4], bool &include);

// monopole observables

int cluster_length(loop *ll);

void cluster_length_recurrent(loop *ll, int &length);

void cluster_sites(loop *ll);

std::vector<int> length_mu(loop *ll);

void length_mu_recurrent(loop *ll, std::vector<int> &lengths_mu);

std::vector<int> currents_directions(loop *ll);

void currents_directions_recurrent(loop *ll, std::vector<int> &directions);

double cluster_variation(loop *loop);

void cluster_variation_recurrent(loop *loop, double &variation,
                                 std::vector<int> distance, const link1 &link);

int site_number(loop *loop);

void site_number_recurrent(loop *loop, int &link_number,
                           std::unordered_map<int, int> &loop_sites);

bool increment_position(std::vector<int> &positions, int max_position);

void iterate_sum_order(std::vector<std::vector<int>> &wrappings, int sum_order);

std::vector<int> group_percolating(std::vector<std::vector<int>> &wrappings);