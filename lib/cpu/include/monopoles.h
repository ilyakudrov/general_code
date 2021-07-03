#pragma once

#ifdef DOUBLE
#define FLOAT double
#else
#define FLOAT float
#endif

#include "../include/link.h"
#include "../include/loop.h"

std::vector<FLOAT> read_angles_float_fortran(std::string &file_path);
std::vector<FLOAT> read_angles_double_fortran(std::string &file_path);
std::vector<FLOAT> read_float_fortran_convet_abelian(std::string &file_path);
std::vector<FLOAT> read_double_fortran_convet_abelian(std::string &file_path);
std::vector<FLOAT> read_double_qc2dstag_convet_abelian(std::string &file_path);

int find_current(link1 &link, std::vector<FLOAT> &J);

std::vector<loop *> find_paths(std::vector<loop *> &neighbours,
                               std::vector<FLOAT> &J);

void find_cluster(loop *ll, std::vector<FLOAT> &J);

std::vector<loop *> calculate_clusters(std::vector<FLOAT> &J);

std::vector<FLOAT> calculate_current(std::vector<FLOAT> &angles);

// monopole observables

void cluster_length(loop *ll, int &length);

void cluster_sites(loop *ll);

void length_mu(loop *ll, std::vector<int> &lengths_mu);