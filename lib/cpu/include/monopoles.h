#pragma once

#ifdef DOUBLE
#define FLOAT double
#else
#define FLOAT float
#endif

#include "../include/link.h"
#include "../include/loop.h"

vector<FLOAT> read_angles_float_fortran(string &file_path);
vector<FLOAT> read_float_fortran_convet_abelian(string &file_path);
vector<FLOAT> read_double_qc2dstag_convet_abelian(string &file_path);

int find_current(link1 &link, vector<FLOAT> &J);

vector<loop *> find_paths(vector<loop *> &neighbours, vector<FLOAT> &J);

void find_cluster(loop *ll, vector<FLOAT> &J);

vector<loop *> calculate_clusters(vector<FLOAT> &J);

vector<FLOAT> calculate_current(vector<FLOAT> &angles);

// monopole observables

void cluster_length(loop *ll, int &length);

void cluster_sites(loop *ll);

void length_mu(loop *ll, vector<int> &lengths_mu);