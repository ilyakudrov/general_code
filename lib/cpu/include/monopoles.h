#pragma once

#ifdef DOUBLE
#define FLOAT double
#else
#define FLOAT float
#endif

#include "../include/link.h"
#include "../include/loop.h"

vector<FLOAT> read_angles_float_fortran(string &file_path);

int find_current(link1 &link, vector<FLOAT> &J);

void find_cluster(link1 &link, vector<FLOAT> &J, loop &ll, int dir_back);
void calculate_clusters(vector<loop *> &LL, vector<FLOAT> &J);

vector<FLOAT> calculate_current(vector<FLOAT> &angles);

// monopole observables

void cluster_length(loop &ll, int &length);

/*template <class T> void length(loop *ll, int &ss1);
template <class T>
result calculate_cluster_lengths(vector<loop *> &LL, int &max_number);
template <class T> void length_mu(loop *ll, int mu, int &s);
template <class T>
void calculate_t_clusters(vector<loop *> &LL, vector<loop *> &t_clusters,
                          int max_number);
template <class T>
void calculate_t_clusters_n(vector<loop *> &LL, vector<loop *>
&t_clusters_n, int max_number, int n); template <class T> void
calculate_s_clusters(vector<loop *> &LL, vector<loop *> &s_clusters, int
max_number); template <class T> FLOAT t_density_n(vector<loop *>
&t_clusters, int n); template <class T> FLOAT
time_length_portion(vector<loop *> &t_clusters); template <class T> void
sites_unique(loop *ll, vector<loop *> &sites); template <class T> FLOAT
*aver_r(vector<loop *> sites); template <class T> FLOAT
distance_shortest(FLOAT a, FLOAT b); template <class T> FLOAT
disp_r(vector<loop *> &sites, FLOAT *aver_coord); template <class T> FLOAT
calculate_disp_r(vector<loop *> &t_clusters); template <class T> bool
sites_close(loop *l, loop *ll); template <class T> FLOAT
dimension(vector<loop *> sites); template <class T> FLOAT
charge_difference(vector<loop *> &t_clusters_1);*/