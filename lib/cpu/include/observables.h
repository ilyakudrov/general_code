#pragma once

#ifdef DOUBLE
#define FLOAT double
#else
#define FLOAT float
#endif

#include "data.h"
#include "link.h"
#include "matrix.h"
#include "result.h"

using namespace std;

// Plaket
template <class T> FLOAT plaket_time(const vector<T> &array);
template <class T> FLOAT plaket_space(const vector<T> &array);
template <class T> FLOAT plaket(const vector<T> &array);

// Wilson loop
template <class T>
vector<FLOAT> wilson(const vector<T> &array, int r_min, int r_max, int time_min,
                     int time_max);
template <class T>
vector<T> wilson_lines(const vector<T> &array, int mu, int length);
template <class T>
vector<T> wilson_line_increase(const vector<T> &array, const vector<T> &lines,
                               int mu, int length);

// Polyakov loop
template <class T> FLOAT polyakov(const vector<T> &array);

// Wilson_plaket_schwinger_correlator
template <class T> vector<T> calculate_plaket_time(const vector<T> array);
template <class T> vector<T> calculate_plaket_space(const vector<T> &array);
template <class T> vector<T> calculate_polyakov_loop(const vector<T> &array);
template <class T>
vector<T> calculate_wilson_loop(const vector<T> &array, int r, int time);
template <class T>
vector<vector<T>> calculate_schwinger_line(const vector<T> &array, int d,
                                           int x_trans);
template <class T>
void wilson_plaket_schwinger_electric(const vector<vector<T>> &schwinger_line,
                                      const vector<T> &plaket,
                                      const vector<T> &polyakov_loop,
                                      vector<vector<result>> &field1,
                                      vector<vector<result>> &field2,
                                      vector<result> &field3, int d, int D,
                                      int x_trans);

// Wilson_plaket_correlator
template <class T>
vector<FLOAT> calculate_plaket_time_tr(const vector<T> &array);
template <class T>
vector<FLOAT> calculate_plaket_space_tr(const vector<T> &array);
template <class T>
FLOAT plaket4_time_optimized(const vector<FLOAT> &plaket_tr, link1<T> &link);
template <class T>
FLOAT plaket4_space_optimized(const vector<FLOAT> &plaket_tr, link1<T> &link,
                              int nu);
template <class T>
vector<FLOAT> calculate_wilson_loop_tr(const vector<T> &array, int r, int time);

result wilson_plaket_correlator_electric(const vector<FLOAT> &wilson_loop_tr,
                                         const vector<FLOAT> &plaket_tr, int r,
                                         int time, int x_trans, int d_min,
                                         int d_max);
result wilson_plaket_correlator_electric_x(const vector<FLOAT> &wilson_loop_tr,
                                           const vector<FLOAT> &plaket_tr,
                                           int r, int time, int x_trans_min,
                                           int x_trans_max, int d);
result wilson_plaket_correlator_magnetic(const vector<FLOAT> &wilson_loop_tr,
                                         const vector<FLOAT> &plaket_tr, int r,
                                         int time, int x_trans, int d_min,
                                         int d_max);
result wilson_plaket_correlator_magnetic_x(const vector<FLOAT> &wilson_loop_tr,
                                           const vector<FLOAT> &plaket_tr,
                                           int r, int time, int x_trans_min,
                                           int x_trans_max, int d);

// Polyakov_correlator
template <class T> FLOAT polyakov_loop_corelator(const vector<T> &array, int D);

FLOAT MAG_functional_su2(const vector<su2> &array);

// monopoles
/*template <class T> void length(loop *ll, int &ss1);
template <class T>
result calculate_cluster_lengths(vector<loop *> &LL, int &max_number);
template <class T> void length_mu(loop *ll, int mu, int &s);
template <class T>
void calculate_t_clusters(vector<loop *> &LL, vector<loop *> &t_clusters,
                          int max_number);
template <class T>
void calculate_t_clusters_n(vector<loop *> &LL, vector<loop *> &t_clusters_n,
                            int max_number, int n);
template <class T>
void calculate_s_clusters(vector<loop *> &LL, vector<loop *> &s_clusters,
                          int max_number);
template <class T> FLOAT t_density_n(vector<loop *> &t_clusters, int n);
template <class T> FLOAT time_length_portion(vector<loop *> &t_clusters);
template <class T> void sites_unique(loop *ll, vector<loop *> &sites);
template <class T> FLOAT *aver_r(vector<loop *> sites);
template <class T> FLOAT distance_shortest(FLOAT a, FLOAT b);
template <class T> FLOAT disp_r(vector<loop *> &sites, FLOAT *aver_coord);
template <class T> FLOAT calculate_disp_r(vector<loop *> &t_clusters);
template <class T> bool sites_close(loop *l, loop *ll);
template <class T> FLOAT dimension(vector<loop *> sites);
template <class T> FLOAT charge_difference(vector<loop *> &t_clusters_1);*/