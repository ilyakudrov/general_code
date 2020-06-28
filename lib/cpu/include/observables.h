#ifndef __OBSERVABLES_H__
#define __OBSERVABLES_H__

#include "data.h"
#include "matrix.h"
#include "link.h"
#include "result.h"
#include "monopoles.h"

using namespace std;

template<typename T>
double plaket_time(const vector<T>& array);
template<typename T>
double plaket_space(const vector<T>& array);
template<typename T>
double plaket(const vector<T>& array);
template<typename T>
double wilson(const vector<T>& array, int r, int time);
// double wilson(const vector<matrix>& array, int r, int t);
// double wilson(const vector<double>& array, int r, int t);
template<typename T>
double polyakov(const vector<T>& array);
double wilson_abelian(const vector<matrix>& array, int R, int T);
void polyakov_abelian(const vector<matrix>& array, double pol[2]);
double wilson_dir(const vector<matrix>& array, int R, int T, int dir);
void fields(const vector<vector<matrix> >& schwinger_line, const vector<matrix>& plaket, const vector<matrix>& polyakov_loop, vector<vector<result> >& field1, vector<vector<result> >& field2, vector<result>& field3, int d, int D, int x_trans);
void field1_average(const vector<vector<matrix> >& schwinger_line, const vector<matrix>& plaket, const vector<matrix>& polyakov_loop, vector<vector<result> >& field1, int d, int D, int x_trans) ;
vector<vector<matrix> > calculate_schwinger_line(vector<matrix>& array, int d, int x_trans);
vector<matrix> calculate_plaket(vector<matrix>& array);
vector<double> calculate_plaket_time_tr(const vector<matrix>& array);
vector<double> calculate_plaket_time_tr(const vector<double>& array);
vector<double> calculate_plaket_space_tr(const vector<matrix>& array);
vector<double> calculate_plaket_space_tr(const vector<double>& array);
double plaket4_time_optimized(const vector<double>& plaket_tr, link1<matrix>& link);
double plaket4_space_optimized(const vector<double>& plaket_tr, link1<matrix>& link, int nu);
vector<matrix> calculate_polyakov_loop(const vector<matrix>& array);
vector<matrix> calculate_wilson_loop(const vector<matrix>& array, int R, int T);
vector<double> calculate_wilson_loop_tr(const vector<matrix>& array, int R, int T);
vector<double> calculate_wilson_loop_tr(const vector<double>& array, int R, int T);
double polyakov_loop_corelator(const vector<matrix>& array, int D);
double plaket_correlator(const vector<matrix>& plaket, int dist);
double plaket_correlator_space(const vector<matrix>& plaket, int dist);
double wilson_plaket_correlator_electric_simple(const vector<matrix>& array, const vector<double>& wilson_loop_tr, const vector<double>& plaket_tr, int R, int T, int d, int dir);
result wilson_plaket_correlator_electric(const vector<double>& wilson_loop_tr, const vector<double>& plaket_tr, int R, int T, int x_trans, int d_min, int d_max);
result wilson_plaket_correlator_electric_x(const vector<matrix>& array, const vector<double>& wilson_loop_tr, const vector<double>& plaket_tr, int R, int T, int x_trans_min, int x_trans_max, int d);
result polyakov_plaket_correlator_electric(const vector<matrix>& array, const vector<matrix>& array_smeared, int R, int x_trans, int d_min, int d_max);
result wilson_plaket_correlator_magnetic(const vector<double>& wilson_loop_tr, const vector<double>& plaket_tr, int R, int T, int x_trans, int d_min, int d_max);
result wilson_plaket_correlator_magnetic_x(vector<matrix>& array, const vector<double>& wilson_loop_tr, const vector<double>& plaket_tr, int R, int T, int x_trans_min, int x_trans_max, int d);
result polyakov_plaket_correlator_magnetic(vector<matrix>& array, const vector<matrix>& array_smeared, int R, int x_trans, int d_min, int d_max);
//monopoles
void length(loop* ll, int& ss1);
result calculate_cluster_lengths(vector<loop*>& LL, int& max_number);
void length_mu(loop* ll, int mu, int& s);
void calculate_t_clusters(vector<loop*>& LL, vector<loop*>& t_clusters, int max_number);
void calculate_t_clusters_n(vector<loop*>& LL, vector<loop*>& t_clusters_n, int max_number, int n);
void calculate_s_clusters(vector<loop*>& LL, vector<loop*>&	s_clusters, int max_number);
double t_density_n(vector<loop*>& t_clusters, int n);
double time_length_portion(vector<loop*>& t_clusters);
void sites_unique(loop* ll, vector<loop*>& sites);
double* aver_r(vector<loop*> sites);
double distance_shortest(double a, double b);
double disp_r(vector<loop*>& sites, double* aver_coord);
double calculate_disp_r(vector<loop*>& t_clusters);
bool sites_close(loop* l, loop* ll);
double dimension(vector<loop*> sites);
double charge_difference(vector<loop*>& t_clusters_1);
#endif
