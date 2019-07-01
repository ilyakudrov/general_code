#ifndef __OBSERVABLES_H__
#define __OBSERVABLES_H__

#include "data.h"
#include "matrix.h"
#include "link.h"
#include "result.h"
#include "monopoles.h"

using namespace std;

double plaket_time(const data& conf);
double plaket_space(const data& conf);
double wilson(const data& conf, int R, int T);
double wilson_abelian(const data& conf, int R, int T);
void polyakov_abelian(const data& conf, double pol[2]);
double wilson_dir(const data& conf, int R, int T, int dir);
void fields(const vector<vector<matrix> >& schwinger_line, const vector<matrix>& plaket, const vector<matrix>& polyakov_loop, vector<vector<result> >& field1, vector<vector<result> >& field2, vector<result>& field3, int d, int D, int x_trans);
void field1_average(const vector<vector<matrix> >& schwinger_line, const vector<matrix>& plaket, const vector<matrix>& polyakov_loop, vector<vector<result> >& field1, int d, int D, int x_trans) ;
vector<vector<matrix> > calculate_schwinger_line(const data& conf, int d, int x_trans);
vector<matrix> calculate_plaket(const data& conf);
vector<matrix> calculate_polyakov_loop(const data& conf);
vector<matrix> calculate_wilson_loop(const data& conf, int R, int T);
vector<double> calculate_wilson_loop_tr(const data& conf, int R, int T);
double polyakov_loop_corelator(const data& conf, int D);
double plaket_correlator(const vector<matrix>& plaket, int dist);
double plaket_correlator_space(const vector<matrix>& plaket, int dist);
result wilson_plaket_correlator_electric_new(const data& conf, vector<double> wilson_loop_tr, int R, int T, int x_trans, int d_min, int d_max);
result wilson_plaket_correlator_electric_x_new(const data& conf, vector<double> wilson_loop_tr, int R, int T, int x_trans_min, int x_trans_max, int d);
result polyakov_plaket_correlator_electric(const data& conf, const data& smeared, int R, int x_trans, int d_min, int d_max);
result wilson_plaket_correlator_magnetic_new(const data& conf, vector<double> wilson_loop_tr, int R, int T, int x_trans, int d_min, int d_max);
result wilson_plaket_correlator_magnetic_x_new(const data& conf, vector<double> wilson_loop_tr, int R, int T, int x_trans_min, int x_trans_max, int d);
result polyakov_plaket_correlator_magnetic(const data& conf, const data& smeared, int R, int x_trans, int d_min, int d_max);
void push_result(result& values1_out, result& values2_out, result& values3_out, vector<vector<result> >& field1_values, vector<vector<result> >& field2_values, vector<result>& field3_values);
void calculate_and_push(data& conf, result& values1_out, result& values2_out, result& values3_out, vector<vector<result> >& field1_values, vector<vector<result> >& field2_values, vector<result>& field3_values, vector<vector<matrix> >& schwinger_line, vector<matrix>& plaket, vector<matrix>& polyakov_loop, int d, int D, int x_trans);
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
