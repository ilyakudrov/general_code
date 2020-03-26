#include "data.h"
#include "observables.h"
#include "link.h"
#include "matrix.h"
#include "result.h"
#include "eigen.h"

#include <stdio.h>
#include <iostream>
#include <ctime>
#include <cstring>
#include <cstdlib>

using namespace std;

int x_size = 32;
int y_size = 32;
int z_size = 32;
int t_size = 32;

int main(int argc, char* argv[]) {
    /*x_size = atof(argv[1]);
    y_size = atof(argv[2]);
    z_size = atof(argv[3]);
    t_size = atof(argv[4]);*/

    link1<matrix> link(x_size, y_size, z_size, t_size);
	data_matrix conf;
	data_double conf_abelian;
	char const *path1 = "../../confs/su2/time_32/mu0.00/conf_0001.fl";
	char const *path_abelian = "../../confs/su2/abelian/CON_MON_MAG_031.LAT";
	conf.read_float(path1);
	conf_abelian.read_float_fortran(path_abelian);
	cout<<wilson(conf_abelian.array, 10, 6)<<endl;
	cout<<plaket_time(conf_abelian.array)<<endl;
	cout<<plaket_space(conf_abelian.array)<<endl;
	result res_plaket;
	result res_wilson;
	result res_correlator;
	res_plaket.array = calculate_plaket_time_tr(conf_abelian.array);
	res_wilson.array = calculate_plaket_time_tr(conf_abelian.array);
	int d_min = 0;
	int d_max = 0;
	int x_trans = 0;
	int R = 10;
	int T = 6;
	res_correlator = wilson_plaket_correlator_electric_optimized(res_wilson.array, res_plaket.array, R, T, x_trans, d_min, d_max);
	cout<<res_correlator.array[0]<<endl;



	cout.precision(10);
	unsigned int start_time =  clock();

	cout<<"test plaket "<<plaket(conf.array)/2<<" right: 0.6769540066"<<endl;
	cout<<"test plaket_time "<<plaket_time(conf.array)/2<<" right: 0.6770628794"<<endl;
	cout<<"test plaket_space "<<plaket_space(conf.array)/2<<" right: 0.6768451339"<<endl;
	cout<<"test polyakov_loop "<<polyakov(conf.array)/2<<" right: -0.004586235468"<<endl;
	cout<<"test wilson_loop_R=10_T=6 "<<wilson(conf.array, 10, 6)<<" right: 0.001178588784"<<endl;

	unsigned int end_time = clock();
    unsigned int search_time = end_time - start_time;
    cout<<"working time: "<<search_time*1./CLOCKS_PER_SEC<<endl;
}
