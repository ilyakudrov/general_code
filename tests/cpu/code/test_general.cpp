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
	unsigned int start_time;
	unsigned int end_time;
    unsigned int search_time;
    link1<matrix> link(x_size, y_size, z_size, t_size);
	data_matrix conf;
	data_double conf_abelian;
	data_matrix conf_offd;
	char const *path1 = "../../confs/su2/time_32/mu0.00/conf_0001.fl";
	char const *path_abelian = "../../confs/su2/abelian/CON_MON_MAG_031.LAT";
	char const *path_offd = "/home/ilya/lattice/conf/offd/CON_OFF_MAG_033.LAT";
	conf.read_float(path1);
	conf_abelian.read_float_fortran(path_abelian);
	conf_offd.read_double_fortran(path_offd);

	for(int i = 0;i < 5;i++){
		cout<<conf_offd.array[i]<<endl;
		cout<<conf_offd.array[i].module()<<endl;
	}


	link1<double> link_abelian(x_size, y_size, z_size, t_size);
	cout<<"first wilson "<<link_abelian.wilson_loop(conf_abelian.array, 1, 1)<<endl;
	start_time =  clock();
	cout<<"wilson_abelian aver "<<wilson(conf_abelian.array, 1, 1)<<endl;
	end_time = clock();
    search_time = end_time - start_time;
    cout<<"wilson time: "<<search_time*1./CLOCKS_PER_SEC<<endl;
	cout<<plaket_time(conf_abelian.array)<<endl;
	cout<<plaket_space(conf_abelian.array)<<endl;

	int R = 1;
	int T = 1;
	result res_plaket_time;
	result res_plaket_space;
	result res_wilson;
	result res_correlator_electric;
	result res_correlator_magnetic;
	res_plaket_time.array = calculate_plaket_time_tr(conf_abelian.array);
	res_plaket_space.array = calculate_plaket_space_tr(conf_abelian.array);
	res_wilson.array = calculate_wilson_loop_tr(conf_abelian.array, R, T);
	double aver[2];
	res_plaket_time.average(aver);
	cout<<"plaket time first "<<res_plaket_time.array[0]<<endl;
	cout<<"plaket_time aver "<<aver[0]<<endl;
	res_wilson.average(aver);
	cout<<"wilson first "<<res_wilson.array[0]<<endl;
	cout<<"wilson aver "<<aver[0]<<endl;

	int d_min = -10;
	int d_max = 10;
	int x_trans = 0;
	res_correlator_electric = wilson_plaket_correlator_electric(res_wilson.array, res_plaket_time.array, R, T, x_trans, d_min, d_max);
	res_correlator_magnetic = wilson_plaket_correlator_magnetic(res_wilson.array, res_plaket_space.array, R, T, x_trans, d_min, d_max);
	// for(int i = 0;i < res_correlator.array.size();i++){
	// 	cout<<res_correlator.array[i]<<endl;
	// }



	cout.precision(10);
	start_time =  clock();

	cout<<"test plaket "<<plaket(conf.array)/2<<" right: 0.6769540066"<<endl;
	cout<<"test plaket_time "<<plaket_time(conf.array)/2<<" right: 0.6770628794"<<endl;
	cout<<"test plaket_space "<<plaket_space(conf.array)/2<<" right: 0.6768451339"<<endl;
	cout<<"test polyakov_loop "<<polyakov(conf.array)/2<<" right: -0.004586235468"<<endl;
	cout<<"test wilson_loop_R=10_T=6 "<<wilson(conf.array, 10, 6)<<" right: 0.001178588784"<<endl;

	end_time = clock();
    search_time = end_time - start_time;
    cout<<"working time: "<<search_time*1./CLOCKS_PER_SEC<<endl;
}
