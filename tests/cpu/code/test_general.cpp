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

    link1 link(x_size, y_size, z_size, t_size);
	data conf;
	char const *path1 = "../../confs/su2/time_32/mu0.00/conf_0001.fl";
	conf.read_float(path1);
	double aver[2];

	cout.precision(10);
	unsigned int start_time =  clock();

	cout<<"test plaket "<<plaket(conf)/2<<" right: 0.6769540066"<<endl;
	cout<<"test plaket_time "<<plaket_time(conf)/2<<" right: 0.6770628794"<<endl;
	cout<<"test plaket_space "<<plaket_space(conf)/2<<" right: 0.6768451339"<<endl;
	cout<<"test polyakov_loop "<<polyakov(conf)/2<<" right: -0.004586235468"<<endl;
	cout<<"test wilson_loop_R=10_T=6 "<<wilson(conf, 10, 6)<<" right: 0.001178588784"<<endl;

	unsigned int end_time = clock();
    	unsigned int search_time = end_time - start_time;
    	cout<<"working time: "<<search_time*1./CLOCKS_PER_SEC<<endl;
}
