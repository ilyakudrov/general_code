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

int x_size;
int y_size;
int z_size;
int t_size;

int main(int argc, char* argv[]) {
    x_size = atof(argv[1]);
    y_size = atof(argv[2]);
    z_size = atof(argv[3]);
    t_size = atof(argv[4]);

    double alpha1 = 1;
	double alpha2 = 1;
    double alpha3 = 0.5;
    double alpha_APE = 0.3;
    double stout_rho = 0.15;

    vector<vector<matrix> > smearing_first(6, vector<matrix>(1));
    vector<vector<matrix> > smearing_second(4, vector<matrix>(1));

    link1 link(x_size, y_size, z_size, t_size);
	data conf;
	data smeared;
	char path1[] = "../../confs/su2/time_32/mu0.00/conf_0001.fl";
	conf.read_float(path1);
	double aver[2];

	cout.precision(10);
	matrix A;
	unsigned int start_time =  clock();

	smearing_first = link.smearing_first_full(conf, alpha3);
	smearing_second = link.smearing_second_full(conf, smearing_first, alpha2);
	smeared.array = link.smearing_HYP(conf, smearing_second, alpha1);

	for(int i = 0;i < 4;i++){
		A = smeared.array[i];
		cout<<"smearing_HYP test "<<A.a0<<" "<<A.a1<<" "<<A.a2<<" "<<A.a3<<endl;
	}
	cout<<"right:"<<endl;
	cout<<"-0.110943 0.629543 -0.583767 0.500582"<<endl;
	cout<<"0.318657 -0.218438 -0.415944 0.823245"<<endl;
	cout<<"0.174894 0.667862 0.674797 0.260808"<<endl;
	cout<<"-0.704066 -0.13241 0.453476 -0.530206"<<endl;

	smeared.array = link.smearing_HYP_refresh(conf, alpha1, alpha2, alpha3);

	for(int i = 0;i < 4;i++){
		A = smeared.array[i];
		cout<<"smearing_HYP_refresh test "<<A.a0<<" "<<A.a1<<" "<<A.a2<<" "<<A.a3<<endl;
	}
	cout<<"right:"<<endl;
	cout<<"-0.110943 0.629543 -0.583767 0.500582"<<endl;
	cout<<"0.318657 -0.218438 -0.415944 0.823245"<<endl;
	cout<<"0.174894 0.667862 0.674797 0.260808"<<endl;
	cout<<"-0.704066 -0.13241 0.453476 -0.530206"<<endl;

	smeared.array = link.smearing_APE(conf, alpha_APE);

	for(int i = 0;i < 4;i++){
		A = smeared.array[i];
		cout<<"smearing_APE test "<<A.a0<<" "<<A.a1<<" "<<A.a2<<" "<<A.a3<<endl;
	}
	cout<<"right:"<<endl;
	cout<<"-0.0833883 0.641689 -0.604086 0.465147"<<endl;
	cout<<"0.336007 -0.177589 -0.426239 0.820903"<<endl;
	cout<<"0.125318 0.626699 0.730398 0.240963"<<endl;
	cout<<"-0.632818 0.13761 0.745239 0.158818"<<endl;

	smeared.array = link.smearing_stout(conf, stout_rho);

	for(int i = 0;i < 4;i++){
		A = smeared.array[i];
		cout<<"smearing_stout test "<<A.a0<<" "<<A.a1<<" "<<A.a2<<" "<<A.a3<<endl;
	}
	cout<<"right:"<<endl;
	cout<<"-0.0185375 0.684286 -0.660024 0.30948"<<endl;
	cout<<"0.363087 -0.0991993 -0.575234 0.726246"<<endl;
	cout<<"0.172742 0.633665 0.73173 0.182208"<<endl;
	cout<<"-0.816237 -0.0404862 0.571983 -0.0703786"<<endl;

	unsigned int end_time = clock();
    unsigned int search_time = end_time - start_time;
    cout<<"working time: "<<search_time*1./CLOCKS_PER_SEC<<endl;
}