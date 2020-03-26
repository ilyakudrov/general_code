#ifndef __LINK_H__
#define __LINK_H__

using namespace std;

extern int x_size;
extern int y_size;
extern int z_size;
extern int t_size;

#include <iostream>
#include "stdlib.h"
#include "math.h"
#include <vector>
#include "matrix.h"
#include "data.h"

class data;

template<typename T>
class link1 {
	public:
	int lattice_size[4];
	int coordinate[4];//0..lattise_size[i]-1
	int direction;

	link1(int lattice_size_x, int lattice_size_y, int lattice_size_z, int lattice_size_t);
	link1(int dir, int lattice_size_x, int lattice_size_y, int lattice_size_z, int lattice_size_t);
	link1(int coord[4], int dir, int lattice_size_x, int lattice_size_y, int lattice_size_z, int lattice_size_t);
	link1(int x, int y, int z, int t, int dir, int lattice_size_x, int lattice_size_y, int lattice_size_z, int lattice_size_t);
	link1(int x, int y, int z, int t, int lattice_size_x, int lattice_size_y, int lattice_size_z, int lattice_size_t);
	link1(const link1& link);
	link1();

	void print_link();
	void move(int dir, int step);
	void go(int x, int y, int z, int t);
	void move_dir(int dir);//changes the direction
	int get_place();
	int get_place1();
	matrix get_matrix(const vector<matrix>& vec); //works with negative directions(takes inverse matrix)
	double get_matrix(const vector<double>& vec);
	matrix get_matrix1(matrix* vec);
	double border_sign(int mu);
	double get_angle_abelian(const vector<matrix>& vec);
	T schwinger_line(const vector<T>& array, int d, int dir, int x);//link is attached to "left" source and directed to the plaket
	T plaket(const vector<T>& array);//Directed to the direction of the field
	T plaket_mu(const vector<T>& array, int mu);//mu is the second direction
	double plaket_abelian_mu(const vector<matrix>& array, int mu);
	T plaket_implement4(const vector<T>& array, int mu);
	double plaket_abelian_implement4(const vector<matrix>& array, int mu);
	T plaket_implement2(const vector<T>& array, int mu);
	T polyakov_loop(const vector<T>& array);//attached to where it loops and directed to the temporal direction
	double polyakov_loop_abelian(const vector<matrix>& array);
	T wilson_loop(const vector<T>& array, int r, int t);
	double wilson_loop_abelian(const vector<matrix>& array, int r, int t);
	//first numeratoк
	//d - distance between "left" source and plaket
	//D - distance between sources
	double field1(const vector<vector<matrix> >& schwinger_line, const vector<matrix>& plaket, const vector<matrix>& polyakov_loop, int d, int D, int dir, int x);//Link is attached to the "left" source, dir points to the plaket
	//second numerator
	double field2(const vector<matrix>& plaket, const vector<matrix>& polyakov_loop, int d, int D, int dir, int x);//attached to the "left" source, dir points to the plaket
	//denominator
	double field3(const vector<matrix>& polyakov_loop, int D, int x);//attached to the "left" source and points to another
	matrix staples_first(const vector<matrix>& vec, int eta);
	matrix staples_second(const vector<vector<matrix> >& smearing_first, int eta, int nu);
	matrix staples_second_refresh(const vector<matrix>& vec, int eta, int nu, double alpha3);//staples for refreshing algorythm(refresh link every step)
	matrix staples_third(const vector<vector<matrix> >& smearing_second, int eta);
	matrix staples_third_refresh(const vector<matrix>& vec, int eta, double alpha2, double alpha3);
	vector<matrix> smearing_first(const data_matrix& conf, double alpha3, int nu, int rho);
	vector<vector<matrix> > smearing_first_full(const data_matrix& conf, double alpha3);
	vector<matrix> smearing_second(const data_matrix& conf, vector<vector<matrix> >& smearing_first, double alpha2, int nu);
	vector<vector<matrix> > smearing_second_full(const data_matrix& conf, vector<vector<matrix> >& smearing_first, double alpha2);
	vector<matrix> smearing_HYP(data_matrix& conf, vector<vector<matrix> >& smearing_second, double alpha1);
	inline int position_first(int a, int b);
	vector<matrix> smearing_APE(data_matrix& conf, double alpha_APE);
	matrix smearing_first_refresh(const vector<matrix>& vec, int nu, int rho, double alpha3);//refresh link every step
	matrix smearing_second_refresh(const vector<matrix>& vec, int nu,  double alpha2, double alpha3);//refresh link every step
	vector<matrix> smearing_HYP_refresh(data_matrix& conf, double alpha1, double alpha2, double alpha3);//refresh link every step
	vector<matrix> smearing_APE_refresh(data_matrix& conf, double alpha_APE);//refresh link every step
	vector<matrix> smearing_stout(data_matrix& conf, double rho);
	matrix stout_factor(data_matrix& conf, double rho);
	matrix stout_omega(data_matrix& conf, double rho);
	void gauge_transform(data_matrix& conf);
	//monopoles
	double monopole_plaket(data_matrix& conf, int i, int j);//monopole_plaket
	double get_current(data_matrix& conf);
	int current_test(double* J);
};
#endif
