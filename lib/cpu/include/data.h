//Data class

#ifndef __DATA_H__
#define __DATA_H__

#include <iostream>
#include <fstream>
#include <vector>
#include "matrix.h"
#include "math.h"

class data_matrix {
	public:
	std::vector<matrix> array;
	data_matrix();
	void read_float(char const* file_name);//read conf file of floats
	void read_float_abelian(char const* file_name);
	void read_double(char const* file_name);//read conf file of double
	void read_double_fortran(char const* file_name);
	void write_double(char const* file_name);//writes in file
	void write_float(char const* file_name);
	int eta(int mu, int x, int y, int z, int t);
	int sign(int x);
};

class data_double {
	public:
	std::vector<double> array;
	data_double();
	void read_float(char const* file_name); //read conf file of floats
	void read_float_fortran(char const* file_name);
	void write_float_fortran(char const* file_name);
	void read_float_fortran_convert_abelian(char const* file_name);
};
#endif
