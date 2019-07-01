#ifndef __RESULT_H__
#define __RESULT_H__

#include <fstream>
#include <vector>
#include <math.h>
#include <iostream>
#include <ctime>
#include <cstdlib>

using namespace std;

class result {
	public:
	vector<double> array;
	result(int size);
	result();
	void average(double a[2]);
	void write(char output_file[]);
	void read(char file_name[], int size1);
	void read_float(char file_name[], int size1);
	double get_min();
	double get_max();
	void get_hist(int n, result& plaket, result& number);
	double average_n(int n);
};
void average_jack(double a[2], result& val1, result& val2, result& val3);
void average_jack_wilson(double a[2], result& val1, result& val2, result& val3);
void average_jack_sum(double a[2], result& val11, result& val12, result& val2, result& val3, result& val4);
void average_jack_difference(double a[2], result& val11, result& val12, result& val2, result& val3, result& val4);
void average_jackknife(double a[2], result& val1);
double bootstrap_wilson(double aver[2], result& val1, result& val2, result& val3);
void average_bootstrap_wilson(double a[2], result& val1, result& val2, result& val3, int k);
#endif
