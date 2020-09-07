#ifndef __RESULT_H__
#define __RESULT_H__

#ifdef DOUBLE
#define FLOAT double
#else
#define FLOAT float
#endif

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

class result {
public:
  vector<FLOAT> array;
  result(int size);
  result();
  void average(FLOAT a[2]);
  void write(char const *file_name);
  void read(char const *file_name, int size1);
  void read_float(char const *file_name, int size1);
  FLOAT get_min();
  FLOAT get_max();
  void get_hist(int n, result &plaket, result &number);
  FLOAT average_n(int n);
};
void average_jack(FLOAT a[2], result &val1, result &val2, result &val3);
void average_jack_wilson(FLOAT a[2], result &val1, result &val2, result &val3);
void average_jack_sum(FLOAT a[2], result &val11, result &val12, result &val2,
                      result &val3, result &val4);
void average_jack_difference(FLOAT a[2], result &val11, result &val12,
                             result &val2, result &val3, result &val4);
void average_jackknife(FLOAT a[2], result &val1);
FLOAT bootstrap_wilson(FLOAT aver[2], result &val1, result &val2, result &val3);
void average_bootstrap_wilson(FLOAT a[2], result &val1, result &val2,
                              result &val3, int k);
#endif
