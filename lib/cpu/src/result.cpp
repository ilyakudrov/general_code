#include "result.h"

result::result(int size) {
	array = vector<double>(size);
}
result::result(){
	array = vector<double>(0);
}
void result::average(double a[2]) {
	if(array.size() == 0){
		a[0] = 0; a[1] = 1;
		return ;
	}
	double aver = 0;
	double err = 0;
	for (int i = 0; i < array.size(); i++) {
		aver += array[i] / array.size();
	}
	for (int i = 0; i < array.size(); i++) {
		err += (array[i] - aver)*(array[i] - aver) / array.size() / (array.size() - 1);
	}
	err = powf(err, 0.5);
	a[0] = aver;
	a[1] = err;
}
void result::read(char const* file_name, int size1){
	array.clear();
	array.reserve(size1);
	ifstream stream(file_name);
	if(!stream.read((char*) &array[0], size1 * sizeof(double))) cout<<"result::read error: "<<file_name<<endl;
	stream.close();
}
void result::read_float(char const* file_name, int size1){
	array.clear();
	ifstream stream(file_name);
	vector<float> vec(size1);
	stream.read((char*) &vec[0], size1 * sizeof(float));
	for(int i = 0;i < size1;i++){
		array.push_back((double)vec[i]);
	}
	stream.close();
}
void result::write(char const* file_name){
	ofstream stream(file_name);
	stream.write((char*) &array[0], array.size() * sizeof(double));
	stream.close();
}
double result::get_min() {
	double a = array[0];
	for (int i = 0; i < array.size(); i++) {
		if (a > array[i]) a = array[i];
	}
	return a;
}
double result::get_max() {
	double a = array[0];
	for (int i = 0; i < array.size(); i++) {
		if (a < array[i]) a = array[i];
	}
	return a;
}
void result::get_hist(int n, result& plaket, result& number) {
	double min = get_min();
	double max = get_max();
	int count = 0;
	double step = (max - min) / n;
	for (int i = 0; i < n; i++) {
		plaket.array.push_back(min + i * step);
		count = 0;
		for (int j = 0; j < array.size(); j++) {
			if (array[j] >= (min + i * step) && array[j] <= (min + (i + 1) * step)) count += 1;
		}
		number.array.push_back(count);
	}
}
double result::average_n(int n){
	int size = array.size();
	double a = 0;
	for(int i = 0;i < size;i++){
		if(i != n) a += array[i] / (size - 1);
	}
	return a;
}
void average_jack(double a[2], result& val1, result& val2, result& val3) {
	int size = val1.array.size();
	vector< double > vec(size);
	double aver1[2];
	double aver2[2];
	double aver3[2];
	val1.average(aver1);
	val2.average(aver2);
	val3.average(aver3);
	double b = (aver1[0] - aver2[0] / 2) / aver3[0];
	double sigma = 0;
	for (int i = 0; i < size; i++) {
		vec[i] = (val1.average_n(i) - val2.average_n(i)/ 2) / val3.average_n(i);
	}
	for (int i = 0; i < size; i++) {
		sigma += powf(vec[i] - b, 2);
	}
	sigma = powf(sigma * (size - 1) / size, 0.5);
	a[0] = b;
	a[1] = sigma;
}

void average_jack_wilson(double a[2], result& val1, result& val2, result& val3) {
	int size = val1.array.size();
	vector<double> vec(size);
	double aver1[2];
	double aver2[2];
	double aver3[2];
	val1.average(aver1);
	val2.average(aver2);
	val3.average(aver3);
	double b = aver1[0]/aver2[0] - aver3[0];
	for (int i = 0; i < size; i++) {
		vec[i] = val1.average_n(i)/val2.average_n(i) - val3.average_n(i);
	}
	double sigma = 0;
	for (int i = 0; i < size; i++) {
		sigma += (vec[i] - b) * (vec[i] - b);
	}
	sigma = powf(sigma * (size - 1) / size, 0.5);
	a[0] = b;
	a[1] = sigma;
}

void average_jack_sum(double a[2], result& val11, result& val12, result& val2, result& val3, result& val4) {
	int size = val11.array.size();
	vector<double> vec(size);
	double aver11[2];
	double aver12[2];
	double aver2[2];
	double aver3[2];
	double aver4[2];
	val11.average(aver11);
	val12.average(aver12);
	val2.average(aver2);
	val3.average(aver3);
	val4.average(aver4);
	double b = (aver11[0] + aver12[0])/aver2[0] - aver3[0] - aver4[0];
	for (int i = 0; i < size; i++) {
		vec[i] = (val11.average_n(i) + val12.average_n(i))/val2.average_n(i) - val3.average_n(i) - val4.average_n(i);
	}
	double sigma = 0;
	for (int i = 0; i < size; i++) {
		sigma += (vec[i] - b) * (vec[i] - b);
	}
	sigma = powf(sigma * (size - 1) / size, 0.5);
	a[0] = b;
	a[1] = sigma;
}

void average_jack_difference(double a[2], result& val11, result& val12, result& val2, result& val3, result& val4) {
	int size = val11.array.size();
	vector<double> vec(size);
	double aver11[2];
	double aver12[2];
	double aver2[2];
	double aver3[2];
	double aver4[2];
	val11.average(aver11);
	val12.average(aver12);
	val2.average(aver2);
	val3.average(aver3);
	val4.average(aver4);
	double b = (aver11[0] - aver12[0])/aver2[0] - aver3[0] + aver4[0];
	for (int i = 0; i < size; i++) {
		vec[i] = (val11.average_n(i) - val12.average_n(i))/val2.average_n(i) - val3.average_n(i) + val4.average_n(i);
	}
	double sigma = 0;
	for (int i = 0; i < size; i++) {
		sigma += (vec[i] - b) * (vec[i] - b);
	}
	sigma = powf(sigma * (size - 1) / size, 0.5);
	a[0] = b;
	a[1] = sigma;
}

void average_jackknife(double a[2], result& val1){
	int size = val1.array.size();
	vector< double > vec(size);
	double aver[2];
	val1.average(aver);
	double b = aver[0];
	double sigma = 0;
	for (int i = 0; i < size; i++) {
		vec[i] = val1.average_n(i);
	}
	for (int i = 0; i < size; i++) {
		sigma += (vec[i] - b) * (vec[i] - b);
	}
	sigma = powf(sigma * (size - 1) / size, 0.5);
	a[0] = b;
	a[1] = sigma;
}

double bootstrap_wilson(double aver[2], result& val1, result& val2, result& val3){
	int size = val1.array.size();
	int rand1 = 0;
	result res1(0);
	result res2(0);
	result res3(0);
	double aver1[2];
	double aver2[2];
	double aver3[2];
	for(int i = 0;i < size;i++){
		rand1 = rand()%size;
		res1.array.push_back(val1.array[rand1]);
		res2.array.push_back(val2.array[rand1]);
		res3.array.push_back(val3.array[rand1]);
	}
	res1.average(aver1);
	res2.average(aver2);
	res3.average(aver3);
	return aver1[0]/aver2[0]-aver3[0];
}

void average_bootstrap_wilson(double a[2], result& val1, result& val2, result& val3, int k) {
	vector<double> vec(k);
	double aver1[2];
	for (int i = 0; i < k; i++){
		vec[i] = bootstrap_wilson(aver1, val1, val2, val3);
	}
	double aver = 0;
	for(int i = 0;i < k;i++){
		aver += vec[i]/k;
	}
	double sigma = 0;
	for (int i = 0; i < k; i++) {
		sigma += (vec[i] - aver) * (vec[i] - aver);
	}
	sigma = powf(sigma / k, 0.5);
	a[0] = aver;
	a[1] = sigma;
}
