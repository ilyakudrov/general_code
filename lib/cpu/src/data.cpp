#include "data.h"
//#include "link.h"

extern int x_size;
extern int y_size;
extern int z_size;
extern int t_size;

using namespace std;

data_matrix::data_matrix() {
	array.reserve(4 * x_size * y_size * z_size * t_size);
}
void data_matrix::read_float(char const* file_name){
	int data_size1 = 4 * x_size * y_size * z_size * t_size;
	array.clear();
	ifstream stream(file_name);
	vector<float> v(data_size1 * 4);
	if(!stream.read((char*) &v[0], data_size1 * 4 * sizeof(float))) cout<<"read_float error: "<<file_name<<endl;
	matrix A;
	for (int i = 0; i < data_size1; i++) {
		A.a0 = (double)v[i * 4];
		A.a1 = (double)v[i * 4 + 1];
		A.a2 = (double)v[i * 4 + 2];
		A.a3 = (double)v[i * 4 + 3];
		array.push_back(A);
	}
	stream.close();
}
void data_matrix::read_float_abelian(char const* file_name){
	int data_size1 = 4 * x_size * y_size * z_size * t_size;
	array.clear();
	ifstream stream(file_name);
	vector<float> v(data_size1 * 4 + 2);
	if(!stream.read((char*) &v[0], (data_size1 * 4 + 2) * sizeof(float))) cout<<"read_float error: "<<file_name<<endl;
	matrix A;
	for (int i = 0; i < data_size1; i++) {
		A.a0 = (double)v[i * 4 + 1];
		A.a1 = (double)v[i * 4 + 2];
		A.a2 = (double)v[i * 4 + 3];
		A.a3 = (double)v[i * 4 + 4];
		array.push_back(A);
	}
	stream.close();
}
void data_matrix::read_double(char const* file_name){
	int data_size1 = 4 * x_size * y_size * z_size * t_size;
	array.clear();
	ifstream stream(file_name);
	vector<double> v(data_size1 * 4);
	if(!stream.read((char*) &v[0], data_size1 * 4 * sizeof(double))) cout<<"read_double error: "<<file_name<<endl;
	matrix A;
	for (int i = 0; i < data_size1; i++) {
		A.a0 = v[i * 4];
		A.a1 = v[i * 4 + 1];
		A.a2 = v[i * 4 + 2];
		A.a3 = v[i * 4 + 3];
		array.push_back(A);
	}
	stream.close();
}
void data_matrix::read_double_fortran(char const* file_name){
	int data_size1 = 4 * x_size * y_size * z_size * t_size;
	array.clear();
	ifstream stream(file_name);
	vector<double> v(data_size1 * 4);
	stream.ignore(4);
	if(!stream.read((char*) &v[0], (data_size1 * 4) * sizeof(double))) cout<<"read_float error: "<<file_name<<endl;
	matrix A;
	for (int i = 0; i < data_size1; i++) {
		A.a0 = v[i * 4];
		A.a1 = v[i * 4 + 1];
		A.a2 = v[i * 4 + 2];
		A.a3 = v[i * 4 + 3];
		array.push_back(A);
	}
	stream.close();
}
void data_matrix::write_float(char const* file_name) {
	int data_size1 = 4 * x_size * y_size * z_size * t_size;
	ofstream stream(file_name);
	vector<float> v(data_size1 * 4);
	for (int i = 0; i < data_size1; i++) {
		v[i * 4] = (float)array[i].a0;
		v[i * 4 + 1] = (float)array[i].a1;
		v[i * 4 + 2] = (float)array[i].a2;
		v[i * 4 + 3] = (float)array[i].a3;
	}
	if(!stream.write((char*) &v[0], data_size1 * 4 * sizeof(float))) cout<<"write_float error: "<<file_name<<endl;
	stream.close();
}
void data_matrix::write_double(char const* file_name) {
	int data_size1 = 4 * x_size * y_size * z_size * t_size;
	ofstream stream(file_name);
	vector<double> v(data_size1 * 4);
	for (int i = 0; i < data_size1; i++) {
		v[i * 4] = array[i].a0;
		v[i * 4 + 1] = array[i].a1;
		v[i * 4 + 2] = array[i].a2;
		v[i * 4 + 3] = array[i].a3;
	}
	if(!stream.write((char*) &v[0], data_size1 * 4 * sizeof(double))) cout<<"write_double error: "<<file_name<<endl;
	stream.close();
}
int data_matrix::eta(int mu, int x, int y, int z, int t) {
	int a = 0;
	int coordinate[4] = { x, y, z, t };
	for (int i = 0; i < mu - 1; i++) {
		a += coordinate[i];
	}
	return sign(a);
}
int data_matrix::sign(int x) {
	if (x % 2 == 0) return 1;
	else if (x % 2 == 1) return -1;
	else return 0;
}

data_double::data_double() {
	array.reserve(4 * x_size * y_size * z_size * t_size);
}

void data_double::read_float(char const* file_name){
	int data_size1 = 4 * x_size * y_size * z_size * t_size;
	array.clear();
	ifstream stream(file_name);
	vector<float> v(data_size1 * 4);
	if(!stream.read((char*) &v[0], data_size1 * sizeof(float))) cout<<"read_float error: "<<file_name<<endl;
	for (int i = 0; i < data_size1; i++) {
		array.push_back((double)v[i]);
	}
	stream.close();
}

void data_double::read_float_fortran(char const* file_name){
	int data_size1 = 4 * x_size * y_size * z_size * t_size;
	array.clear();
	ifstream stream(file_name);
	vector<float> v(data_size1 + 2);
	if(!stream.read((char*) &v[0], (data_size1 + 2) * sizeof(float))) cout<<"read_float error: "<<file_name<<endl;
	for (int i = 0; i < data_size1; i++) {
		array.push_back((double)v[i + 1]);
	}
	stream.close();
}

void data_double::write_float_fortran(char const* file_name) {
	int data_size1 = 4 * x_size * y_size * z_size * t_size;
	ofstream stream(file_name);
	vector<float> v(data_size1 + 2);
	for (int i = 0; i < data_size1; i++) {
		v[i + 1] = (float) array[i];
	}
	if(!stream.write((char*) &v[0], (data_size1 + 2) * sizeof(float))) cout<<"write_double error: "<<file_name<<endl;
	stream.close();
}

void data_double::read_float_fortran_convert_abelian(char const* file_name){
	int data_size1 = 4 * x_size * y_size * z_size * t_size;
	array.clear();
	ifstream stream(file_name);
	vector<float> v(data_size1 * 4 + 2);
	if(!stream.read((char*) &v[0], (data_size1 * 4 + 2) * sizeof(float))) cout<<"read_float error: "<<file_name<<endl;
	for (int i = 0; i < data_size1; i++) {
		array.push_back((double)atan2(v[i * 4 + 4], v[i * 4 + 1]));
	}
	stream.close();
}
