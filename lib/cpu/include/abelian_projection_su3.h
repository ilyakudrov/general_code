#pragma once

#include "../include/matrix.h"
#include <vector>

using namespace std;

vector<vector<double>> make_angles_SU3(vector<su3> &conf);

void angles_project(std::vector<std::vector<double>> &angles);

void make_unitary(std::vector<std::vector<double>> &conf);