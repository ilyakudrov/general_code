#include "../include/abelian_projection_su3.h"
#include "../include/matrix.h"
#include <cmath>
#include <vector>

using namespace std;

vector<vector<FLOAT>> make_angles_SU3(vector<su3_full> &conf) {

  vector<vector<FLOAT>> angles(3, vector<FLOAT>(conf.size()));

  vector<FLOAT> angle_tmp(3);

  for (int i = 0; i < conf.size(); i++) {

    double sum = 0;

    for (int j = 0; j < 3; j++) {
      angle_tmp[j] =
          atan2(conf[i].matrix[j][j].imag(), conf[i].matrix[j][j].real());
      sum += angle_tmp[j];
    }

    // while (sum >= M_PI) {
    //   sum -= 2 * M_PI;
    // }

    // while (sum < -M_PI) {
    //   sum += 2 * M_PI;
    // }

    for (int j = 0; j < 3; j++) {
      angles[j][i] = angle_tmp[j] - sum / 3;
    }
  }

  return angles;
}