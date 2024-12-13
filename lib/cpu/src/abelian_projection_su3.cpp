#include "../include/abelian_projection_su3.h"
#include "../include/matrix.h"

#include <cmath>
#include <vector>

std::vector<std::vector<double>> make_angles_SU3(std::vector<su3> &conf) {

  std::vector<std::vector<double>> angles(3, std::vector<double>(conf.size()));

  std::vector<double> angle_tmp(3);

  for (int i = 0; i < conf.size(); i++) {

    double sum = 0;

    for (int j = 0; j < 3; j++) {
      angle_tmp[j] =
          atan2(conf[i].matrix[j][j].imag, conf[i].matrix[j][j].real);
      sum += angle_tmp[j];
    }

    while (sum >= M_PI) {
      sum -= 2 * M_PI;
    }

    while (sum < -M_PI) {
      sum += 2 * M_PI;
    }

    for (int j = 0; j < 3; j++) {
      angles[j][i] = angle_tmp[j] - sum / 3;
    }
  }

  return angles;
}

void angles_project(std::vector<std::vector<double>> &angles) {

  std::vector<double> angle_tmp(3);

  for (int i = 0; i < angles[0].size(); i++) {

    double sum = 0;

    for (int j = 0; j < 3; j++) {
      sum += angles[j][i];
    }

    while (sum >= M_PI) {
      sum -= 2 * M_PI;
    }

    while (sum < -M_PI) {
      sum += 2 * M_PI;
    }

    for (int j = 0; j < 3; j++) {
      angles[j][i] -= sum / 3;
    }
  }
}

void make_unitary(std::vector<std::vector<double>> &angles) {
  double sum;
  for (int i = 0; i < angles[0].size(); i++) {
    sum = 0;
    for (int j = 0; j < 3; j++) {
      sum += angles[j][i];
    }

    while (sum >= M_PI) {
      sum -= 2 * M_PI;
    }

    while (sum < -M_PI) {
      sum += 2 * M_PI;
    }

    for (int j = 0; j < 3; j++) {
      angles[j][i] -= sum / 3;
    }
  }
}