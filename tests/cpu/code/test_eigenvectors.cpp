#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/eigen.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"

#include <cstring>
#include <ctime>
#include <iostream>
#include <vector>

int x_size;
int y_size;
int z_size;
int t_size;
int size1;
int size2;

int main(int argc, char *argv[]) {
  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  std::cout.precision(17);

  x_size = 40;
  y_size = 40;
  z_size = 40;
  t_size = 40;
  size1 = x_size * y_size;
  size2 = x_size * y_size * z_size;

  link1 link(x_size, y_size, z_size, t_size);
  Data::data<su2> conf_qc2dstag;
  std::string path_qc2dstag = "../../confs/qc2dstag/40^4/mu0.05/s0/CONF0201";

  double mu_q = 0.05;
  int vector_size = x_size * y_size * z_size * t_size * 2;

  conf_qc2dstag.read_double_qc2dstag(path_qc2dstag);

  std::vector<su2> matrix_staggered =
      make_matrix_staggered(conf_qc2dstag.array, mu_q);

  std::vector<complex> vector_input(vector_size);

  int place;
  for (int t = 0; t < t_size; t++) {
    for (int z = 0; z < z_size; z++) {
      for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
          place = t * 2 * x_size * y_size * z_size + z * 2 * x_size * y_size +
                  y * 2 * x_size + x * 2;
          vector_input[place].re =
              0.1 * (x + 1) + 0.2 * (y + 1) + 0.4 * (z + 1) + 0.5 * (t + 1);
          vector_input[place].im =
              1 + 0.6 * (x + 1) + 0.7 * (y + 1) + 0.8 * (z + 1) + 0.9 * (t + 1);
          vector_input[place + 1].re =
              1 + 0.1 * (x + 1) + 0.2 * (y + 1) + 0.4 * (z + 1) + 0.5 * (t + 1);
          vector_input[place + 1].im =
              0.6 * (x + 1) + 0.7 * (y + 1) + 0.8 * (z + 1) + 0.9 * (t + 1);
        }
      }
    }
  }

  std::vector<complex> vector_output =
      matrix_multiplication_staggered(matrix_staggered, vector_input);

  for (int i = 0; i < 10; i++) {
    std::cout << vector_output[i].re << " i(" << vector_output[i].im << ")"
              << std::endl;
  }
  double norm = 0;
  for (int i = 0; i < vector_output.size(); i++) {
    norm += vector_output[i].re * vector_output[i].re +
            vector_output[i].im * vector_output[i].im;
  }
  std::cout << "norm " << norm << std::endl;
}