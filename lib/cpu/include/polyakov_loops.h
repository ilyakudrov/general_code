#pragma once

#include "data.h"
#include "matrix.h"

#include <complex>
#include <map>
#include <vector>

template <class DataPattern, class MatrixType>
MatrixType inline get_polyakov_loop(
    const Data::LatticeData<DataPattern, MatrixType> &conf,
    DataPattern &data_pattern) {
  MatrixType A;
  int tmp = data_pattern.lat_coord[3];
  for (int i = 0; i < data_pattern.lat_dim[3]; i++) {
    A = A * conf[data_pattern.get_index_link(3)];
    data_pattern.move_forward(1, 3);
  }
  data_pattern.lat_coord[3] = tmp;
  return A;
}

template <class DataPattern, class MatrixType>
double polyakov_loop(const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  double polyakov_loop;
  MatrixType A;
#pragma omp parallel for collapse(3) private(A) firstprivate(data_pattern)     \
    reduction(+ : polyakov_loop)
  for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
    for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
      for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
        data_pattern.lat_coord = {x, y, z, 0};
        polyakov_loop += get_polyakov_loop(conf, data_pattern).tr_real();
      }
    }
  }
  return polyakov_loop / (data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                          data_pattern.lat_dim[2]);
}

template <class DataPattern, class MatrixType>
std::vector<typename MatrixType::trace_type> calculate_polyakov_loops_tr(
    const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<typename MatrixType::trace_type> polyakov_loops(
      data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
      data_pattern.lat_dim[2]);
#pragma omp parallel for collapse(3) firstprivate(data_pattern)
  for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
    for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
      for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
        data_pattern.lat_coord = {x, y, z, 0};
        polyakov_loops[data_pattern.get_index_site_spacial()] =
            get_polyakov_loop(conf, data_pattern).tr();
      }
    }
  }
  return polyakov_loops;
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType> calculate_polyakov_loops(
    const Data::LatticeData<DataPattern, MatrixType> &conf) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> polyakov_loops(data_pattern.lat_dim[0] *
                                         data_pattern.lat_dim[1] *
                                         data_pattern.lat_dim[2]);
#pragma omp parallel for collapse(3) firstprivate(data_pattern)
  for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
    for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
      for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
        data_pattern.lat_coord = {x, y, z, 0};
        polyakov_loops[data_pattern.get_index_site_spacial()] =
            get_polyakov_loop(conf, data_pattern);
      }
    }
  }
  return polyakov_loops;
}

template <class DataPattern, class MatrixType>
std::vector<double> polyakov_loop_correlator_color_average(
    const Data::LatticeData<DataPattern, MatrixType> &conf, int D_max) {
  std::vector<typename MatrixType::trace_type> polyakov_loops =
      calculate_polyakov_loops_tr(conf);
  DataPattern data_pattern(conf.lat_dim);
  int result_size = 2 * D_max + 1;
  int result_size1 = result_size * result_size;
  std::vector<double> correlator(result_size * result_size * result_size);
  double distance;
  int Dx_min;
  int place;
  int correlator_place;

  for (int Dz = -D_max; Dz <= D_max; Dz++) {
    for (int Dy = -D_max; Dy <= D_max; Dy++) {
      if (Dy != 0 && Dz != 0) {
        Dx_min = 0;
      } else {
        Dx_min = -D_max;
      }
      for (int Dx = Dx_min; Dx <= D_max; Dx++) {
        distance = sqrt(Dx * Dx + Dy * Dy + Dz * Dz);
        if (distance <= D_max) {
          correlator_place = (Dz + D_max) * result_size1 +
                             (Dy + D_max) * result_size + Dx + D_max;
          double correlator_tmp = 0;
#pragma omp parallel for collapse(2) private(place)                            \
    firstprivate(data_pattern, Dx, Dy, Dz) reduction(+ : correlator_tmp)
          for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
            for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
              for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
                data_pattern.lat_coord = {x, y, z, 0};
                place = data_pattern.get_index_site_spacial();
                data_pattern.move_forward(Dx, 0);
                data_pattern.move_forward(Dy, 1);
                data_pattern.move_forward(Dz, 2);
                correlator_tmp += std::real(
                    polyakov_loops[place] *
                    std::conj(
                        polyakov_loops[data_pattern.get_index_site_spacial()]));
              }
            }
          }
          correlator[correlator_place] = correlator_tmp;
        }
      }
    }
  }
  int size = data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
             data_pattern.lat_dim[2];
#pragma omp parallel for collapse(3)                                           \
    firstprivate(size, D_max, result_size1, result_size)
  for (int Dz = -D_max; Dz <= D_max; Dz++) {
    for (int Dy = -D_max; Dy <= D_max; Dy++) {
      for (int Dx = -D_max; Dx <= D_max; Dx++) {
        correlator[(Dz + D_max) * result_size1 + (Dy + D_max) * result_size +
                   Dx + D_max] /= size;
      }
    }
  }
  return correlator;
}

template <class DataPattern, class MatrixType>
std::vector<double> polyakov_loop_correlator_singlet(
    const Data::LatticeData<DataPattern, MatrixType> &conf, int D_max) {
  std::vector<MatrixType> polyakov_loops = calculate_polyakov_loops(conf);
  DataPattern data_pattern(conf.lat_dim);
  int result_size = 2 * D_max + 1;
  int result_size1 = result_size * result_size;
  std::vector<double> correlator(result_size * result_size * result_size);
  double distance;
  int Dx_min;
  int place;
  int correlator_place;

  for (int Dz = -D_max; Dz <= D_max; Dz++) {
    for (int Dy = -D_max; Dy <= D_max; Dy++) {
      if (Dy != 0 && Dz != 0) {
        Dx_min = 0;
      } else {
        Dx_min = -D_max;
      }
      for (int Dx = Dx_min; Dx <= D_max; Dx++) {
        distance = sqrt(Dx * Dx + Dy * Dy + Dz * Dz);
        if (distance <= D_max) {
          correlator_place = (Dz + D_max) * result_size1 +
                             (Dy + D_max) * result_size + Dx + D_max;
          double correlator_tmp = 0;
#pragma omp parallel for collapse(2) private(place)                            \
    firstprivate(data_pattern, Dx, Dy, Dz) reduction(+ : correlator_tmp)
          for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
            for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
              for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
                data_pattern.lat_coord = {x, y, z, 0};
                place = data_pattern.get_index_site_spacial();
                data_pattern.move(Dx, 0);
                data_pattern.move(Dy, 1);
                data_pattern.move(Dz, 2);
                correlator_tmp += polyakov_loops[place].multiply_conj_tr(
                    polyakov_loops[data_pattern.get_index_site_spacial()]);
              }
            }
          }
          correlator[correlator_place] = correlator_tmp;
        }
      }
    }
  }
  int size = data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
             data_pattern.lat_dim[2];
#pragma omp parallel for collapse(3)                                           \
    firstprivate(size, D_max, result_size1, result_size)
  for (int Dz = -D_max; Dz <= D_max; Dz++) {
    for (int Dy = -D_max; Dy <= D_max; Dy++) {
      for (int Dx = -D_max; Dx <= D_max; Dx++) {
        correlator[(Dz + D_max) * result_size1 + (Dy + D_max) * result_size +
                   Dx + D_max] /= size;
      }
    }
  }
  return correlator;
}

std::map<double, double>
polyakov_average_directions(const std::vector<double> &correlators, int D_max);