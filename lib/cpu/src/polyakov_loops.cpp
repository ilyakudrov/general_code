#include "../include/polyakov_loops.h"
#include "../include/indexing.h"
#include "../include/link.h"
#include "../include/matrix.h"

#include <algorithm>
#include <complex>
#include <map>
#include <omp.h>
#include <vector>

extern int x_size;
extern int y_size;
extern int z_size;
extern int t_size;
extern int size1;
extern int size2;

#define DATA_SIZE 4 * x_size *y_size *z_size *t_size
#define PLACE3_DIR                                                             \
  (t) * 3 * x_size *y_size *z_size + (z) * 3 * x_size *y_size +                \
      (y) * 3 * x_size + (x) * 3 + dir
#define PLACE3_LINK_NODIR                                                      \
  (link.coordinate[3]) * 3 * x_size *y_size *z_size +                          \
      (link.coordinate[2]) * 3 * x_size *y_size +                              \
      (link.coordinate[1]) * 3 * x_size + (link.coordinate[0]) * 3
#define PLACE1_NODIR                                                           \
  (t) * x_size *y_size *z_size + (z) * x_size *y_size + (y) * x_size + (x)
#define PLACE_PLAKET_TIME                                                      \
  (link.coordinate[3]) * 3 * x_size *y_size *z_size +                          \
      (link.coordinate[2]) * 3 * x_size *y_size +                              \
      (link.coordinate[1]) * 3 * x_size + (link.coordinate[0]) * 3 +           \
      link.direction
#define PLACE_PLAKET_SPACE                                                     \
  (link.coordinate[3]) * 3 * x_size *y_size *z_size +                          \
      (link.coordinate[2]) * 3 * x_size *y_size +                              \
      (link.coordinate[1]) * 3 * x_size + (link.coordinate[0]) * 3
#define PLACE1_LINK_NODIR                                                      \
  (link.coordinate[3]) * x_size *y_size *z_size +                              \
      (link.coordinate[2]) * x_size *y_size + (link.coordinate[1]) * x_size +  \
      (link.coordinate[0])

#define SPACE_ITER_START                                                       \
  for (int t = 0; t < t_size; t++) {                                           \
    for (int z = 0; z < z_size; z++) {                                         \
      for (int y = 0; y < y_size; y++) {                                       \
        for (int x = 0; x < x_size; x++) {                                     \
          link.go(x, y, z, t);                                                 \
          link.update(0);                                                      \
          link.update(1);                                                      \
          link.update(2);                                                      \
          link.update(3);

#define SPACE_ITER_END                                                         \
  }                                                                            \
  }                                                                            \
  }                                                                            \
  }

#define ITER_START_ZYX                                                         \
  for (int z = 0; z < z_size; z++) {                                           \
    for (int y = 0; y < y_size; y++) {                                         \
      for (int x = 0; x < x_size; x++) {                                       \
        link.go(x, y, z, 0);                                                   \
        link.update(0);                                                        \
        link.update(1);                                                        \
        link.update(2);

#define ITER_END_3                                                             \
  }                                                                            \
  }                                                                            \
  }

#define SPACE_ITER_START_3D                                                    \
  for (int z = 0; z < z_size; z++) {                                           \
    for (int y = 0; y < y_size; y++) {                                         \
      for (int x = 0; x < x_size; x++) {                                       \
        link.go(x, y, z, 0);                                                   \
        link.update(0);                                                        \
        link.update(1);                                                        \
        link.update(2);                                                        \
        link.update(3);

#define SPACE_ITER_END_3D                                                      \
  }                                                                            \
  }                                                                            \
  }

#pragma omp declare reduction(                                                 \
        vec_double_plus : std::vector<double> : std::transform(                \
                omp_out.begin(), omp_out.end(), omp_in.begin(),                \
                    omp_out.begin(), std::plus<double>()))                     \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

template <class T>
T get_polyakov_loop(const std::vector<T> &array, int index, int size2,
                    int length) {
  T A;
  for (int i = 0; i < length; i++) {
    A = A * array[index];
    index += size2;
  }
  return A;
}

template <class T> double polyakov_loop(const std::vector<T> &array) {
  std::vector<int> lat_dim = {x_size, y_size, z_size, t_size};
  std::vector<int> lat_coord = {0, 0, 0, 0};
  double polyakov_loop;
  T A;
  int index;
#pragma omp parallel for collapse(3) private(A, index)                         \
    firstprivate(lat_coord, lat_dim, size2) reduction(+ : polyakov_loop)
  for (int z = 0; z < lat_dim[2]; z++) {
    for (int y = 0; y < lat_dim[1]; y++) {
      for (int x = 0; x < lat_dim[0]; x++) {
        lat_coord[0] = x;
        lat_coord[1] = y;
        lat_coord[2] = z;
        polyakov_loop +=
            get_polyakov_loop(array, get_index_matrix(lat_coord, 3), size2 * 4,
                              lat_dim[3])
                .tr();
      }
    }
  }
  return polyakov_loop / (x_size * y_size * z_size);
}

template <class T>
double polyakov_loop_parallel(const std::vector<std::vector<T>> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  int place;
  double polyakov_loop = 0;
#pragma omp parallel for collapse(3) private(link, place)                      \
    reduction(+ : polyakov_loop)
  ITER_START_ZYX
  link.move_dir(3);
  place = link.place / 4;
  polyakov_loop += link.polyakov_loop(array).tr();
  ITER_END_3
  return polyakov_loop / (x_size * y_size * z_size);
}

template <class T>
std::vector<double> calculate_polyakov_loops_tr(const std::vector<T> &array) {
  std::vector<int> lat_dim = {x_size, y_size, z_size, t_size};
  std::vector<int> lat_coord = {0, 0, 0, 0};
  std::vector<double> polyakov_loops(x_size * y_size * z_size);

#pragma omp parallel for collapse(3) firstprivate(lat_coord, lat_dim, size2)
  for (int z = 0; z < lat_dim[2]; z++) {
    for (int y = 0; y < lat_dim[1]; y++) {
      for (int x = 0; x < lat_dim[0]; x++) {
        lat_coord[0] = x;
        lat_coord[1] = y;
        lat_coord[2] = z;
        polyakov_loops[get_index_site(lat_coord)] =
            get_polyakov_loop(array, get_index_matrix(lat_coord, 3), size2 * 4,
                              lat_dim[3])
                .tr();
      }
    }
  }
  return polyakov_loops;
}

std::vector<std::complex<double>>
calculate_polyakov_loops_tr(const std::vector<su3> &array) {
  std::vector<int> lat_dim = {x_size, y_size, z_size, t_size};
  std::vector<int> lat_coord = {0, 0, 0, 0};
  std::vector<std::complex<double>> polyakov_loops(x_size * y_size * z_size);

#pragma omp parallel for collapse(3) firstprivate(lat_coord, lat_dim, size2)
  for (int z = 0; z < lat_dim[2]; z++) {
    for (int y = 0; y < lat_dim[1]; y++) {
      for (int x = 0; x < lat_dim[0]; x++) {
        lat_coord[0] = x;
        lat_coord[1] = y;
        lat_coord[2] = z;
        polyakov_loops[get_index_site(lat_coord)] =
            get_polyakov_loop(array, get_index_matrix(lat_coord, 3), size2 * 4,
                              lat_dim[3])
                .tr_complex();
      }
    }
  }
  return polyakov_loops;
}

std::vector<std::complex<double>>
calculate_polyakov_loops_tr(const std::vector<su3_abelian> &array) {
  std::vector<int> lat_dim = {x_size, y_size, z_size, t_size};
  std::vector<int> lat_coord = {0, 0, 0, 0};
  std::vector<std::complex<double>> polyakov_loops(x_size * y_size * z_size);

#pragma omp parallel for collapse(3) firstprivate(lat_coord, lat_dim, size2)
  for (int z = 0; z < lat_dim[2]; z++) {
    for (int y = 0; y < lat_dim[1]; y++) {
      for (int x = 0; x < lat_dim[0]; x++) {
        lat_coord[0] = x;
        lat_coord[1] = y;
        lat_coord[2] = z;
        polyakov_loops[get_index_site(lat_coord)] =
            get_polyakov_loop(array, get_index_matrix(lat_coord, 3), size2 * 4,
                              lat_dim[3])
                .tr_complex();
      }
    }
  }
  return polyakov_loops;
}

template <class T>
std::vector<double>
calculate_polyakov_loops_tr(const std::vector<std::vector<T>> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<double> polyakov_loops(x_size * y_size * z_size);
  link.move_dir(3);
  int place;

  ITER_START_ZYX

  place = link.place / 4;

  polyakov_loops[place] = link.polyakov_loop(array).tr();

  ITER_END_3

  return polyakov_loops;
}

std::vector<std::complex<double>>
calculate_polyakov_loops_tr(const std::vector<std::vector<su3>> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<std::complex<double>> polyakov_loops(x_size * y_size * z_size);
  link.move_dir(3);
  int place;

  ITER_START_ZYX

  place = link.place / 4;

  polyakov_loops[place] = link.polyakov_loop(array).tr_complex();

  ITER_END_3

  return polyakov_loops;
}

std::vector<std::complex<double>> calculate_polyakov_loops_tr(
    const std::vector<std::vector<su3_abelian>> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<std::complex<double>> polyakov_loops(x_size * y_size * z_size);
  link.move_dir(3);
  int place;

  ITER_START_ZYX

  place = link.place / 4;

  polyakov_loops[place] = link.polyakov_loop(array).tr_complex();

  ITER_END_3

  return polyakov_loops;
}

template <class T>
std::vector<double> polyakov_loop_correlator(const std::vector<T> &conf,
                                             int D_max) {
  std::vector<double> polyakov_loops = calculate_polyakov_loops_tr(conf);
  std::vector<int> lat_coord = {0, 0, 0, 0};
  std::vector<int> lat_dim = {x_size, y_size, z_size, t_size};

  int result_size = 2 * D_max + 1;
  int result_size1 = result_size * result_size;

  std::vector<double> correlator(result_size * result_size * result_size);

  double polyakov_tmp;
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
#pragma omp parallel for collapse(2) private(polyakov_tmp, place)              \
    firstprivate(lat_coord, lat_dim, Dx, Dy, Dz) reduction(+ : correlator_tmp)
          for (int z = 0; z < lat_dim[2]; z++) {
            for (int y = 0; y < lat_dim[1]; y++) {
              for (int x = 0; x < lat_dim[0]; x++) {
                lat_coord[0] = x;
                lat_coord[1] = y;
                lat_coord[2] = z;
                polyakov_tmp = polyakov_loops[get_index_site(lat_coord)];
                lat_coord[0] = (lat_coord[0] + lat_dim[0] + Dx) % lat_dim[0];
                lat_coord[1] = (lat_coord[1] + lat_dim[1] + Dy) % lat_dim[1];
                lat_coord[2] = (lat_coord[2] + lat_dim[2] + Dz) % lat_dim[2];
                place = get_index_site(lat_coord);
                correlator_tmp += polyakov_tmp * polyakov_loops[place];
              }
            }
          }
          correlator[correlator_place] = correlator_tmp;
        }
      }
    }
  }
  int size = x_size * y_size * z_size;
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

template <>
std::vector<double> polyakov_loop_correlator(const std::vector<su3> &conf,
                                             int D_max) {
  std::vector<std::complex<double>> polyakov_loops =
      calculate_polyakov_loops_tr(conf);
  std::vector<int> lat_coord = {0, 0, 0, 0};
  std::vector<int> lat_dim = {x_size, y_size, z_size, t_size};

  int result_size = 2 * D_max + 1;
  int result_size1 = result_size * result_size;

  std::vector<double> correlator(result_size * result_size * result_size);

  std::complex<double> polyakov_tmp;
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
#pragma omp parallel for collapse(2) private(polyakov_tmp, place)              \
    firstprivate(lat_coord, lat_dim, Dx, Dy, Dz) reduction(+ : correlator_tmp)
          for (int z = 0; z < lat_dim[2]; z++) {
            for (int y = 0; y < lat_dim[1]; y++) {
              for (int x = 0; x < lat_dim[0]; x++) {
                lat_coord[0] = x;
                lat_coord[1] = y;
                lat_coord[2] = z;
                polyakov_tmp = polyakov_loops[get_index_site(lat_coord)];
                lat_coord[0] = (lat_coord[0] + lat_dim[0] + Dx) % lat_dim[0];
                lat_coord[1] = (lat_coord[1] + lat_dim[1] + Dy) % lat_dim[1];
                lat_coord[2] = (lat_coord[2] + lat_dim[2] + Dz) % lat_dim[2];
                place = get_index_site(lat_coord);
                correlator_tmp +=
                    polyakov_tmp.real() * polyakov_loops[place].real() +
                    polyakov_tmp.imag() * polyakov_loops[place].imag();
              }
            }
          }
          correlator[correlator_place] = correlator_tmp;
        }
      }
    }
  }
  int size = x_size * y_size * z_size;
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

template <>
std::vector<double>
polyakov_loop_correlator(const std::vector<su3_abelian> &conf, int D_max) {
  std::vector<std::complex<double>> polyakov_loops =
      calculate_polyakov_loops_tr(conf);
  std::vector<int> lat_coord = {0, 0, 0, 0};
  std::vector<int> lat_dim = {x_size, y_size, z_size, t_size};

  int result_size = 2 * D_max + 1;
  int result_size1 = result_size * result_size;

  std::vector<double> correlator(result_size * result_size * result_size);

  std::complex<double> polyakov_tmp;
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
#pragma omp parallel for collapse(2) private(polyakov_tmp, place)              \
    firstprivate(lat_coord, lat_dim, Dx, Dy, Dz) reduction(+ : correlator_tmp)
          for (int z = 0; z < lat_dim[2]; z++) {
            for (int y = 0; y < lat_dim[1]; y++) {
              for (int x = 0; x < lat_dim[0]; x++) {
                lat_coord[0] = x;
                lat_coord[1] = y;
                lat_coord[2] = z;
                polyakov_tmp = polyakov_loops[get_index_site(lat_coord)];
                lat_coord[0] = (lat_coord[0] + lat_dim[0] + Dx) % lat_dim[0];
                lat_coord[1] = (lat_coord[1] + lat_dim[1] + Dy) % lat_dim[1];
                lat_coord[2] = (lat_coord[2] + lat_dim[2] + Dz) % lat_dim[2];
                place = get_index_site(lat_coord);
                correlator_tmp +=
                    polyakov_tmp.real() * polyakov_loops[place].real() +
                    polyakov_tmp.imag() * polyakov_loops[place].imag();
              }
            }
          }
          correlator[correlator_place] = correlator_tmp;
        }
      }
    }
  }
  int size = x_size * y_size * z_size;
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

template <class T>
std::vector<double>
polyakov_loop_correlator(const std::vector<std::vector<T>> &conf, int D_max) {
  std::vector<double> polyakov_loops = calculate_polyakov_loops_tr(conf);
  link1 link(x_size, y_size, z_size, t_size);

  int result_size = 2 * D_max + 1;
  int result_size1 = result_size * result_size;

  std::vector<double> correlator(result_size * result_size * result_size);

  double polyakov_tmp;
  double distance;
  int Dx_min;

#pragma omp parallel for collapse(3) private(polyakov_tmp, link, distance,     \
                                                 Dx_min)                       \
    firstprivate(result_size, result_size1)                                    \
    reduction(vec_double_plus : correlator)
  ITER_START_ZYX

  polyakov_tmp = polyakov_loops[link.place / 4];
  link.move(2, -D_max);
  for (int Dz = -D_max; Dz <= D_max; Dz++) {
    link.move(1, -D_max);
    for (int Dy = -D_max; Dy <= D_max; Dy++) {
      if (Dy != 0 && Dz != 0)
        Dx_min = 0;
      else
        Dx_min = -D_max;

      link.move(0, Dx_min);

      for (int Dx = Dx_min; Dx <= D_max; Dx++) {

        distance = sqrt(Dx * Dx + Dy * Dy + Dz * Dz);

        if (distance <= D_max) {

          correlator[(Dz + D_max) * result_size1 + (Dy + D_max) * result_size +
                     Dx + D_max] +=
              polyakov_tmp * polyakov_loops[link.place / 4];
        }

        link.move(0, 1);
      }
      link.move(1, 1);
      link.move(0, -D_max - 1);
    }
    link.move(2, 1);
    link.move(1, -D_max - 1);
  }

  ITER_END_3

  int size = x_size * y_size * z_size;

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

template <>
std::vector<double>
polyakov_loop_correlator(const std::vector<std::vector<su3>> &conf, int D_max) {
  std::vector<std::complex<double>> polyakov_loops =
      calculate_polyakov_loops_tr(conf);
  link1 link(x_size, y_size, z_size, t_size);

  int result_size = 2 * D_max + 1;
  int result_size1 = result_size * result_size;

  std::vector<double> correlator(result_size * result_size * result_size);

  std::complex<double> polyakov_tmp;
  double distance;
  int Dx_min;

#pragma omp parallel for collapse(3) private(polyakov_tmp, link, distance,     \
                                                 Dx_min)                       \
    firstprivate(result_size, result_size1)                                    \
    reduction(vec_double_plus : correlator)
  ITER_START_ZYX

  polyakov_tmp = polyakov_loops[link.place / 4];
  link.move(2, -D_max);
  for (int Dz = -D_max; Dz <= D_max; Dz++) {
    link.move(1, -D_max);
    for (int Dy = -D_max; Dy <= D_max; Dy++) {
      if (Dy != 0 && Dz != 0)
        Dx_min = 0;
      else
        Dx_min = -D_max;

      link.move(0, Dx_min);

      for (int Dx = Dx_min; Dx <= D_max; Dx++) {

        distance = sqrt(Dx * Dx + Dy * Dy + Dz * Dz);

        if (distance <= D_max) {

          correlator[(Dz + D_max) * result_size1 + (Dy + D_max) * result_size +
                     Dx + D_max] +=
              polyakov_tmp.real() * polyakov_loops[link.place / 4].real() +
              polyakov_tmp.imag() * polyakov_loops[link.place / 4].imag();
        }

        link.move(0, 1);
      }
      link.move(1, 1);
      link.move(0, -D_max - 1);
    }
    link.move(2, 1);
    link.move(1, -D_max - 1);
  }

  ITER_END_3

  int size = x_size * y_size * z_size;

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

template <>
std::vector<double>
polyakov_loop_correlator(const std::vector<std::vector<su3_abelian>> &conf,
                         int D_max) {
  std::vector<std::complex<double>> polyakov_loops =
      calculate_polyakov_loops_tr(conf);
  link1 link(x_size, y_size, z_size, t_size);

  int result_size = 2 * D_max + 1;
  int result_size1 = result_size * result_size;

  std::vector<double> correlator(result_size * result_size * result_size);

  std::complex<double> polyakov_tmp;
  double distance;
  int Dx_min;

#pragma omp parallel for collapse(3) private(polyakov_tmp, link, distance,     \
                                                 Dx_min)                       \
    firstprivate(result_size, result_size1)                                    \
    reduction(vec_double_plus : correlator)
  ITER_START_ZYX

  polyakov_tmp = polyakov_loops[link.place / 4];
  link.move(2, -D_max);
  for (int Dz = -D_max; Dz <= D_max; Dz++) {
    link.move(1, -D_max);
    for (int Dy = -D_max; Dy <= D_max; Dy++) {
      if (Dy != 0 && Dz != 0)
        Dx_min = 0;
      else
        Dx_min = -D_max;

      link.move(0, Dx_min);

      for (int Dx = Dx_min; Dx <= D_max; Dx++) {

        distance = sqrt(Dx * Dx + Dy * Dy + Dz * Dz);

        if (distance <= D_max) {

          correlator[(Dz + D_max) * result_size1 + (Dy + D_max) * result_size +
                     Dx + D_max] +=
              polyakov_tmp.real() * polyakov_loops[link.place / 4].real() +
              polyakov_tmp.imag() * polyakov_loops[link.place / 4].imag();
        }

        link.move(0, 1);
      }
      link.move(1, 1);
      link.move(0, -D_max - 1);
    }
    link.move(2, 1);
    link.move(1, -D_max - 1);
  }

  ITER_END_3

  int size = x_size * y_size * z_size;

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

template <class T>
std::vector<T> calculate_polyakov_loops(const std::vector<T> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<T> polyakov_loops(x_size * y_size * z_size);
  link.move_dir(3);
  int place;

  ITER_START_ZYX

  place = link.place / 4;

  polyakov_loops[place] = link.polyakov_loop(array);

  ITER_END_3

  return polyakov_loops;
}

template <class T>
std::vector<T>
calculate_polyakov_loops(const std::vector<std::vector<T>> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<T> polyakov_loops(x_size * y_size * z_size);
  link.move_dir(3);
  int place;

  ITER_START_ZYX

  place = link.place / 4;

  polyakov_loops[place] = link.polyakov_loop(array);

  ITER_END_3

  return polyakov_loops;
}

std::map<double, double>
polyakov_average_directions(const std::vector<double> &correlators, int D_max) {
  std::map<double, double> result;
  std::map<double, int> corr_num;

  double distance;
  int result_size = 2 * D_max + 1;
  int result_size1 = result_size * result_size;
  int Dx_min;

  for (int Dz = -D_max; Dz <= D_max; Dz++) {
    for (int Dy = -D_max; Dy <= D_max; Dy++) {
      if (Dy != 0 && Dz != 0)
        Dx_min = 0;
      else
        Dx_min = -D_max;
      for (int Dx = Dx_min; Dx <= D_max; Dx++) {

        distance = sqrt(Dx * Dx + Dy * Dy + Dz * Dz);

        if (!(Dy != 0 && Dz != 0 && Dx < 0)) {
          if (distance <= D_max && !(Dx == 0 && Dy == 0 && Dz == 0)) {
            result[distance] +=
                correlators[(Dz + D_max) * result_size1 +
                            (Dy + D_max) * result_size + Dx + D_max];
            corr_num[distance]++;
          }
        }
      }
    }
  }

  for (auto it = result.begin(); it != result.end(); ++it) {
    it->second /= corr_num[it->first];
  }

  return result;
}

template <class T>
std::vector<double> polyakov_loop_correlator_singlet(const std::vector<T> &conf,
                                                     int D_max) {
  std::vector<T> polyakov_loops = calculate_polyakov_loops(conf);
  std::map<double, double> polyakov_loop_correlator;
  link1 link(x_size, y_size, z_size, t_size);

  int result_size = 2 * D_max + 1;
  int result_size1 = result_size * result_size;

  std::vector<double> correlator(result_size * result_size * result_size);

  T polyakov_tmp;
  double distance;
  int Dx_min;

#pragma omp parallel for collapse(3) private(polyakov_tmp, link, distance,     \
                                                 Dx_min)                       \
    firstprivate(result_size, result_size1)                                    \
    reduction(vec_double_plus : correlator)
  ITER_START_ZYX

  polyakov_tmp = polyakov_loops[link.place / 4];
  link.move(2, -D_max);
  for (int Dz = -D_max; Dz <= D_max; Dz++) {
    link.move(1, -D_max);
    for (int Dy = -D_max; Dy <= D_max; Dy++) {
      if (Dy != 0 && Dz != 0)
        Dx_min = 0;
      else
        Dx_min = -D_max;

      link.move(0, Dx_min);

      for (int Dx = Dx_min; Dx <= D_max; Dx++) {

        distance = sqrt(Dx * Dx + Dy * Dy + Dz * Dz);

        if (distance <= D_max) {

          correlator[(Dz + D_max) * result_size1 + (Dy + D_max) * result_size +
                     Dx + D_max] +=
              polyakov_tmp.multiply_conj_tr(polyakov_loops[link.place / 4]);
        }

        link.move(0, 1);
      }
      link.move(1, 1);
      link.move(0, -D_max - 1);
    }
    link.move(2, 1);
    link.move(1, -D_max - 1);
  }

  ITER_END_3

  int size = x_size * y_size * z_size;

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

template <class T>
std::vector<double>
polyakov_loop_correlator_singlet(const std::vector<std::vector<T>> &conf,
                                 int D_max) {
  std::vector<T> polyakov_loops = calculate_polyakov_loops(conf);
  std::map<double, double> polyakov_loop_correlator;
  link1 link(x_size, y_size, z_size, t_size);

  int result_size = 2 * D_max + 1;
  int result_size1 = result_size * result_size;

  std::vector<double> correlator(result_size * result_size * result_size);

  T polyakov_tmp;
  double distance;
  int Dx_min;

#pragma omp parallel for collapse(3) private(polyakov_tmp, link, distance,     \
                                                 Dx_min)                       \
    firstprivate(result_size, result_size1)                                    \
    reduction(vec_double_plus : correlator)
  ITER_START_ZYX

  polyakov_tmp = polyakov_loops[link.place / 4];
  link.move(2, -D_max);
  for (int Dz = -D_max; Dz <= D_max; Dz++) {
    link.move(1, -D_max);
    for (int Dy = -D_max; Dy <= D_max; Dy++) {
      if (Dy != 0 && Dz != 0)
        Dx_min = 0;
      else
        Dx_min = -D_max;

      link.move(0, Dx_min);

      for (int Dx = Dx_min; Dx <= D_max; Dx++) {

        distance = sqrt(Dx * Dx + Dy * Dy + Dz * Dz);

        if (distance <= D_max) {

          correlator[(Dz + D_max) * result_size1 + (Dy + D_max) * result_size +
                     Dx + D_max] +=
              polyakov_tmp.multiply_conj_tr(polyakov_loops[link.place / 4]);
        }

        link.move(0, 1);
      }
      link.move(1, 1);
      link.move(0, -D_max - 1);
    }
    link.move(2, 1);
    link.move(1, -D_max - 1);
  }

  ITER_END_3

  int size = x_size * y_size * z_size;

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

// su2
template double polyakov_loop(const std::vector<su2> &array);
template double
polyakov_loop_parallel(const std::vector<std::vector<su2>> &array);
template std::vector<double>
calculate_polyakov_loops_tr(const std::vector<su2> &array);
template std::vector<double>
calculate_polyakov_loops_tr(const std::vector<std::vector<su2>> &array);
template std::vector<double>
polyakov_loop_correlator(const std::vector<su2> &conf, int D_max);
template std::vector<double>
polyakov_loop_correlator(const std::vector<std::vector<su2>> &conf, int D_max);
template std::vector<su2>
calculate_polyakov_loops(const std::vector<su2> &array);
template std::vector<su2>
calculate_polyakov_loops(const std::vector<std::vector<su2>> &array);
template std::vector<double>
polyakov_loop_correlator_singlet(const std::vector<su2> &conf, int D_max);
template std::vector<double>
polyakov_loop_correlator_singlet(const std::vector<std::vector<su2>> &conf,
                                 int D_max);

// abelian
template double polyakov_loop(const std::vector<abelian> &array);
template double
polyakov_loop_parallel(const std::vector<std::vector<abelian>> &array);
template std::vector<double>
calculate_polyakov_loops_tr(const std::vector<abelian> &array);
template std::vector<double>
calculate_polyakov_loops_tr(const std::vector<std::vector<abelian>> &array);
template std::vector<double>
polyakov_loop_correlator(const std::vector<abelian> &conf, int D_max);
template std::vector<double>
polyakov_loop_correlator(const std::vector<std::vector<abelian>> &conf,
                         int D_max);
template std::vector<abelian>
calculate_polyakov_loops(const std::vector<abelian> &array);
template std::vector<abelian>
calculate_polyakov_loops(const std::vector<std::vector<abelian>> &array);
template std::vector<double>
polyakov_loop_correlator_singlet(const std::vector<abelian> &conf, int D_max);
template std::vector<double>
polyakov_loop_correlator_singlet(const std::vector<std::vector<abelian>> &conf,
                                 int D_max);

// su3
template double polyakov_loop(const std::vector<su3> &array);
template double
polyakov_loop_parallel(const std::vector<std::vector<su3>> &array);
template std::vector<su3>
calculate_polyakov_loops(const std::vector<su3> &array);
template std::vector<su3>
calculate_polyakov_loops(const std::vector<std::vector<su3>> &array);
template std::vector<double>
polyakov_loop_correlator_singlet(const std::vector<su3> &conf, int D_max);
template std::vector<double>
polyakov_loop_correlator_singlet(const std::vector<std::vector<su3>> &conf,
                                 int D_max);

// su3_abelian
template double polyakov_loop(const std::vector<su3_abelian> &array);
template double
polyakov_loop_parallel(const std::vector<std::vector<su3_abelian>> &array);
template std::vector<su3_abelian>
calculate_polyakov_loops(const std::vector<su3_abelian> &array);
template std::vector<su3_abelian>
calculate_polyakov_loops(const std::vector<std::vector<su3_abelian>> &array);
template std::vector<double>
polyakov_loop_correlator_singlet(const std::vector<su3_abelian> &conf,
                                 int D_max);
template std::vector<double> polyakov_loop_correlator_singlet(
    const std::vector<std::vector<su3_abelian>> &conf, int D_max);

// su3_angles
template double polyakov_loop(const std::vector<su3_angles> &array);
template double
polyakov_loop_parallel(const std::vector<std::vector<su3_angles>> &array);
template std::vector<su3_angles>
calculate_polyakov_loops(const std::vector<su3_angles> &array);
template std::vector<su3_angles>
calculate_polyakov_loops(const std::vector<std::vector<su3_angles>> &array);
template std::vector<double>
polyakov_loop_correlator_singlet(const std::vector<su3_angles> &conf,
                                 int D_max);
template std::vector<double> polyakov_loop_correlator_singlet(
    const std::vector<std::vector<su3_angles>> &conf, int D_max);