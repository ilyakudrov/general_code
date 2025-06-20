#pragma once

#include "../include/matrix.h"

#include <array>
#include <vector>

extern int x_size;
extern int y_size;
extern int z_size;
extern int t_size;
extern int size1;
extern int size2;

// namespace MainLatticeParameters {
// extern int x_size;
// extern int y_size;
// extern int z_size;
// extern int t_size;
// } // namespace MainLatticeParameters

inline int get_index_site(std::array<int, 4> &lat_coord) {
  return size2 * lat_coord[3] + size1 * lat_coord[2] + x_size * lat_coord[1] +
         lat_coord[0];
}

inline int get_index_matrix(std::array<int, 4> &lat_coord, int mu) {
  return get_index_site(lat_coord) * 4 + mu;
}

template <int dim> struct multi_index_t {
  std::array<int, dim> size_array;
  template <typename... Args>
  multi_index_t(Args &&...args) : size_array(std::forward<Args>(args)...) {}
  // multi_index_t(std::array<int, dim> &&array) :
  // size_array(std::forward<std::array<int, dim>>(array)) {}

  struct iterator {
    struct sentinel_t {};

    std::array<int, dim> index_array = {};
    std::array<int, dim> const &size_array;
    bool _end = false;

    iterator(std::array<int, dim> const &size_array) : size_array(size_array) {}

    auto &operator++() {
      for (int i = 0; i < dim; ++i) {
        if (index_array[i] < size_array[i] - 1) {
          ++index_array[i];
          for (int j = 0; j < i; ++j) {
            index_array[j] = 0;
          }
          return *this;
        }
      }
      _end = true;
      return *this;
    }
    auto &operator*() { return index_array; }
    bool operator!=(sentinel_t) const { return !_end; }
  };

  auto begin() const { return iterator{size_array}; }
  auto end() const { return typename iterator::sentinel_t{}; }
};

// class for indexing place of matrix in vector of data in four dimensions
// order of coordinates from the fastest running to the slowest is
// mu, x, y, z, t
class DataPatternLexicographical {
public:
  std::array<int, 4> lat_dim;
  std::array<int, 4> lat_coord;
  int size1;
  int size2;

  DataPatternLexicographical(const std::array<int, 4> &_lat_dim)
      : lat_dim(_lat_dim) {
    lat_coord = {0, 0, 0, 0};
    size1 = lat_dim[0] * lat_dim[1];
    size2 = lat_dim[0] * lat_dim[1] * lat_dim[2];
  }

  DataPatternLexicographical(const DataPatternLexicographical &_data_pattern)
      : lat_dim(_data_pattern.lat_dim), lat_coord(_data_pattern.lat_coord),
        size1(_data_pattern.size1), size2(_data_pattern.size2) {}

  int get_data_size() const {
    return lat_dim[0] * lat_dim[1] * lat_dim[2] * lat_dim[3] * 4;
  }

  int get_lattice_size() const {
    return lat_dim[0] * lat_dim[1] * lat_dim[2] * lat_dim[3];
  }

  void move_forward(int length, int mu) {
    lat_coord[mu] = (lat_coord[mu] + length) % lat_dim[mu];
  }

  void move_backward(int length, int mu) {
    lat_coord[mu] = (lat_coord[mu] - length + lat_dim[mu]) % lat_dim[mu];
  }

  // move when length can be negative
  inline void move(int length, int mu) {
    lat_coord[mu] = (lat_coord[mu] + length + lat_dim[mu]) % lat_dim[mu];
  }

  // place of site if link direction is neglected, indices run as x, y, z ,t
  int get_index_site() const {
    return size2 * lat_coord[3] + size1 * lat_coord[2] +
           lat_dim[0] * lat_coord[1] + lat_coord[0];
  }

  // place of site if link direction is neglected and there were only x, y, z
  // coordinates
  int get_index_site_spacial() const {
    return size1 * lat_coord[2] + lat_dim[0] * lat_coord[1] + lat_coord[0];
  }

  // place of link in vector of data
  int get_index_link(int mu) const { return get_index_site() * 4 + mu; }

  multi_index_t<4> get_multi_index() { return multi_index_t<4>(lat_dim); }
};

// class for indexing place of matrix in vector of data in four dimensions
// order of coordinates from the fastest running to the slowest is
// x, y, z, t, mu
class DataPatternSeparateDir {
public:
  std::array<int, 4> lat_dim;
  std::array<int, 4> lat_coord;
  int size1;
  int size2;
  int size3;

  DataPatternSeparateDir(const std::array<int, 4> &_lat_dim)
      : lat_dim(_lat_dim) {
    lat_coord = {0, 0, 0, 0};
    size1 = lat_dim[0] * lat_dim[1];
    size2 = lat_dim[0] * lat_dim[1] * lat_dim[2];
    size3 = lat_dim[0] * lat_dim[1] * lat_dim[2] * lat_dim[3];
  }

  DataPatternSeparateDir(const DataPatternSeparateDir &_data_pattern)
      : lat_dim(_data_pattern.lat_dim), lat_coord(_data_pattern.lat_coord),
        size1(_data_pattern.size1), size2(_data_pattern.size2) {}

  int get_data_size() const {
    return lat_dim[0] * lat_dim[1] * lat_dim[2] * lat_dim[3] * 4;
  }

  int get_lattice_size() const {
    return lat_dim[0] * lat_dim[1] * lat_dim[2] * lat_dim[3];
  }

  inline void move_forward(int length, int mu) {
    lat_coord[mu] = (lat_coord[mu] + length) % lat_dim[mu];
  }

  inline void move_backward(int length, int mu) {
    lat_coord[mu] = (lat_coord[mu] - length + lat_dim[mu]) % lat_dim[mu];
  }

  // move when length can be negative
  inline void move(int length, int mu) {
    lat_coord[mu] = (lat_coord[mu] + length + lat_dim[mu]) % lat_dim[mu];
  }

  // place of site if link direction is neglected, indices run as x, y, z, t
  int get_index_site() const {
    return size2 * lat_coord[3] + size1 * lat_coord[2] +
           lat_dim[0] * lat_coord[1] + lat_coord[0];
  }

  // place of site if link direction is neglected and there were only x, y, z
  // coordinates
  int get_index_site_spacial() const {
    return size1 * lat_coord[2] + lat_dim[0] * lat_coord[1] + lat_coord[0];
  }

  // place of link in vector of data
  int get_index_link(int mu) const { return get_index_site() + mu * size3; }

  multi_index_t<4> get_multi_index() { return multi_index_t<4>(lat_dim); }
};

template <int N_dim, class MatrixType> class FilePatternLexicographical {
public:
  typedef MatrixType matrix_type;
  FilePatternLexicographical() {}
  int get_index_site(std::array<int, N_dim> &lat_dim,
                     std::array<int, N_dim> &lat_coord) const {
    int index = 0;
    int size = 1;
    for (int i = 0; i < N_dim; i++) {
      index += lat_coord[i] * size;
      size *= lat_dim[i];
    }
    return index;
  }

  int get_index_link(std::array<int, N_dim> &lat_dim,
                     std::array<int, N_dim> &lat_coord, int mu) const {
    return get_index_site(lat_dim, lat_coord) * N_dim + mu;
  }

  int get_index_matrix_data(std::array<int, N_dim> &lat_dim,
                            std::array<int, N_dim> &lat_coord, int mu,
                            int element_num, int matrix_data_size) const {
    return get_index_link(lat_dim, lat_coord, mu) * matrix_data_size +
           element_num;
  }
};

template <int N_dim> class FilePatternLexicographical<N_dim, su3_abelian> {
public:
  typedef su3_abelian matrix_type;
  FilePatternLexicographical() {}
  int get_index_site(std::array<int, N_dim> &lat_dim,
                     std::array<int, N_dim> &lat_coord) const {
    int index = 0;
    int size = 1;
    for (int i = 0; i < N_dim; i++) {
      index += lat_coord[i] * size;
      size *= lat_dim[i];
    }
    return index;
  }

  int get_index_link(std::array<int, N_dim> &lat_dim,
                     std::array<int, N_dim> &lat_coord, int mu) const {
    return get_index_site(lat_dim, lat_coord) * N_dim + mu;
  }

  int get_index_matrix_data(std::array<int, N_dim> &lat_dim,
                            std::array<int, N_dim> &lat_coord, int mu,
                            int element_num, int matrix_data_size) const {
    int size = 1;
    for (int i = 0; i < N_dim; i++) {
      size *= lat_dim[i];
    }
    return get_index_link(lat_dim, lat_coord, mu) + element_num * size * N_dim;
  }
};

template <int N_dim, class MatrixType> class FilePatternQCDSTAG {
public:
  typedef MatrixType matrix_type;
  FilePatternQCDSTAG() {}
  int get_index_site(std::array<int, N_dim> &lat_dim,
                     std::array<int, N_dim> &lat_coord) const {
    int index = 0;
    int size = 1;
    for (int i = 0; i < N_dim; i++) {
      index += lat_coord[i] * size;
      size *= lat_dim[i];
    }
    return index;
  }

  int get_index_link(std::array<int, N_dim> &lat_dim,
                     std::array<int, N_dim> &lat_coord, int mu) const {
    int mu1;
    if (mu == lat_dim.size() - 1)
      mu1 = 0;
    else
      mu1 = mu + 1;
    int size = 1;
    for (int i = 0; i < N_dim; i++) {
      size *= lat_dim[i];
    }
    return get_index_site(lat_dim, lat_coord) + mu1 * size;
  }

  int get_index_matrix_data(std::array<int, N_dim> &lat_dim,
                            std::array<int, N_dim> &lat_coord, int mu,
                            int element_num, int matrix_data_size) const {
    return get_index_link(lat_dim, lat_coord, mu) * matrix_data_size +
           element_num;
  }
};

template <int N_dim> class FilePatternQCDSTAG<N_dim, su2> {
public:
  typedef su2 matrix_type;
  FilePatternQCDSTAG() {}
  int get_index_site(std::array<int, N_dim> &lat_dim,
                     std::array<int, N_dim> &lat_coord) const {
    int index = 0;
    int size = 1;
    for (int i = 0; i < N_dim; i++) {
      index += lat_coord[i] * size;
      size *= lat_dim[i];
    }
    return index;
  }

  int get_index_link(std::array<int, N_dim> &lat_dim,
                     std::array<int, N_dim> &lat_coord, int mu) const {
    int mu1;
    if (mu == lat_dim.size() - 1)
      mu1 = 0;
    else
      mu1 = mu + 1;
    int size = 1;
    for (int i = 0; i < N_dim; i++) {
      size *= lat_dim[i];
    }
    return get_index_site(lat_dim, lat_coord) + mu1 * size;
  }

  int get_index_matrix_data(std::array<int, N_dim> &lat_dim,
                            std::array<int, N_dim> &lat_coord, int mu,
                            int element_num, int matrix_data_size) const {
    if (element_num > 0)
      element_num = N_dim - element_num;
    return get_index_link(lat_dim, lat_coord, mu) * matrix_data_size +
           element_num;
  }
};

template <int N_dim, class MatrixType> class FilePatternILDG {
public:
  typedef MatrixType matrix_type;
  FilePatternILDG() {}
  int get_index_site(std::array<int, N_dim> &lat_dim,
                     std::array<int, N_dim> &lat_coord) const {
    int index = 0;
    int size = 1;
    for (int i = 0; i < N_dim; i++) {
      index += lat_coord[i] * size;
      size *= lat_dim[i];
    }
    return index;
  }

  int get_index_link(std::array<int, N_dim> &lat_dim,
                     std::array<int, N_dim> &lat_coord, int mu) const {
    return get_index_site(lat_dim, lat_coord) * N_dim + mu;
  }

  int get_index_matrix_data(std::array<int, N_dim> &lat_dim,
                            std::array<int, N_dim> &lat_coord, int mu,
                            int element_num, int matrix_data_size) const {
    return get_index_link(lat_dim, lat_coord, mu) * matrix_data_size +
           element_num;
  }
};

// template <typename... index_t> auto multi_index(index_t &&...index) {
//   static constexpr int size = sizeof...(index_t);
//   auto ar = std::array<int, size>{std::forward<index_t>(index)...};
//   return multi_index_t<size>(ar);
// }