#include <vector>

extern int x_size;
extern int y_size;
extern int z_size;
extern int t_size;
extern int size1;
extern int size2;

// int get_index_site(std::vector<int> &lat_coord);
// int get_index_matrix(std::vector<int> &lat_coord, int mu);

inline int get_index_site(std::vector<int> &lat_coord) {
  return size2 * lat_coord[3] + size1 * lat_coord[2] + x_size * lat_coord[1] +
         lat_coord[0];
}

inline int get_index_matrix(std::vector<int> &lat_coord, int mu) {
  return get_index_site(lat_coord) * 4 + mu;
}