#include "../include/monopoles.h"

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

// read configuration of angles
vector<FLOAT> read_angles_float_fortran(string &file_path) {
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  vector<FLOAT> angles;
  angles.reserve(data_size1);
  ifstream stream(file_path);
  vector<float> v(data_size1 + 2);
  if (!stream.read((char *)&v[0], (data_size1 + 2) * sizeof(float)))
    cout << "read_angles_float_fortran error: " << file_path << endl;
  for (int i = 0; i < data_size1; i++) {
    angles.push_back((FLOAT)v[i + 1]);
  }
  stream.close();
  return angles;
}

// calculate monopole_plaketes on the lattice
vector<vector<FLOAT>> calculate_monopole_plaket(vector<FLOAT> &angles) {
  int data_size = x_size * y_size * z_size * t_size;
  vector<vector<FLOAT>> plakets(6, vector<FLOAT>(data_size));
  link1 link(x_size, y_size, z_size, t_size);

  int count = 0;
  for (int mu = 0; mu < 4; mu++) {
    link.move_dir(mu);
    for (int nu = mu + 1; nu < 4; nu++) {

      SPACE_ITER_START

      plakets[count][link.place / 4] = link.monopole_plaket_mu(angles, nu);

      SPACE_ITER_END

      count++;
    }
  }
  return plakets;
}

// int loop::get_dir(int i) {
//   loop *tmp = link[i - 1];
//   if (tmp == NULL)
//     return 0;

//   int delta = (tmp->node.coordinate[0]) - node.coordinate[0];
//   if (delta == (1 - x_size))
//     return 1;
//   if (delta == (x_size - 1))
//     return -1;
//   if (delta != 0)
//     return 1 * delta;

//   delta = (tmp->node.coordinate[1]) - node.coordinate[1];
//   if (delta == (1 - y_size))
//     return 2;
//   if (delta == (y_size - 1))
//     return -2;
//   if (delta != 0)
//     return delta * 2;

//   delta = (tmp->node.coordinate[2]) - node.coordinate[2];
//   if (delta == (1 - z_size))
//     return 3;
//   if (delta == (z_size - 1))
//     return -3;
//   if (delta != 0)
//     return delta * 3;

//   delta = (tmp->node.coordinate[3]) - node.coordinate[3];
//   if (delta == (1 - t_size))
//     return 4;
//   if (delta == (t_size - 1))
//     return -4;
//   if (delta != 0)
//     return delta * 4;
//   return 0;
// }

vector<FLOAT> calculate_current(vector<FLOAT> &angles) {
  int data_size = 4 * x_size * y_size * z_size * t_size;
  link1 link(x_size, y_size, z_size, t_size);

  vector<FLOAT> J(data_size);

  vector<vector<FLOAT>> monopole_plaket = calculate_monopole_plaket(angles);

  SPACE_ITER_START

  link.get_current(monopole_plaket, &J[link.place]);

  SPACE_ITER_END

  return J;
}

int find_current(link1 &link, vector<FLOAT> &J) {
  for (int mu = 0; mu < 4; mu++) {
    if ((J[link.place + mu] > 0.3) || (J[link.place + mu] < -0.3))
      return mu + 1;
  }
  for (int mu = 0; mu < 4; mu++) {
    link.move(mu, -1);
    if ((J[link.place + mu] > 0.3) || (J[link.place + mu] < -0.3))
      return -mu - 1;
    link.move(mu, 1);
  }
  return 0;
}

void find_cluster(link1 &link, vector<FLOAT> &J, loop &ll, int dir_back) {
  int dir;
  int dir1;
  do {
    dir = find_current(link, J);
    if (dir > 0) {
      J[link.place + dir - 1] = 0.;
      link.move(dir - 1, 1);
      ll.link.push_back(new loop(link));
    }
    if (dir < 0) {
      link.move(-dir - 1, -1);
      J[link.place - dir - 1] = 0.;
      ll.link.push_back(new loop(link));
    }

    dir1 = find_current(link, J);
    // cout << dir1 << endl;
    if (dir1 != 0)
      find_cluster(link, J, *ll.link[ll.link.size() - 1], dir);
  } while (find_current(link, J) != 0);
  if (dir_back > 0)
    link.move(dir_back - 1, -1);
  else
    link.move(-dir_back - 1, 1);
};

void calculate_clusters(vector<loop *> &LL, vector<FLOAT> &J) {
  int dir1;

  link1 link(x_size, y_size, z_size, t_size);

  SPACE_ITER_START

  // if (x == 0 && y == 0 && z == 0 && t == 0) {
  //   cout << link.place << endl;
  //   cout << J[0] << endl;
  //   cout << J[link.place] << endl;
  // }
  dir1 = find_current(link, J);
  if (dir1 != 0) {
    LL.push_back(new loop(link));
    find_cluster(link, J, *LL[LL.size() - 1], dir1);
  }

  SPACE_ITER_END
}

// monopole observables

void cluster_length(loop &ll, int &length) {
  for (int i = 0; i < ll.link.size(); i++) {
    length++;
    cluster_length(*ll.link[i], length);
  }
}

// result calculate_cluster_lengths(vector<loop *> &LL, int &max_number) {
//   int n = 0;
//   result res(0);
//   int count = 0;
//   for (int i = 0; i < LL.size(); i++) {
//     n = 0;
//     length(LL[i], n);
//     res.array.push_back(n);
//     if (n > count) {
//       count = n;
//       max_number = i;
//     }
//   }
//   return res;
// }

// void length_mu(loop *ll, int mu, int &s) {
//   if (ll->link[0] == NULL)
//     return;
//   int i = 0;
//   int dir = 0;
//   do {
//     length_mu(ll->link[i], mu, s);
//     dir = ll->get_dir(i + 1);
//     if (dir == mu)
//       s += 1;
//     if (dir == -mu)
//       s -= 1;
//     i++;
//   } while ((ll->link[i] != NULL) && (i <= 6));
// }

// void calculate_t_clusters(vector<loop *> &LL, vector<loop *> &t_clusters,
//                           int max_number) {
//   int s = 0;
//   for (int i = 0; i < LL.size(); i++) {
//     if (i != max_number) {
//       s = 0;
//       length_mu(LL[i], 4, s);
//       if (s != 0)
//         t_clusters.push_back(LL[i]);
//     }
//   }
// }

// void calculate_t_clusters_n(vector<loop *> &LL, vector<loop *> &t_clusters_n,
//                             int max_number, int n) {
//   int s = 0;
//   for (int i = 0; i < LL.size(); i++) {
//     if (i != max_number) {
//       s = 0;
//       length_mu(LL[i], 4, s);
//       if (abs(s / t_size) == n)
//         t_clusters_n.push_back(LL[i]);
//     }
//   }
// }

// void calculate_s_clusters(vector<loop *> &LL, vector<loop *> &s_clusters,
//                           int max_number) {
//   int s = 0;
//   for (int i = 0; i < LL.size(); i++) {
//     if (i != max_number) {
//       for (int j = 1; j < 4; j++) {
//         s = 0;
//         length_mu(LL[i], j, s);
//         if (s != 0)
//           s_clusters.push_back(LL[i]);
//       }
//     }
//   }
// }

// FLOAT t_density_n(vector<loop *> &t_clusters, int n) {
//   int s = 0;
//   int count = 0;
//   for (int i = 0; i < t_clusters.size(); i++) {
//     s = 0;
//     length_mu(t_clusters[i], 4, s);
//     if (abs(s / t_size) == n)
//       count++;
//   }
//   return (FLOAT)count;
// }

// FLOAT time_length_portion(vector<loop *> &t_clusters) {
//   result res(0);
//   int s1 = 0;
//   int s2 = 0;
//   for (int i = 0; i < t_clusters.size(); i++) {
//     s1 = 0;
//     s2 = 0;
//     length(t_clusters[i], s1);
//     length_mu(t_clusters[i], 4, s2);
//     res.array.push_back(fabs(1. * s1 / s2));
//   }
//   FLOAT aver[2];
//   res.average(aver);
//   return aver[0];
// }

// void sites_unique(loop *ll, vector<loop *> &sites) {
//   int a = 0;
//   for (int r = 0; r < sites.size(); r++) {
//     if (sites[r]->node.coordinate[0] == ll->node.coordinate[0] &&
//         sites[r]->node.coordinate[1] == ll->node.coordinate[1] &&
//         sites[r]->node.coordinate[2] == ll->node.coordinate[2] &&
//         sites[r]->node.coordinate[3] == ll->node.coordinate[3])
//       a = 1;
//   }
//   if (a != 1)
//     sites.push_back(ll);
//   int i = 0;
//   while (ll->link[i] != NULL && i <= 6) {
//     sites_unique(ll->link[i], sites);
//     i++;
//   }
// }

// void aver_r(vector<loop *> sites, FLOAT *aver_coord) {
//   int size = sites.size();
//   aver_coord[0] = 0;
//   aver_coord[1] = 0;
//   aver_coord[2] = 0;
//   for (int k = 0; k < size; k++) {
//     aver_coord[0] += 1. * sites[k]->node.coordinate[0] / size;
//     aver_coord[1] += 1. * sites[k]->node.coordinate[1] / size;
//     aver_coord[2] += 1. * sites[k]->node.coordinate[2] / size;
//   }
// }

// FLOAT distance_shortest(FLOAT a, FLOAT b) {
//   if (fabs(a - b) <= (t_size - fabs(a - b)))
//     return fabs(a - b);
//   else
//     return (t_size - fabs(a - b));
// }

// FLOAT disp_r(vector<loop *> &sites, FLOAT *aver_coord) {
//   FLOAT disp = 0;
//   FLOAT dist_x = 0;
//   FLOAT dist_y = 0;
//   FLOAT dist_z = 0;
//   for (int k = 0; k < sites.size(); k++) {
//     dist_x = distance_shortest(sites[k]->node.coordinate[0], aver_coord[0]);
//     dist_y = distance_shortest(sites[k]->node.coordinate[1], aver_coord[1]);
//     dist_z = distance_shortest(sites[k]->node.coordinate[2], aver_coord[2]);
//     disp += dist_x * dist_x + dist_y * dist_y + dist_z * dist_z;
//   }
//   disp = disp / sites.size();
//   return disp;
// }

// FLOAT calculate_disp_r(vector<loop *> &t_clusters) {
//   result res(0);
//   vector<loop *> sites(0);
//   FLOAT aver_coord[3];
//   for (int i = 0; i < t_clusters.size(); i++) {
//     sites_unique(t_clusters[i], sites);
//     aver_r(sites, aver_coord);
//     res.array.push_back(1. / disp_r(sites, aver_coord));
//   }
//   sites.clear();
//   FLOAT aver[2];
//   res.average(aver);
//   return aver[0];
// }

// bool sites_close(loop *l, loop *ll) {
//   int x = distance_shortest(ll->node.coordinate[0], l->node.coordinate[0]);
//   int y = distance_shortest(ll->node.coordinate[1], l->node.coordinate[1]);
//   int z = distance_shortest(ll->node.coordinate[2], l->node.coordinate[2]);
//   int t = distance_shortest(ll->node.coordinate[3], l->node.coordinate[3]);
//   return ((x * x + y * y + z * z + t * t) == 1);
// }

// FLOAT dimension(vector<loop *> sites) {
//   int count = 0;
//   for (int i = 0; i < sites.size(); i++) {
//     for (int j = 0; j < sites.size(); j++) {
//       if (sites_close(sites[i], sites[j]))
//         count++;
//     }
//   }
//   return 1. * count / sites.size();
// }

// FLOAT charge_difference(vector<loop *> &t_clusters_1) {
//   int count1 = 0;
//   int count2 = 0;
//   int t_length = 0;
//   for (int i = 0; i < t_clusters_1.size(); i++) {
//     t_length = 0;
//     length_mu(t_clusters_1[i], 4, t_length);
//     if (t_length > 0)
//       count1++;
//     if (t_length < 0)
//       count2++;
//   }
//   return (FLOAT)(count1 - count2);