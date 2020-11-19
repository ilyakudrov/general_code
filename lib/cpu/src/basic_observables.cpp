#define data_size 4 * x_size *y_size *z_size *t_size
#define PLACE3_DIR                                                             \
  (t) * 3 * x_size *y_size *z_size + (z)*3 * x_size *y_size + (y)*3 * x_size + \
      (x)*3 + dir
#define PLACE3_LINK_NODIR                                                      \
  (link.coordinate[3]) * 3 * x_size *y_size *z_size +                          \
      (link.coordinate[2]) * 3 * x_size *y_size +                              \
      (link.coordinate[1]) * 3 * x_size + (link.coordinate[0]) * 3
#define PLACE1_NODIR                                                           \
  (t) * x_size *y_size *z_size + (z)*x_size *y_size + (y)*x_size + (x)
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
        link.go(x, y, z, t);                                                   \
        link.update(0);                                                        \
        link.update(1);                                                        \
        link.update(2);

#define ITER_END_3                                                             \
  }                                                                            \
  }                                                                            \
  }

#include "../include/basic_observables.h"
#include "../include/flux_tube.h"
#include "../include/link.h"

template <class T> FLOAT plaket_time(const vector<T> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  result res(0);
  FLOAT aver[2];
  for (int dir = 0; dir < 3; dir++) {
    link.move_dir(dir);
    SPACE_ITER_START;
    res.array.push_back(link.plaket_mu(array, 3).tr());
    SPACE_ITER_END;
  }
  res.average(aver);
  return aver[0];
}

template <class T> FLOAT plaket_space(const vector<T> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  result res(0);
  FLOAT aver[2];
  SPACE_ITER_START;
  for (int mu = 0; mu < 3; mu++) {
    for (int nu = mu + 1; nu < 3; nu++) {
      link.move_dir(nu);
      res.array.push_back(link.plaket_mu(array, mu).tr());
    }
  }
  SPACE_ITER_END;
  res.average(aver);
  return aver[0];
}

template <class T> FLOAT plaket(const vector<T> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  result res(data_size / 4 * 6);
  FLOAT aver[2];
  SPACE_ITER_START;
  for (int mu = 0; mu < 4; mu++) {
    for (int nu = mu + 1; nu < 4; nu++) {
      link.move_dir(nu);
      res.array.push_back(link.plaket_mu(array, mu).tr());
    }
  }
  SPACE_ITER_END;
  res.average(aver);
  return aver[0];
}

// fast wilson_loop
// TODO: optimize (in direction 0 wilson lines are calculated way faster, than
// in 3-th)
template <class T>
vector<FLOAT> wilson(const vector<T> &array, int r_min, int r_max, int time_min,
                     int time_max) {
  link1 link(x_size, y_size, z_size, t_size);
  vector<FLOAT> wilson((time_max - time_min + 1) * (r_max - r_min + 1));
  vector<vector<T>> time_lines(time_max - time_min + 1);
  vector<T> space_lines;
  for (int i = time_min; i <= time_max; i++) {
    time_lines[i - time_min] = wilson_lines(array, 3, i);
  }
  T A;
  for (int dir = 0; dir < 3; dir++) {
    for (int r = r_min; r <= r_max; r++) {
      if (r == r_min)
        space_lines = wilson_lines(array, dir, r);
      else
        space_lines = wilson_line_increase(array, space_lines, dir, r - 1);
      for (int time = time_min; time <= time_max; time++) {
        // TODO: elaborate
        for (int t = 0; t < t_size; t++) {
          ITER_START_ZYX
          link.go(x, y, z, t);
          A = time_lines[time - time_min][PLACE1_LINK_NODIR];
          link.move(3, time);
          A = A * space_lines[PLACE1_LINK_NODIR];
          link.move(3, -time);
          link.move(dir, r);
          A = A * time_lines[time - time_min][PLACE1_LINK_NODIR].conj();
          link.move(dir, -r);
          A = A * space_lines[PLACE1_LINK_NODIR].conj();
          wilson[(r - r_min) + (time - time_min) * (r_max - r_min + 1)] +=
              A.tr();
          ITER_END_3
        }
      }
    }
  }
  for (int i = 0; i < (time_max - time_min + 1) * (r_max - r_min + 1); i++) {
    wilson[i] = wilson[i] / (data_size / 4 * 3);
  }
  return wilson;
}

template <class T>
vector<T> wilson_lines(const vector<T> &array, int mu, int length) {
  link1 link(x_size, y_size, z_size, t_size);
  link.move_dir(mu);
  vector<T> vec;
  vec.reserve(data_size / 4);
  SPACE_ITER_START;
  vec[PLACE1_LINK_NODIR] = link.wilson_line(array, length);
  SPACE_ITER_END;
  return vec;
}

// off-axis wilson loop
template <class T>
vector<FLOAT> wilson_offaxis(const vector<T> &array, int r_min, int r_max,
                             int time_min, int time_max, int ratio) {
  link1 link(x_size, y_size, z_size, t_size);
  vector<FLOAT> wilson((time_max - time_min + 1) * (r_max - r_min + 1));
  vector<vector<T>> time_lines(time_max - time_min + 1);
  vector<T> space_lines;
  for (int i = time_min; i <= time_max; i++) {
    time_lines[i - time_min] = wilson_lines(array, 3, i);
  }
  T A;
  for (int dir = 0; dir < 3; dir++) {
    for (int r = r_min; r <= r_max; r++) {
      if (r == r_min)
        space_lines = wilson_lines(array, dir, r);
      else
        space_lines = wilson_line_increase(array, space_lines, dir, r - 1);
      for (int time = time_min; time <= time_max; time++) {
        for (int t = 0; t < t_size; t++) {
          ITER_START_ZYX
          link.go(x, y, z, t);
          A = time_lines[time - time_min][PLACE1_LINK_NODIR];
          link.move(3, time);
          A = A * space_lines[PLACE1_LINK_NODIR];
          link.move(3, -time);
          link.move(dir, r);
          A = A * time_lines[time - time_min][PLACE1_LINK_NODIR].conj();
          link.move(dir, -r);
          A = A * space_lines[PLACE1_LINK_NODIR].conj();
          wilson[(r - r_min) + (time - time_min) * (r_max - r_min + 1)] +=
              A.tr();
          ITER_END_3
        }
      }
    }
  }
  for (int i = 0; i < (time_max - time_min + 1) * (r_max - r_min + 1); i++) {
    wilson[i] = wilson[i] / (data_size / 4 * 3);
  }
  return wilson;
}

template <class T>
vector<T> wilson_line_increase(const vector<T> &array, const vector<T> &lines,
                               int mu, int length) {
  link1 link(x_size, y_size, z_size, t_size);
  link.move_dir(mu);
  vector<T> lines_new(data_size / 4);
  SPACE_ITER_START;
  link.move(mu, length);
  lines_new[PLACE1_NODIR] = lines[PLACE1_NODIR] * link.get_matrix(array);
  SPACE_ITER_END;
  return lines_new;
}

template <class T> FLOAT polyakov(const vector<T> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  result res(0);
  FLOAT aver[2];
  link.move_dir(3);
  SPACE_ITER_START;
  res.array.push_back(link.polyakov_loop(array).tr());
  SPACE_ITER_END;
  res.average(aver);
  return aver[0];
}

template <class T>
FLOAT polyakov_loop_corelator(const vector<T> &array, int D) {
  vector<T> polyakov_loop = calculate_polyakov_loop(array);
  link1 link(x_size, y_size, z_size, t_size);
  result vec(0);
  FLOAT aver[2];
  FLOAT a;
  for (int mu = 0; mu < 3; mu++) {
    SPACE_ITER_START;
    a = polyakov_loop[PLACE1_LINK_NODIR].tr();
    link.move(mu, D);
    a *= polyakov_loop[PLACE1_LINK_NODIR].conj().tr();
    vec.array.push_back(a);
    SPACE_ITER_END;
  }
  vec.average(aver);
  return aver[0];
}

FLOAT MAG_functional_su2(const vector<su2> &array) {
  FLOAT result = 0;
  for (int i = 0; i < data_size; i++) {
    result += (array[i].sigma3_mult() * array[i].conj()).tr();
  }
  return result / (2 * data_size);
}

// monopoles

/*void length(loop* ll, int& ss1){
    if(ll->link[0]==NULL) return ;
    int i=0;
    int Dir=0;
    do {
       length(ll->link[i], ss1);
       Dir=ll->get_dir(i+1);

       if( Dir!=0 ) ss1=ss1+1;

       i++;
    } while((ll->link[i]!=NULL)&&(i<=6));
}

result calculate_cluster_lengths(vector<loop*>& LL, int& max_number){
        int n = 0;
        result res(0);
        int count = 0;
        for(int i = 0;i < LL.size();i++){
                n = 0;
                length(LL[i], n);
                res.array.push_back(n);
                if(n > count) {
                        count = n;
                        max_number = i;
                }
        }
        return res;
}

void length_mu(loop* ll, int mu, int& s){
        if(ll->link[0]==NULL) return ;
        int i=0;
        int dir=0;
        do {
                length_mu(ll->link[i], mu, s);
                dir=ll->get_dir(i+1);
                if(dir == mu) s+=1;
                if(dir == -mu) s-=1;
                i++;
        } while((ll->link[i]!=NULL)&&(i<=6));
}

void calculate_t_clusters(vector<loop*>& LL, vector<loop*>& t_clusters,
int max_number){ int s = 0; for(int i = 0;i < LL.size();i++){ if(i !=
max_number){ s = 0; length_mu(LL[i], 4, s); if(s != 0)
t_clusters.push_back(LL[i]);
                }
        }
}

void calculate_t_clusters_n(vector<loop*>& LL, vector<loop*>&
t_clusters_n, int max_number, int n){ int s = 0; for(int i = 0;i <
LL.size();i++){ if(i != max_number){ s = 0; length_mu(LL[i], 4, s);
if(abs(s/t_size) == n) t_clusters_n.push_back(LL[i]);
                }
        }
}

void calculate_s_clusters(vector<loop*>& LL, vector<loop*>& s_clusters,
int max_number){ int s = 0; for(int i = 0;i < LL.size();i++){ if(i !=
max_number){ for(int j = 1;j < 4;j++){ s = 0; length_mu(LL[i], j, s);
                                if(s != 0) s_clusters.push_back(LL[i]);
                        }
                }
        }
}

FLOAT t_density_n(vector<loop*>& t_clusters, int n){
        int s = 0;
        int count = 0;
        for(int i = 0;i < t_clusters.size();i++){
                s = 0;
                length_mu(t_clusters[i], 4, s);
                if(abs(s/t_size) == n) count++;
        }
        return (FLOAT)count;
}

FLOAT time_length_portion(vector<loop*>& t_clusters){
        result res(0);
        int s1 = 0;
        int s2 = 0;
        for(int i = 0;i < t_clusters.size();i++){
                s1 = 0;
                s2 = 0;
                length(t_clusters[i], s1);
                length_mu(t_clusters[i], 4, s2);
                res.array.push_back(fabs(1.*s1/s2));
        }
        FLOAT aver[2];
        res.average(aver);
        return aver[0];
}

void sites_unique(loop* ll, vector<loop*>& sites){
        int a = 0;
        for(int r=0;r < sites.size();r++){
                if(sites[r]->node.coordinate[0] ==
ll->node.coordinate[0]
                        && sites[r]->node.coordinate[1] ==
ll->node.coordinate[1]
                        && sites[r]->node.coordinate[2] ==
ll->node.coordinate[2]
                        && sites[r]->node.coordinate[3] ==
ll->node.coordinate[3]) a = 1;
        }
        if(a != 1) sites.push_back(ll);
        int i = 0;
        while (ll->link[i]!=NULL && i<=6){
                sites_unique(ll->link[i], sites);
                i++;
        }
}

void aver_r(vector<loop*> sites, FLOAT* aver_coord){
        int size = sites.size();
        aver_coord[0] = 0; aver_coord[1] = 0; aver_coord[2] = 0;
        for(int k = 0;k < size;k++){
                aver_coord[0] += 1.*sites[k]->node.coordinate[0]/size;
                aver_coord[1] += 1.*sites[k]->node.coordinate[1]/size;
                aver_coord[2] += 1.*sites[k]->node.coordinate[2]/size;
        }
}

FLOAT distance_shortest(FLOAT a, FLOAT b){
        if(fabs(a - b) <= (t_size - fabs(a - b))) return fabs(a - b);
        else return (t_size - fabs(a - b));
}

FLOAT disp_r(vector<loop*>& sites, FLOAT* aver_coord){
        FLOAT disp = 0;
        FLOAT dist_x = 0; FLOAT dist_y = 0; FLOAT dist_z = 0;
        for(int k = 0;k < sites.size();k++){
                dist_x = distance_shortest(sites[k]->node.coordinate[0],
aver_coord[0]); dist_y = distance_shortest(sites[k]->node.coordinate[1],
aver_coord[1]); dist_z = distance_shortest(sites[k]->node.coordinate[2],
aver_coord[2]); disp += dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
        }
        disp = disp/sites.size();
        return disp;
}

FLOAT calculate_disp_r(vector<loop*>& t_clusters){
        result res(0);
        vector<loop*> sites(0);
        FLOAT aver_coord[3];
        for(int i = 0;i < t_clusters.size();i++){
                sites_unique(t_clusters[i], sites);
                aver_r(sites, aver_coord);
                res.array.push_back(1./disp_r(sites, aver_coord));
        }
        sites.clear();
        FLOAT aver[2];
        res.average(aver);
        return aver[0];
}

bool sites_close(loop* l, loop* ll){
        int x = distance_shortest(ll->node.coordinate[0],
l->node.coordinate[0]); int y =
distance_shortest(ll->node.coordinate[1], l->node.coordinate[1]); int z
= distance_shortest(ll->node.coordinate[2], l->node.coordinate[2]); int
t = distance_shortest(ll->node.coordinate[3], l->node.coordinate[3]);
return ((x*x + y*y + z*z + t*t) == 1);
}

FLOAT dimension(vector<loop*> sites) {
        int count = 0;
        for(int i = 0;i < sites.size();i++){
                for(int j = 0;j < sites.size();j++){
                        if(sites_close(sites[i], sites[j])) count++;
                }
        }
        return 1.*count/sites.size();
}

FLOAT charge_difference(vector<loop*>& t_clusters_1){
        int count1 = 0;
        int count2 = 0;
        int t_length = 0;
        for(int i = 0;i < t_clusters_1.size();i++){
                t_length = 0;
                length_mu(t_clusters_1[i], 4, t_length);
                if(t_length > 0) count1++;
                if(t_length < 0) count2++;
        }
        return (FLOAT)(count1 - count2);
}*/

// su2
template FLOAT plaket_time(const vector<su2> &array);
template FLOAT plaket_space(const vector<su2> &array);
template FLOAT plaket(const vector<su2> &array);
template vector<FLOAT> wilson(const vector<su2> &array, int r_min, int r_max,
                              int time_min, int time_max);
template vector<su2> wilson_lines(const vector<su2> &array, int mu, int length);
template vector<su2> wilson_line_increase(const vector<su2> &array,
                                          const vector<su2> &lines, int mu,
                                          int length);
template FLOAT polyakov(const vector<su2> &array);
template FLOAT polyakov_loop_corelator(const vector<su2> &array, int D);

// abelian
template FLOAT plaket_time(const vector<abelian> &array);
template FLOAT plaket_space(const vector<abelian> &array);
template FLOAT plaket(const vector<abelian> &array);
template vector<FLOAT> wilson(const vector<abelian> &array, int r_min,
                              int r_max, int time_min, int time_max);
template vector<abelian> wilson_lines(const vector<abelian> &array, int mu,
                                      int length);
template vector<abelian> wilson_line_increase(const vector<abelian> &array,
                                              const vector<abelian> &lines,
                                              int mu, int length);
template FLOAT polyakov(const vector<abelian> &array);
template FLOAT polyakov_loop_corelator(const vector<abelian> &array, int D);