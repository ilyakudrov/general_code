#ifndef __MONOPOLES_H__
#define __MONOPOLES_H__

#include "link.h"

int test(link1 l, double* J);

class loop {
  public:
    link1 node;
    loop* link[7];
  
    loop(const link1 n1);
    loop();
    loop(const loop& l);

//  function calculates where to make the next step
    int get_dir(int i);
};

void find_cluster(loop* ll, double* J);
void calculate_clusters(vector<loop*>& LL, double* J);

void calculate_current(data& conf, double* J);
#endif