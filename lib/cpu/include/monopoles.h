#ifndef __MONOPOLES_H__
#define __MONOPOLES_H__

#include "link.h"

int test(link1<matrix> l, FLOAT* J);

class loop {
  public:
    link1<matrix> node;
    loop* link[7];

    loop(const link1<matrix> n1);
    loop();
    loop(const loop& l);

//  function calculates where to make the next step
    int get_dir(int i);
};

void find_cluster(loop* ll, FLOAT* J);
void calculate_clusters(vector<loop*>& LL, FLOAT* J);

void calculate_current(data& conf, FLOAT* J);
#endif