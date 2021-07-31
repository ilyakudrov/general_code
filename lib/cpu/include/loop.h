#pragma once

#include "../include/link.h"

class loop {
public:
  // link1 node;
  std::vector<int> coordinate;
  std::vector<int> charge;
  std::vector<loop *> link;

  loop(const link1 &node1);
  loop(const loop &l);

  //  function calculates where to make the next step
  int get_dir(int i);
};