#pragma once

#include "../include/link.h"

#include <array>

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

class loop_new {
public:
  std::array<int, 4> coordinate;
  std::vector<int> charge;
  std::vector<loop_new *> nodes;

  loop_new(const std::array<int, 4> &_node_coordinate)
      : coordinate(_node_coordinate) {}
  // loop_new(const loop_new &_loop) : coordinate(_loop.coordinate),
  // charge(l.charge), nodes(l.nodes) {
  // }
};