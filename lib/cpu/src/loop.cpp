#include "../include/loop.h"

loop::loop(const link1 &node1) {
  std::vector<loop *> link;
  coordinate.reserve(4);
  coordinate[0] = node1.coordinate[0];
  coordinate[1] = node1.coordinate[1];
  coordinate[2] = node1.coordinate[2];
  coordinate[3] = node1.coordinate[3];
};

loop::loop(const loop &l) {
  std::vector<loop *> link(l.link);
  coordinate.reserve(4);
  coordinate[0] = l.coordinate[0];
  coordinate[1] = l.coordinate[1];
  coordinate[2] = l.coordinate[2];
  coordinate[3] = l.coordinate[3];
  link = l.link;
};