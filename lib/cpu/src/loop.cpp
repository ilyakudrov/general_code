#include "../include/loop.h"

loop::loop(const link1 &node1) {
  // node = link1(node1);
  coordinate[0] = node1.coordinate[0];
  coordinate[1] = node1.coordinate[1];
  coordinate[2] = node1.coordinate[2];
  coordinate[3] = node1.coordinate[3];
  // for (int i = 0; i <= 6; i++)
  //   link[i] = nullptr;
};

loop::loop(const loop &l) {
  // node = link1(l.node);
  coordinate[0] = l.coordinate[0];
  coordinate[1] = l.coordinate[1];
  coordinate[2] = l.coordinate[2];
  coordinate[3] = l.coordinate[3];
  link = l.link;
  // for (int i = 0; i <= 6; i++)
  //   link[i] = l.link[i];
};