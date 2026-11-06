// Run from the package root with
// c++ -std=c++17 -pthread tests/native/bh.cpp -o /tmp/smallvis-bh
// /tmp/smallvis-bh
#include "../../src/bh.h"
#include <cassert>
#include <limits>

std::size_t check_counts(const smallvis::QuadTree::Node &node) {
  std::size_t children = 0;
  for (const auto &child : node.children) {
    if (child) {
      children += check_counts(*child);
    }
  }
  assert(children == 0 || children == node.num_points);
  return node.num_points;
}

int main() {
  {
    smallvis::QuadTree tree({0, 0, 0, 0, 0, 0, 2, 1, -3, 4, 1, -2});
    assert(check_counts(*tree.root) == 6);
  }
  {
    smallvis::QuadTree tree({-1, 0, 1, 0, 0, 0, 3, 1});
    assert(check_counts(*tree.root) == 4);
  }
  for (double bad : {std::numeric_limits<double>::quiet_NaN(),
                     std::numeric_limits<double>::infinity()}) {
    try {
      smallvis::QuadTree tree({0, 0, bad, 1});
      assert(false);
    } catch (const std::invalid_argument &) {
    }
  }
}
