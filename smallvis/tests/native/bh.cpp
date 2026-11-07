// Run from the package root with
// c++ -std=c++17 -pthread tests/native/bh.cpp -o /tmp/smallvis-bh
// /tmp/smallvis-bh
#include "../../src/bh.h"
#include <cassert>
#include <limits>

std::size_t check_counts(const smallvis::QuadTree::Node &node) {
  if (node.terminal_points) {
    std::size_t points = 0;
    for (const auto &point : *node.terminal_points) {
      points += point.count;
    }
    assert(points == node.num_points);
  }
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
    const double eps = std::numeric_limits<double>::epsilon();
    smallvis::QuadTree tree({1 + eps, 0, 1 + 2 * eps, 0, 1 + eps, 0});
    assert(check_counts(*tree.root) == 3);
  }
  // Exercise bucket destruction during exception unwinding as well as success.
  try {
    const double eps = std::numeric_limits<double>::epsilon();
    smallvis::QuadTree tree({1 + eps, 0, 1 + 2 * eps, 0});
    throw std::runtime_error("unwind tree");
  } catch (const std::runtime_error &) {
  }
  {
    const double largest = std::numeric_limits<double>::max();
    smallvis::QuadTree tree({-largest, 0, largest, 0, 0, 0});
    assert(check_counts(*tree.root) == 3);
  }
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
