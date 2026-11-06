// Run from the package root with
// c++ -std=c++17 -pthread tests/native/threads.cpp -o /tmp/smallvis-threads
// /tmp/smallvis-threads
#include "../../src/threads.h"
#include <atomic>
#include <cassert>
#include <stdexcept>
#include <string>

int main() {
  bool called = false;
  smallvis::parallel_for(
      0, 3, [&](std::size_t, std::size_t, std::size_t) { called = true; });
  assert(!called);

  for (auto requested : {0, 1, 2, 4, 20}) {
    std::vector<int> visits(11, 0);
    std::vector<std::size_t> lengths(smallvis::worker_count(11, requested), 0);
    smallvis::parallel_for(
        11, requested, [&](std::size_t start, std::size_t end, std::size_t id) {
          lengths[id] = end - start;
          for (auto i = start; i < end; ++i) {
            ++visits[i];
          }
        });
    for (auto count : visits) {
      assert(count == 1);
    }
    const auto bounds = std::minmax_element(lengths.begin(), lengths.end());
    assert(*bounds.first > 0 && *bounds.second - *bounds.first <= 1);
  }

  std::atomic<int> completed{0};
  std::size_t launched = 0;
  auto fail_second_launch = [&](auto worker, auto start, auto end, auto id) {
    if (++launched == 2) {
      throw std::runtime_error("startup failure");
    }
    return std::thread(worker, start, end, id);
  };
  try {
    smallvis::parallel_for(
        4, 2, [&](std::size_t, std::size_t, std::size_t) { ++completed; },
        fail_second_launch);
    assert(false);
  } catch (const std::runtime_error &error) {
    assert(std::string(error.what()) == "startup failure");
    assert(completed == 1);
  }

  for (auto requested : {1, 2}) {
    completed = 0;
    try {
      smallvis::parallel_for(4, requested,
                             [&](std::size_t, std::size_t, std::size_t id) {
                               if (id == 0) {
                                 throw std::runtime_error("worker failure");
                               }
                               ++completed;
                             });
      assert(false);
    } catch (const std::runtime_error &error) {
      assert(std::string(error.what()) == "worker failure");
      assert(completed == requested - 1);
    }
  }
}
