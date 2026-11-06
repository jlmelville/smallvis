// BSD 3-Clause License
//
// Copyright (c) 2024, James Melville
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef SMALLVIS_THREADS_H
#define SMALLVIS_THREADS_H

#include <algorithm>
#include <cstddef>
#include <exception>
#include <thread>
#include <vector>

namespace smallvis {

inline std::size_t worker_count(std::size_t n, std::size_t requested) {
  return std::min(n, std::max(std::size_t{1}, requested));
}

class ThreadJoiner {
public:
  explicit ThreadJoiner(std::vector<std::thread> &threads) : threads(threads) {}
  ThreadJoiner(const ThreadJoiner &) = delete;
  ThreadJoiner &operator=(const ThreadJoiner &) = delete;

  ~ThreadJoiner() {
    for (auto &thread : threads) {
      thread.join();
    }
  }

private:
  std::vector<std::thread> &threads;
};

struct ThreadLauncher {
  template <typename Worker>
  std::thread operator()(Worker worker, std::size_t start, std::size_t end,
                         std::size_t id) const {
    return std::thread(worker, start, end, id);
  }
};

// The launcher is replaceable by standalone native tests of partial startup.
template <typename Worker, typename Launcher = ThreadLauncher>
void parallel_for(std::size_t n, std::size_t requested, Worker worker,
                  Launcher launch = {}) {
  const auto count = worker_count(n, requested);
  if (count == 0) {
    return;
  }
  if (count == 1) {
    worker(0, n, 0);
    return;
  }

  std::vector<std::exception_ptr> failures(count);
  auto guarded_worker = [&](std::size_t start, std::size_t end,
                            std::size_t id) {
    try {
      worker(start, end, id);
    } catch (...) {
      failures[id] = std::current_exception();
    }
  };

  std::vector<std::thread> threads;
  threads.reserve(count);
  {
    ThreadJoiner joiner(threads);
    const auto chunk = n / count;
    const auto remainder = n % count;
    for (std::size_t id = 0; id < count; ++id) {
      const auto start = id * chunk + std::min(id, remainder);
      const auto end = start + chunk + (id < remainder ? 1 : 0);
      threads.push_back(launch(guarded_worker, start, end, id));
    }
  }

  for (const auto &failure : failures) {
    if (failure) {
      std::rethrow_exception(failure);
    }
  }
}
} // namespace smallvis

#endif // SMALLVIS_THREADS_H
