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

#include <functional>
#include <thread>

namespace smallvis {

// Trait to detect the arity of the worker (2- or 3- arg)
template <typename T> struct function_traits;

template <typename ReturnType, typename... Args>
struct function_traits<std::function<ReturnType(Args...)>> {
  static constexpr std::size_t arity = sizeof...(Args);
};

template <typename T> auto make_function(T &&t) {
  return std::function{std::forward<T>(t)};
}

// Enable if the callable has 3 arguments (start, end, thread_id)
template <typename WorkerFunc>
std::enable_if_t<function_traits<decltype(make_function(
                     std::declval<WorkerFunc>()))>::arity == 3>
parallel_for(std::size_t N, std::size_t n_threads, WorkerFunc worker) {
  if (n_threads > 1) {
    std::size_t chunk_size = N / n_threads;
    std::vector<std::thread> threads;

    for (std::size_t t = 0; t < n_threads; ++t) {
      std::size_t start = t * chunk_size;
      std::size_t end = (t == n_threads - 1) ? N : (t + 1) * chunk_size;
      threads.emplace_back(worker, start, end, t);
    }

    for (auto &t : threads) {
      t.join();
    }
  } else {
    worker(0, N, 0);
  }
}

// Enable if the callable has 2 arguments (start, end)
template <typename WorkerFunc>
std::enable_if_t<function_traits<decltype(make_function(
                     std::declval<WorkerFunc>()))>::arity == 2>
parallel_for(std::size_t N, std::size_t n_threads, WorkerFunc worker) {
  if (n_threads > 1) {
    std::size_t chunk_size = N / n_threads;
    std::vector<std::thread> threads;

    for (std::size_t t = 0; t < n_threads; ++t) {
      std::size_t start = t * chunk_size;
      std::size_t end = (t == n_threads - 1) ? N : (t + 1) * chunk_size;
      threads.emplace_back(worker, start, end);
    }

    for (auto &t : threads) {
      t.join();
    }
  } else {
    worker(0, N);
  }
}
} // namespace smallvis

#endif // SMALLVIS_THREADS_H