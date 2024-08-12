#ifndef SMALLVIS_THREADS_H
#define SMALLVIS_THREADS_H

#include <functional>
#include <thread>

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

#endif // SMALLVIS_THREADS_H