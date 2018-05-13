
#pragma once

#include <geogram/basic/counted.h>
#include <geogram/basic/memory.h>
#include <geogram/basic/numeric.h>
#include <geogram/basic/thread_sync.h>
#include <algorithm>
#ifdef WIN32
#include <execution>
#endif
#include <numeric>
#include <cassert>

// nTopology replacement for Geogram's process.h

namespace GEO {

namespace Process {
  inline  
  void initialize() {}

  inline
  void show_stats() {}

  inline
  void terminate() { assert(0); std::terminate(); }

  inline
  void brute_force_kill() { std::abort(); }

  inline
  bool is_running_threads() { return false; }

  inline
  index_t maximum_concurrent_threads() { 
    return 8; // ?
  }
};


template <class Functor>
void parallel_for(Functor&& ff,
                  index_t   st,
                  index_t   end,
                  index_t   threads_per_core = 1,
                  bool      interleaved = false)
{
  if (st == end) {
    return;
  }

  if (interleaved) {
    std::cout << "interleaved!" << std::endl;
  }

  if (threads_per_core != 1) {
    std::cout << "threads_per_core != 1!" << std::endl;
  }

  geo_debug_assert(end > st);
  std::vector<int> vv(end - st);
  std::iota(vv.begin(), vv.end(), st);
  std::for_each(
#ifdef WIN32
    std::execution::par_unseq,
#endif
    vv.begin(),
    vv.end(),
    [&](auto&& index) { ff.operator()(index); });
}

}

