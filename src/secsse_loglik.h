// Copyright 2023 Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <cassert>
#include <vector>
#include <memory>
#include <string>
#include <utility>
#include <algorithm>
#include "odeint.h"             // NOLINT [build/include_subdir]
#include "DAISIE_DE_rhs.h"         // NOLINT [build/include_subdir]


// retreives value set by RcppParallel::setThreadOptions(numThreads)
// or tbb::task_arena::automatic if missing.
size_t get_rcpp_num_threads();


using state_ptr = std::vector<double>*;

// Models of 'integration_node`
//
//  struct dnode_t {
//    state_ptr state;  // pointer to state
//    double time;      // branch length to ancestor
//    double loglik;    // calculatet loglik
//    ...
// };
//
// struct inode_t {
//   state_ptr state;       // pointer to state
//   dnode_t desc[2];    // descendants
//   double loglik;         // calculated loglik
//   ...
// };

namespace terse {

struct dnode_t {
  state_ptr state = nullptr;
  double time = 0;   // branch length to ancestor
  double loglik = 0.0;
};

struct inode_t {
  state_ptr state = nullptr;
  dnode_t desc[2];
  double loglik = 0.0;
};

}    // namespace terse

namespace storing {

struct storage_t {
  storage_t(double T, const std::vector<double>& State) :
    t(T),
    state(State) {}
  double t;
  std::vector<double> state;
};

struct dnode_t {
  dnode_t() noexcept = default;
  dnode_t(const terse::dnode_t& rhs) noexcept :    // NOLINT [runtime/explicit]
    state(rhs.state),
    time(rhs.time) {}
  state_ptr state = nullptr;
  double time = 0.0;   // branch length to ancestor
  std::vector<storage_t> storage;
};

struct inode_t {
  inode_t() noexcept = default;
  inode_t(const terse::inode_t& rhs) :             // NOLINT [runtime/explicit]
    state(rhs.state),
    desc{rhs.desc[0],
         rhs.desc[1]} {}
  state_ptr state = nullptr;
  dnode_t desc[2];
};

}  // namespace storing

template <typename INODE>
using inodes_t = std::vector<INODE>;


struct phy_edge_t {
  size_t n = 0;
  size_t m = 0;
  double time = 0.0;    // branch length n <-> m
};


// returns phy_edge_t vector sorted by 'N'
inline std::vector<phy_edge_t> make_phy_edge_vector(
    loglik::rmatrix<const double> forTime) {
  auto res = std::vector<phy_edge_t>{forTime.nrow()};
  for (size_t i = 0; i < forTime.nrow(); ++i) {
    auto row = forTime.row(i);
    res[i] = { static_cast<size_t>(row[0]),
               static_cast<size_t>(row[1]),
               row[2] };
  }
  std::sort(std::begin(res), std::end(res), [](auto& a, auto& b) {
    return a.n < b.n;
  });
  return res;
}


inline inodes_t<terse::inode_t> find_inte_nodes(
    const std::vector<phy_edge_t>& phy_edge,
    loglik::rvector<const int> ances,
    std::vector<std::vector<double>>& states) {
  auto res = inodes_t<terse::inode_t>{ances.size()};
  auto comp = [](auto& edge, size_t val) { return edge.n < val; };
  tbb::parallel_for<int>(0, ances.size(), 1, [&](int i) {
    const auto focal = ances[i];
    auto& inode = res[i];
    inode.state = &states[focal - 1];
    inode.state->clear();   // 'dirty' condition
    auto it0 = std::lower_bound(std::begin(phy_edge),
                                std::end(phy_edge),
                                focal, comp);
    auto it1 = std::lower_bound(it0 + 1, std::end(phy_edge), focal, comp);
    // the next thingy is easy to overlook: the sequence matters for creating
    // the 'merged' branch. imposes some pre-condition that is nowere to find :(
    if (it0->m > it1->m) {
      std::swap(it0, it1);
    }
    inode.desc[0] = { &states[it0->m - 1], it0->time };
    inode.desc[1] = { &states[it1->m - 1], it1->time };
  });
  return res;
}


template <typename RaIt>
inline double normalize_loglik(RaIt first, RaIt last) {
  return 0.0;

  const auto sabs = std::accumulate(first, last, 0.0,
                                    [](const auto& s, const auto& x) {
                                      return s + std::abs(x);
                                    });
  if (sabs <= 0.0) return 0.0;
  const auto fact = 1.0 / sabs;
  for (; first != last; ++first) *first *= fact;
  return std::log(sabs);
}

template <typename ODE,
          typename NORMALIZER>
class Integrator {
 public:
  using ode_type = ODE;

  Integrator(std::unique_ptr<ode_type>&& od,
             const std::string& method,
             double atol,
             double rtol) :
    od_(std::move(od)),
    method_(method),
    atol_(atol),
    rtol_(rtol)
  {}

  size_t size() const noexcept { return od_->size(); }

  void operator()(terse::inode_t& inode) const {
    auto s = size();
    std::vector<double> y[2] = { std::vector<double>(s),
                                 std::vector<double>(s) };
#ifdef SECSSE_NESTED_PARALLELISM
    tbb::parallel_for(0, 2, [&](size_t i) {
#else
      for (size_t i = 0; i < 2; ++i) {
#endif
        auto& dnode = inode.desc[i];
        std::copy_n(std::begin(*dnode.state), s, std::begin(y[i]));


        NORMALIZER norm;
        do_integrate(y[i], 0.0, dnode.time, SECSSE_DEFAULT_DTF, norm);
        dnode.loglik = norm.loglik + normalize_loglik(std::begin(y[i]),
                                                      std::end(y[i]));
      }
#ifdef SECSSE_NESTED_PARALLELISM
      );                        // NOLINT [build/include_subdir]
#endif
    inode.state->resize(s);
    od_->mergebranch(y[0], y[1], *inode.state);
    inode.loglik = inode.desc[0].loglik
    + inode.desc[1].loglik
    + normalize_loglik(std::begin(*inode.state),
                       std::end(*inode.state));
    }

    void operator()(std::vector<double>& state, double t0, double t1,
                  NORMALIZER& norm) const {
      do_integrate(state, t0, t1, SECSSE_DEFAULT_DTF, norm);
    }

    void operator()(std::vector<double>& state, double t0, double t1) const {
      odeintcpp::no_normalization no_norm;
      do_integrate(state, t0, t1, SECSSE_DEFAULT_DTF, no_norm);
    }

 private:
    template <typename N>
    void do_integrate(std::vector<double>& state,
                      double t0,
                      double t1,
                      double dtf,
                      N& norm) const {
      odeintcpp::integrate(method_,
                           od_.get(),
                           &state,
                           t0,
                           t1,
                           dtf * (t1 - t0),
                           atol_,
                           rtol_,
                           norm);
    }

    std::unique_ptr<ODE> od_;
    const std::string method_;
    const double atol_;
    const double rtol_;
  };                               // NOLINT [whitespace/indent]


  struct calc_ll_res {
    double loglik;
    std::vector<double> node_M;         // last/root M node
    std::vector<double> merge_branch;   // last/root merged branch
  };

  // generic loglik function
  template <typename INTEGRATOR>
  inline calc_ll_res calc_ll(
      const INTEGRATOR& integrator,
      inodes_t<terse::inode_t>& inodes,
      std::vector<std::vector<double>>& /* in/out */ states) {
    auto is_dirty = [](const auto& inode) {
      return inode.state->empty() &&
        (inode.desc[0].state->empty() || inode.desc[1].state->empty());
    };

    for (auto first = std::begin(inodes); first != std::end(inodes) ;) {
      auto last = std::partition(first,
                                 std::end(inodes),
                                 std::not_fn(is_dirty));
      tbb::parallel_for_each(first, last, [&](auto& inode) {
        integrator(inode);
      });
      first = last;
    }
    // collect output
    const auto& root_node = inodes.back();    // the last calculated
    const auto merge_branch =
      std::vector<double>(std::begin(*root_node.state),
                          std::end(*root_node.state));
    std::vector<double> node_M{ *root_node.desc[1].state };


    integrator(node_M, 0.0, root_node.desc[1].time);

    normalize_loglik(std::begin(node_M), std::end(node_M));
    const auto tot_loglik =
      std::accumulate(std::begin(inodes), std::end(inodes), 0.0,
                      [](auto& sum, const auto& node) {
                        return sum + node.loglik; });
    return { tot_loglik, std::move(node_M), std::move(merge_branch) };
  }


