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
#include "DAISIE_DE_odeint.h"             // NOLINT [build/include_subdir]
#include "DAISIE_DE_rhs.h"                // NOLINT [build/include_subdir]


using state_ptr = std::vector<double>*;

// Models of 'integration_node`
//
//  struct dnode_t {
//    state_ptr state;  // pointer to state
//    double time;      // branch length to ancestor
//    double loglik;    // calculate loglik
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

}  // namespace storing



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

  void operator()(std::vector<double>& state, double t0, double t1,
                NORMALIZER& norm) const {
    do_integrate(state, t0, t1, SECSSE_DEFAULT_DTF, norm);
  }

  void operator()(std::vector<double>& state, double t0, double t1) const {
    odeintcpp::no_normalization no_norm;
    do_integrate(state, t0, t1, SECSSE_DEFAULT_DTF, no_norm);
  }

  void operator()(std::vector<double> init_state,
                  std::vector<double>& times,
                  std::vector< std::vector<double>>* states_out) const {
    odeintcpp::integrate(method_,
                         od_.get(),
                         init_state,
                         times,
                         SECSSE_DEFAULT_DTF,
                         atol_,
                         rtol_,
                         states_out);
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

