//  Copyright (c) 2021 - 2025, Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#pragma once

// [[Rcpp::depends(BH)]]
#include <utility>   // std::move
#include <memory>    // std::unique_ptr
#include <string>
#include <vector>
#include <type_traits>
#include <algorithm>

#include "config.h"
#include "Rcpp.h"                     // NOLINT [build/include_subdir]
#include "boost/numeric/odeint.hpp"   // NOLINT [build/include_subdir]


#ifdef USE_BULRISCH_STOER_PATCH

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si/dimensionless.hpp>

using bstime_t = boost::units::quantity<boost::units::si::dimensionless,double>;

#else   // USE_BULRISCH_STOER_PATCH

// The default. Causes unitialized member m_last_dt in
// boost::odeint::bulrisch_stoer<>, declared in
// boost/numreic/odeint/stepper/bulrisch_stoer.hpp
using bstime_t = double;

#endif   // USE_BULRISCH_STOER_PATCH

// forward declare
template <typename RaIt>
double normalize_loglik(RaIt first, RaIt last);


namespace odeintcpp {
namespace bno = boost::numeric::odeint;

struct normalize{
  double loglik = 0.0;
};

struct no_normalization{
  double loglik = 0.0; // placeholder
};

template <
  typename STEPPER,
  typename ODE,
  typename STATE,
  typename NORMALIZER
>
void integrate(STEPPER&& stepper, ODE& ode, STATE* y,
               double t0, double t1, double dt,
               NORMALIZER& norm) {


  using time_type = typename STEPPER::time_type;

  if constexpr (std::is_same<NORMALIZER, normalize>::value) {

    auto observer = [&norm](STATE &x, double t) {
      //   auto d = x.size() / 2;
      norm.loglik += 0.0;
    };

    bno::integrate_adaptive(stepper, std::ref(ode), (*y),
                            time_type{t0}, time_type{t1}, time_type{dt},
                            observer);
  } else {
    bno::integrate_adaptive(stepper, std::ref(ode), (*y),
                            time_type{t0}, time_type{t1}, time_type{dt});
  }
}


namespace {

template <typename T>
struct is_unique_ptr : std::false_type {};

                     template <typename T, typename D>
                     struct is_unique_ptr<std::unique_ptr<T, D>> : std::true_type {};

}

template <
  typename STATE,
  typename ODE,
  typename NORMALIZER
>
void integrate(const std::string& stepper_name,
               ODE ode,
               STATE* y,
               double t0,
               double t1,
               double dt,
               double atol,
               double rtol,
               NORMALIZER&  norm) {
  static_assert(is_unique_ptr<ODE>::value ||
                std::is_pointer_v<ODE>,
                "ODE shall be pointer or unique_ptr type");
  if ("odeint::runge_kutta_cash_karp54" == stepper_name) {
    integrate(bno::make_controlled<bno::runge_kutta_cash_karp54<STATE>>(atol,
                                                                        rtol),
                                                                        *ode, y, t0, t1, dt, norm);
  } else if ("odeint::runge_kutta_fehlberg78" == stepper_name) {
    integrate(bno::make_controlled<bno::runge_kutta_fehlberg78<STATE>>(atol,
                                                                       rtol),
                                                                       *ode, y, t0, t1,
                                                                       dt, norm);
  } else if ("odeint::runge_kutta_dopri5" == stepper_name) {
    integrate(bno::make_controlled<bno::runge_kutta_dopri5<STATE>>(atol,
                                                                   rtol),
                                                                   *ode, y, t0, t1,
                                                                   dt, norm);
  } else if ("odeint::bulirsch_stoer" == stepper_name) {
    // no controlled stepper for bulirsch stoer
    integrate(bno::bulirsch_stoer<STATE, double, STATE, bstime_t>(atol,
                                                                  rtol),
                                                                  *ode, y, t0, t1,
                                                                  dt, norm);
  } else if ("odeint::runge_kutta4" == stepper_name) {
    integrate(bno::runge_kutta4<STATE>(), *ode, y, t0, t1, dt, norm);
  } else {
    throw std::runtime_error("odeintcpp::integrate: unknown stepper");
  }
}

}   // namespace odeintcpp
