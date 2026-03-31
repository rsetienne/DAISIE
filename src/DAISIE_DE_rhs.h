//  Copyright (c) 2021 - 2023, Thijs Janzen
//  Copyright (c) 2023, Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#pragma once
#include <Rcpp.h>
#include <type_traits>
#include <vector>
#include <utility>
#include <string>

namespace loglik {

struct interval {
  const double lc_;   // cladogenesis rate
  const double m_;    // extinction rate
  const double la_;   // anagenesis rate
  const double g_;    // colonisation rate

  // constructor
  interval(const double& lc,
           const double& la,
           const double& m,
           const double& g)
    : lc_(lc),
      m_(m),
      la_(la),
      g_(g){
  }
};

struct interval2_NE : public interval {
  using interval::interval;

  size_t size() const noexcept {
    return 2;
  }

  void operator()(const std::vector<double>& x,
                  std::vector<double>& dxdt,
                  const double /* t */) const {
    auto DM2  = x[0];
    auto E    = x[1];

    // DA1
    dxdt[0] = -(lc_ + m_ + g_ + la_) * DM2;

    // E
    dxdt[1] = m_ - (m_ + lc_) * E + lc_ * E * E;
  }
};

struct interval2_ES : public interval {
  using interval::interval;

  size_t size() const noexcept {
    return 5;
  }

  void operator()(const std::vector<double>& x,
                  std::vector<double>& dxdt,
                  const double /* t */) const {
    auto DE   = x[0];
    auto DM2  = x[1];
    auto DM3  = x[2];
    auto E    = x[3];
    auto DA3  = x[4];

    // DE1
    dxdt[0] = -(lc_ + m_) * DE +
      2 * lc_ * DE * E;
    // DM2
    dxdt[1] =  -(lc_ + m_ + g_ + la_) * DM2 +
      (la_ * DE + 2 * lc_ * DE * E) * DA3;
    // DM3
    dxdt[2] =  -(lc_ + m_ + la_) * DM3 +
      (m_ + la_ * E + lc_ * E * E) * DA3;
    // E
    dxdt[3] =  m_ - (m_ + lc_) * E +
      lc_ * E * E;
    // DA3
    dxdt[4] =  -g_ * DA3 + g_ * DM3;
  }
};


struct interval2_EC : public interval {
  using interval::interval;

  size_t size() const noexcept {
    return 4;
  }

  void operator()(const std::vector<double>& x,
                  std::vector<double>& dxdt,
                  const double /* t */) const {
    auto DE   = x[0];
    auto DM3  = x[1];
    auto E    = x[2];
    auto DA3  = x[3];

    // DE
    dxdt[0] = -(lc_ + m_) * DE +
      2 * lc_ * DE * E;
    // DM3
    dxdt[1] =  -(lc_ + m_ + la_) * DM3 +
      (m_ + la_ * E + lc_ * E * E) * DA3;
    // E
    dxdt[2] =  m_ - (m_ + lc_) * E +
      lc_ * E * E;
    // DA3
    dxdt[3] =  -g_ * DA3 + g_ * DM3;
  }
};




struct interval3_ES : public interval {
  using interval::interval;

  size_t size() const noexcept {
    return 7;
  }

  // this is the dx/dt calculation // true rhs that gets integrated
  // along the branches
  void operator()(const std::vector<double>& x,
                  std::vector<double>& dxdt,
                  const double /* t */) const {
    auto DE  = x[0];
    auto DM1 = x[1];
    auto DM2 = x[2];
    auto DM3 = x[3];
    auto E   = x[4];
    auto DA2 = x[5];
    auto DA3 = x[6];

    // DE
    dxdt[0] = -(lc_ + m_) * DE +
      2 * lc_ * DE * E;
    // DM1
    dxdt[1] = -(lc_ + m_ + la_ + g_) * DM1 +
      g_* DM2 +
      (m_ + la_ * E + lc_ * E * E) * DA2;
    // DM2
    dxdt[2] =  -(lc_ + m_ + la_) * DM2 +
      (m_+ la_ * E + lc_ * E * E) * DA2 +
      (la_ * DE + 2 * lc_ * DE*E) * DA3;
    // DM3
    dxdt[3] =  -(lc_ + m_ + la_) * DM3 +
      (m_ + la_ * E + lc_ * E * E) * DA3;
    // E
    dxdt[4] =  m_ - (m_ + lc_) * E +
      lc_ * E * E;
    // DA2
    dxdt[5] =  -g_ * DA2 + g_ * DM2;

    // DA3
    dxdt[6] =  -g_ * DA3 + g_ * DM3;
  }
};

struct interval3_NE : public interval {
  using interval::interval;

  size_t size() const noexcept {
    // (DE + DM3 + E) * n + DA3
    return 4;
  }
  void operator()(const std::vector<double>& x,
                  std::vector<double>& dxdt,
                  const double /* t */) const {
    auto DM1 = x[0];
    auto DM2 = x[1];
    auto E   = x[2];
    auto DA2 = x[3];

    dxdt[0] = -(lc_ + m_ + la_ + g_) * DM1 +
      (m_ + la_ * E + lc_ * E * E) * DA2 +
      g_ * DM2;

    dxdt[1] = -(lc_ + m_ + la_) * DM2 +
      (m_ + la_ * E + lc_ * E * E) * DA2;


    dxdt[2] = m_ - (m_ + lc_) * E +
      lc_ * E * E;

    dxdt[3] = -g_ * DA2 + g_ * DM2;
  }
};

struct interval4 : public interval {
  using interval::interval;

  size_t size() const noexcept {
    return 3;
  }

  void operator()(const std::vector<double>& x,
                  std::vector<double>& dxdt,
                  const double /* t */) const {
    auto DA1 = x[0];

    auto DM1  = x[1];
    auto E    = x[2];

    // DA1
    dxdt[0] = -g_ * DA1 + g_ * DM1;

    // DM1
    dxdt[1] = -(la_ + m_ + lc_) * DM1 +
      (m_ + la_ * E + lc_ * E * E) * DA1;

    // E
    dxdt[2] = m_ - (m_ + lc_) * E +
      lc_ * E * E;
  }
};

}   // namespace loglik
