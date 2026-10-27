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
#include <complex>

namespace DAISIEDE {

  using complex = std::complex<double>;

struct interval {
  const double lc_;   // cladogenesis rate
  const double mu_E_;    // extinction rate
  const double mu_NE_;
  const double la_;   // anagenesis rate
  const double g_;    // colonisation rate

  // constructor
  interval(double lc,
           double la,
           double mu_E,
           double mu_NE,
           double g) :
      lc_(lc),
      mu_E_(mu_E),
      mu_NE_(mu_NE),
      la_(la),
      g_(g)
    {}
};

template < typename datatype >
struct interval2_NE : public interval {
  using interval::interval;

  size_t size() const noexcept {
    return 2;
  }

  void operator()(const std::vector<datatype>& x,
                        std::vector<datatype>& dxdt,
                        const double /* t */) const {
    auto DM2  = x[0];
    auto E    = x[1];

    // DA1
    dxdt[0] = -(lc_ + mu_NE_ + g_ + la_) * DM2;

    // E
    dxdt[1] = mu_E_ - (mu_E_ + lc_) * E + lc_ * E * E;
  }
};

template <typename datatype>
struct interval2_ES : public interval {
  using interval::interval;

  size_t size() const noexcept {
    return 5;
  }

  void operator()(const std::vector<datatype>& x,
                  std::vector<datatype>& dxdt,
                  const double /* t */) const {
    auto DE   = x[0];
    auto DM2  = x[1];
    auto DM3  = x[2];
    auto E    = x[3];
    auto DA3  = x[4];

    // DE1
    dxdt[0] = -(lc_ + mu_E_) * DE +
      2.0 * lc_ * DE * E;
    // DM2
    dxdt[1] =  -(lc_ + mu_NE_ + g_ + la_) * DM2 +
      (la_ * DE + 2.0 * lc_ * DE * E) * DA3;
    // DM3
    dxdt[2] =  -(lc_ + mu_NE_ + la_) * DM3 +
      (mu_NE_ + la_ * E + lc_ * E * E) * DA3;
    // E
    dxdt[3] =  mu_E_ - (mu_E_ + lc_) * E +
      lc_ * E * E;
    // DA3
    dxdt[4] =  -g_ * DA3 + g_ * DM3;
  }
};

template <typename datatype>
struct interval2_EC : public interval {
  using interval::interval;

  size_t size() const noexcept {
    return 4;
  }

  void operator()(const std::vector< datatype >& x,
                  std::vector< datatype >& dxdt,
                  const double /* t */) const {
    auto DE   = x[0];
    auto DM3  = x[1];
    auto E    = x[2];
    auto DA3  = x[3];

    // DE
    dxdt[0] = -(lc_ + mu_E_) * DE +
      2.0 * lc_ * DE * E;
    // DM3
    dxdt[1] =  -(lc_ + mu_NE_ + la_) * DM3 +
      (mu_NE_ + la_ * E + lc_ * E * E) * DA3;
    // E
    dxdt[2] =  mu_E_ - (mu_E_ + lc_) * E +
      lc_ * E * E;
    // DA3
    dxdt[3] =  -g_ * DA3 + g_ * DM3;
  }
};



template <typename datatype>
struct interval3_ES : public interval {
  using interval::interval;

  size_t size() const noexcept {
    return 7;
  }

  // this is the dx/dt calculation // true rhs that gets integrated
  // along the branches
  void operator()(const std::vector<datatype>& x,
                  std::vector<datatype>& dxdt,
                  const double /* t */) const {
    auto DE  = x[0];
    auto DM1 = x[1];
    auto DM2 = x[2];
    auto DM3 = x[3];
    auto E   = x[4];
    auto DA2 = x[5];
    auto DA3 = x[6];

    // DE
    dxdt[0] = -(lc_ + mu_E_) * DE +
      2.0 * lc_ * DE * E;
    // DM1
    dxdt[1] = -(lc_ + mu_NE_ + la_ + g_) * DM1 +
      g_* DM2 +
      (mu_NE_ + la_ * E + lc_ * E * E) * DA2;
    // DM2
    dxdt[2] =  -(lc_ + mu_NE_ + la_) * DM2 +
      (mu_NE_+ la_ * E + lc_ * E * E) * DA2 +
      (la_ * DE + 2.0 * lc_ * DE*E) * DA3;
    // DM3
    dxdt[3] =  -(lc_ + mu_NE_ + la_) * DM3 +
      (mu_NE_ + la_ * E + lc_ * E * E) * DA3;
    // E
    dxdt[4] =  mu_E_ - (mu_E_ + lc_) * E +
      lc_ * E * E;
    // DA2
    dxdt[5] =  -g_ * DA2 + g_ * DM2;

    // DA3
    dxdt[6] =  -g_ * DA3 + g_ * DM3;
  }
};

template <typename datatype>
struct interval3_NE : public interval {
  using interval::interval;

  size_t size() const noexcept {
    // (DE + DM3 + E) * n + DA3
    return 4;
  }
  void operator()(const std::vector<datatype>& x,
                  std::vector<datatype>& dxdt,
                  const double /* t */) const {
    auto DM1 = x[0];
    auto DM2 = x[1];
    auto E   = x[2];
    auto DA2 = x[3];

    dxdt[0] = -(lc_ + mu_NE_ + la_ + g_) * DM1 +
      (mu_NE_ + la_ * E + lc_ * E * E) * DA2 +
      g_ * DM2;

    dxdt[1] = -(lc_ + mu_NE_ + la_) * DM2 +
      (mu_NE_ + la_ * E + lc_ * E * E) * DA2;

    dxdt[2] = mu_E_ - (mu_E_ + lc_) * E +
      lc_ * E * E;

    dxdt[3] = -g_ * DA2 + g_ * DM2;
  }
};

template <typename datatype>
struct interval4 : public interval {
  using interval::interval;

  size_t size() const noexcept {
    return 3;
  }

  void operator()(const std::vector<datatype>& x,
                        std::vector<datatype>& dxdt,
                        const double /* t */) const {
    auto DA1 = x[0];
    auto DM1  = x[1];
    auto E    = x[2];

    // DA1
    dxdt[0] = -g_ * DA1 + g_ * DM1;

    // DM1
    dxdt[1] = -(la_ + mu_NE_ + lc_) * DM1 +
      (mu_NE_ + la_ * E + lc_ * E * E) * DA1;

    // E
    dxdt[2] = mu_E_ - (mu_E_ + lc_) * E +
      lc_ * E * E;
  }
};

// Truncated power series for sampling = "n" (missing species): every state
// variable is a power series in z, with E(0) = z, truncated after degree
// P - 1, and stored variable by variable: x[v * P + k] is the coefficient of
// z^k of variable v. The variables are in the same order as in the
// corresponding interval structs above. The equations for the coefficients
// follow from the Cauchy product; the coefficient of degree k depends only on
// coefficients of degree <= k, so truncation is exact.
struct series_interval : public interval {
  enum kind_t { k2_NE, k2_ES, k2_EC, k3_ES, k3_NE, k4 };

  series_interval(double lc, double la, double mu_E, double mu_NE, double g,
                  int P, kind_t kind) :
    interval(lc, la, mu_E, mu_NE, g), P_(P), kind_(kind),
    EE_(P), DEE_(P), f_(P), h_(P), fa_(P), fb_(P), ha_(P) {}

  size_t size() const noexcept {
    switch (kind_) {
      case k2_NE: return 2 * P_;
      case k2_ES: return 5 * P_;
      case k2_EC: return 4 * P_;
      case k3_ES: return 7 * P_;
      case k3_NE: return 4 * P_;
      case k4:    return 3 * P_;
    }
    return 0;
  }

  // out = a * b (truncated Cauchy product)
  void conv(const double* a, const double* b, double* out) const {
    for (int k = 0; k < P_; ++k) {
      double s = 0.0;
      for (int j = 0; j <= k; ++j) s += a[j] * b[k - j];
      out[k] = s;
    }
  }

  // f = mu_NE + la E + lc E^2 (from E and EE_), stored in f_
  void make_f(const double* E) const {
    for (int k = 0; k < P_; ++k) {
      f_[k] = (k == 0 ? mu_NE_ : 0.0) + la_ * E[k] + lc_ * EE_[k];
    }
  }

  // dE = mu_E - (mu_E + lc) E + lc E^2 (from E and EE_)
  void dE(const double* E, double* d) const {
    for (int k = 0; k < P_; ++k) {
      d[k] = (k == 0 ? mu_E_ : 0.0) - (mu_E_ + lc_) * E[k] + lc_ * EE_[k];
    }
  }

  void operator()(const std::vector<double>& x,
                  std::vector<double>& dxdt,
                  const double /* t */) const {
    const int P = P_;
    auto v = [&](int i) { return x.data() + static_cast<size_t>(i) * P; };
    auto d = [&](int i) { return dxdt.data() + static_cast<size_t>(i) * P; };
    switch (kind_) {
      case k2_NE: {   // DM2, E
        const double *DM2 = v(0), *E = v(1);
        conv(E, E, EE_.data());
        for (int k = 0; k < P; ++k) d(0)[k] = -(lc_ + mu_NE_ + g_ + la_) * DM2[k];
        dE(E, d(1));
        break;
      }
      case k2_ES: {   // DE, DM2, DM3, E, DA3
        const double *DE = v(0), *DM2 = v(1), *DM3 = v(2), *E = v(3), *DA3 = v(4);
        conv(E, E, EE_.data());
        conv(DE, E, DEE_.data());
        make_f(E);
        for (int k = 0; k < P; ++k) h_[k] = la_ * DE[k] + 2.0 * lc_ * DEE_[k];
        conv(h_.data(), DA3, ha_.data());
        conv(f_.data(), DA3, fa_.data());
        for (int k = 0; k < P; ++k) {
          d(0)[k] = -(lc_ + mu_E_) * DE[k] + 2.0 * lc_ * DEE_[k];
          d(1)[k] = -(lc_ + mu_NE_ + g_ + la_) * DM2[k] + ha_[k];
          d(2)[k] = -(lc_ + mu_NE_ + la_) * DM3[k] + fa_[k];
          d(4)[k] = -g_ * DA3[k] + g_ * DM3[k];
        }
        dE(E, d(3));
        break;
      }
      case k2_EC: {   // DE, DM3, E, DA3
        const double *DE = v(0), *DM3 = v(1), *E = v(2), *DA3 = v(3);
        conv(E, E, EE_.data());
        conv(DE, E, DEE_.data());
        make_f(E);
        conv(f_.data(), DA3, fa_.data());
        for (int k = 0; k < P; ++k) {
          d(0)[k] = -(lc_ + mu_E_) * DE[k] + 2.0 * lc_ * DEE_[k];
          d(1)[k] = -(lc_ + mu_NE_ + la_) * DM3[k] + fa_[k];
          d(3)[k] = -g_ * DA3[k] + g_ * DM3[k];
        }
        dE(E, d(2));
        break;
      }
      case k3_ES: {   // DE, DM1, DM2, DM3, E, DA2, DA3
        const double *DE = v(0), *DM1 = v(1), *DM2 = v(2), *DM3 = v(3),
                     *E = v(4), *DA2 = v(5), *DA3 = v(6);
        conv(E, E, EE_.data());
        conv(DE, E, DEE_.data());
        make_f(E);
        for (int k = 0; k < P; ++k) h_[k] = la_ * DE[k] + 2.0 * lc_ * DEE_[k];
        conv(f_.data(), DA2, fa_.data());   // f DA2
        conv(f_.data(), DA3, fb_.data());   // f DA3
        conv(h_.data(), DA3, ha_.data());   // h DA3
        for (int k = 0; k < P; ++k) {
          d(0)[k] = -(lc_ + mu_E_) * DE[k] + 2.0 * lc_ * DEE_[k];
          d(1)[k] = -(lc_ + mu_NE_ + la_ + g_) * DM1[k] + g_ * DM2[k] + fa_[k];
          d(2)[k] = -(lc_ + mu_NE_ + la_) * DM2[k] + fa_[k] + ha_[k];
          d(3)[k] = -(lc_ + mu_NE_ + la_) * DM3[k] + fb_[k];
          d(5)[k] = -g_ * DA2[k] + g_ * DM2[k];
          d(6)[k] = -g_ * DA3[k] + g_ * DM3[k];
        }
        dE(E, d(4));
        break;
      }
      case k3_NE: {   // DM1, DM2, E, DA2
        const double *DM1 = v(0), *DM2 = v(1), *E = v(2), *DA2 = v(3);
        conv(E, E, EE_.data());
        make_f(E);
        conv(f_.data(), DA2, fa_.data());
        for (int k = 0; k < P; ++k) {
          d(0)[k] = -(lc_ + mu_NE_ + la_ + g_) * DM1[k] + fa_[k] + g_ * DM2[k];
          d(1)[k] = -(lc_ + mu_NE_ + la_) * DM2[k] + fa_[k];
          d(3)[k] = -g_ * DA2[k] + g_ * DM2[k];
        }
        dE(E, d(2));
        break;
      }
      case k4: {      // DA1, DM1, E
        const double *DA1 = v(0), *DM1 = v(1), *E = v(2);
        conv(E, E, EE_.data());
        make_f(E);
        conv(f_.data(), DA1, fa_.data());
        for (int k = 0; k < P; ++k) {
          d(0)[k] = -g_ * DA1[k] + g_ * DM1[k];
          d(1)[k] = -(la_ + mu_NE_ + lc_) * DM1[k] + fa_[k];
        }
        dE(E, d(2));
        break;
      }
    }
  }

  const int P_;
  const kind_t kind_;
  mutable std::vector<double> EE_, DEE_, f_, h_, fa_, fb_, ha_;
};

}   // namespace loglik
