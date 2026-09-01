//  Copyright (c) 2021 - 2023, Thijs Janzen
//  Copyright (c) 2023, Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#pragma once
#include <Rcpp.h>
#include <RcppParallel.h>
#include <type_traits>
#include <vector>
#include <utility>
#include <string>
#include "TRAISIE_sq_matrix.h"                  // NOLINT [build/include_subdir]



namespace TRAISIE {

template <typename T>
using rvector = RcppParallel::RVector<T>;

template <typename T>
using rmatrix = RcppParallel::RMatrix<T>;

inline std::vector<double> make_dist_g(const Rcpp::NumericVector& tma,
                                       const Rcpp::NumericVector& gamma,
                                       const size_t num_unique_states) {
  std::vector<double> dist_gamma(gamma.begin(), gamma.end());

  if (num_unique_states == 1) return(dist_gamma);

  if (Rcpp::NumericVector::is_na(tma[0])) {
    // tma not known, please note that entire vector is NA in this case
    return(dist_gamma);
  } else {
    // no NAs, we are good:
    for (size_t i = 0; i < dist_gamma.size(); ++i) {
      dist_gamma[i] *= tma[i];
    }
  }

  return dist_gamma;
}

inline double calc_sum_dist_g_(const std::vector<double>& dist_gamma) {
  return std::accumulate(dist_gamma.begin(), dist_gamma.end(), 0.0);
}

inline double calc_sum(const std::vector<double>& dist_g_,
                       const vector_view_t<const double>& DM3) {
  double s = 0.0;
  for (size_t i = 0; i < dist_g_.size(); ++i) {
    s += dist_g_[i] * DM3[i];
  }
  return s;
}

struct interval {
  const rvector<const double> lc_;   // cladogenesis rates
  const rvector<const double> m_;    //  extinction rates

  const rvector<const double> la_;   // anagenesis rates
  const sq_matrix q_;                // transition rates
  const sq_matrix p_;

  const size_t n_;                  // number of unique states

  const std::vector<double> t_vec;

  const std::vector<double> dist_g_;

  const double sum_dist_g_;

  // constructor
  interval(const Rcpp::NumericVector& lc,
            const Rcpp::NumericVector& la,
            const Rcpp::NumericVector& m,
            const Rcpp::NumericVector& g,
            const Rcpp::NumericMatrix& q,
            const Rcpp::NumericMatrix& p,
            const Rcpp::NumericVector& tma,
            const size_t n)
    : lc_(lc),
      m_(m),
      la_(la),
      q_(q),
      p_(p),
      n_(n),
      t_vec(q_.row_sums()),
      dist_g_(make_dist_g(tma, g, n)),
      sum_dist_g_(calc_sum_dist_g_(dist_g_)) {

     if (p_.size() != q_.size())
       throw "p and q are not same dimensions";
  }
};


struct interval1 : public interval {
  using interval::interval;

  size_t size() const noexcept {
    // (DE + DM3 + E) * n + DA3
    return 3 * n_ + 1;
  }

  void mergebranch(const std::vector<double>& N,
                   const std::vector<double>& M,
                   std::vector<double>& out) const {
    out = N;

    for (size_t i = 0; i < n_; ++i) {
      // out[DE_0]  = lc_[0] * N[DE_0] * M[DE_0]
      out[i] = lc_[i] * N[i] * M[i];
    }
  }

  // this is the dx/dt calculation // true rhs that gets integrated
  // along the branches
  void operator()(const std::vector<double>& x,
                std::vector<double>& dxdt,
                const double /* t */) const {
    auto DA3 = x.back();

    auto DE  = vector_view_t<const double>(x.data() , n_);
    auto DM3 = vector_view_t<const double>(x.data() + 1 * n_, n_);
    auto E   = vector_view_t<const double>(x.data() + 2 * n_, n_);

    auto q_mult_E   = q_ * E;
    auto q_mult_DE  = q_ * DE;
    auto q_mult_DM3 = q_ * DM3;

    double s_g_DM3 = calc_sum(dist_g_, DM3);

    static sq_matrix pq(n_);
    static sq_matrix opq(n_);

    element_mult(p_, q_, &pq);
    element_mult_one_minus(p_, q_, &opq);

    static std::vector<double> pq_mult_E(n_);     vector_view_mult(pq, E, &pq_mult_E);
    static std::vector<double> opq_mult_DM3(n_);  vector_view_mult(opq, DM3, &opq_mult_DM3);

    for (size_t i = 0; i < n_; ++i) {
      auto lambda_c_mu_t_vec_sum = lc_[i] + m_[i] + t_vec[i];

      // DE
      dxdt[i] = -(lambda_c_mu_t_vec_sum) * DE[i] +
        2 * lc_[i] * DE[i] * E[i] +
        q_mult_DE[i];
      // DM3
      dxdt[i + n_] =
       -(lambda_c_mu_t_vec_sum + sum_dist_g_ + la_[i]) * DM3[i] +
       (m_[i] + la_[i] * E[i] + lc_[i] * E[i] * E[i] + pq_mult_E[i]) * DA3 +
       opq_mult_DM3[i] +
       s_g_DM3;
      // E
      dxdt[i + n_ + n_] = m_[i] - (lambda_c_mu_t_vec_sum) * E[i] +
        lc_[i] * E[i] * E[i] +
        q_mult_E[i];
    }

    // DA3
    dxdt.back() = -sum_dist_g_ * DA3 + s_g_DM3;
  }
};

struct interval2 : public interval {
  using interval::interval;

  size_t size() const noexcept {
    // (DE + DM3 + E) * n + DA3
    return 4 * n_ + 1;
  }

  // this is the dx/dt calculation // true rhs that gets integrated
  // along the branches
  void operator()(const std::vector<double>& x,
                std::vector<double>& dxdt,
                const double /* t */) const {
    auto DA3 = x.back();

    auto DE  = vector_view_t<const double>(x.data() + 0 * n_, n_);
    auto DM2 = vector_view_t<const double>(x.data() + 1 * n_, n_);
    auto DM3 = vector_view_t<const double>(x.data() + 2 * n_, n_);
    auto E   = vector_view_t<const double>(x.data() + 3 * n_, n_);

    auto q_mult_E   = q_ * E;
    auto q_mult_DE  = q_ * DE;
    auto q_mult_DM2 = q_ * DM2;
    auto q_mult_DM3 = q_ * DM3;

    double s_g_DM3 = calc_sum(dist_g_, DM3);

    static sq_matrix pq(n_);
    static sq_matrix opq(n_);

    element_mult(p_, q_, &pq);
    element_mult_one_minus(p_, q_, &opq);

    static std::vector<double> pq_mult_E(n_);     vector_view_mult(pq, E, &pq_mult_E);
    static std::vector<double> pq_mult_DE(n_);    vector_view_mult(pq, DE, &pq_mult_DE);
    static std::vector<double> opq_mult_DM2(n_);  vector_view_mult(opq, DM2, &opq_mult_DM2);
    static std::vector<double> opq_mult_DM3(n_);  vector_view_mult(opq, DM3, &opq_mult_DM3);

    for (size_t i = 0; i < n_; ++i) {
      auto lambda_c_mu_t_vec_sum = lc_[i] + m_[i] + t_vec[i];

      // DE
      dxdt[i] = -(lambda_c_mu_t_vec_sum) * DE[i] +
        2 * lc_[i] * DE[i] * E[i] +
        q_mult_DE[i];

      // DM2
      dxdt[i + n_] =
       -(lambda_c_mu_t_vec_sum + sum_dist_g_ + la_[i]) * DM2[i] +
       (la_[i] * DE[i] + 2 * lc_[i] * DE[i] * E[i] + pq_mult_DE[i]) * DA3 +
       opq_mult_DM2[i];

      // DM3
      dxdt[i + n_ + n_] =
      -(lambda_c_mu_t_vec_sum + sum_dist_g_ + la_[i]) * DM3[i] +
       (m_[i] + la_[i] * E[i] + lc_[i] * E[i] * E[i] + pq_mult_E[i]) * DA3 +
       opq_mult_DM3[i] +
       s_g_DM3;

      // E
      dxdt[i + n_ + n_ + n_] = m_[i] - (lambda_c_mu_t_vec_sum) * E[i] +
        lc_[i] * E[i] * E[i] +
        q_mult_E[i];
    }

    // DA3
    dxdt.back() = -sum_dist_g_ * DA3 + s_g_DM3;
  }
};





struct interval3 : public interval {
  using interval::interval;



  size_t size() const noexcept {
    // (DE + DM3 + E) * n + DA3
    return 5 * n_ + 2;
  }
  void operator()(const std::vector<double>& x,
                  std::vector<double>& dxdt,
                  const double /* t */) {
    auto DA3 = x.back();
    auto DA2 = x[x.size() - 2];

    auto DE  = vector_view_t<const double>(x.data() , n_);
    auto DM1 = vector_view_t<const double>(x.data() + 1 * n_, n_);
    auto DM2 = vector_view_t<const double>(x.data() + 2 * n_, n_);
    auto DM3 = vector_view_t<const double>(x.data() + 3 * n_, n_);
    auto E   = vector_view_t<const double>(x.data() + 4 * n_, n_);

    auto q_mult_E   = q_ * E;
    auto q_mult_DE  = q_ * DE;
    auto q_mult_DM1 = q_ * DM1;
    auto q_mult_DM2 = q_ * DM2;
    auto q_mult_DM3 = q_ * DM3;

    double s_g_DM3 = calc_sum(dist_g_, DM3);
    double s_g_DM2 = calc_sum(dist_g_, DM2);
    
    static sq_matrix pq(n_);
    static sq_matrix opq(n_);

    element_mult(p_, q_, &pq);
    element_mult_one_minus(p_, q_, &opq);

    static std::vector<double> pq_mult_E(n_);     vector_view_mult(pq, E, &pq_mult_E);
    static std::vector<double> pq_mult_DE(n_);    vector_view_mult(pq, DE, &pq_mult_DE);
    static std::vector<double> opq_mult_DM1(n_);  vector_view_mult(opq, DM1, &opq_mult_DM1);
    static std::vector<double> opq_mult_DM2(n_);  vector_view_mult(opq, DM2, &opq_mult_DM2);
    static std::vector<double> opq_mult_DM3(n_);  vector_view_mult(opq, DM3, &opq_mult_DM3);

    for (size_t i = 0; i < n_; ++i) {
      auto lambda_c_mu_t_vec_sum = lc_[i] + m_[i] + t_vec[i];

      // DE
      dxdt[i] = -(lambda_c_mu_t_vec_sum) * DE[i] +
        2 * lc_[i] * DE[i] * E[i] +
        q_mult_DE[i];

      // DM1
      dxdt[i + 1 * n_] =
       -(lambda_c_mu_t_vec_sum + sum_dist_g_ + la_[i]) * DM1[i] +
       (m_[i] + la_[i] * E[i] + lc_[i] * E[i] * E[i] + pq_mult_E[i]) * DA2 +
       opq_mult_DM1[i] +
       s_g_DM2;

      // DM2
      dxdt[i + 2 * n_] =
       -(lambda_c_mu_t_vec_sum + sum_dist_g_ + la_[i]) * DM2[i] +
       (m_[i] + la_[i] * E[i] + lc_[i] * E[i] * E[i] + pq_mult_E[i]) * DA2 +
       (la_[i] * DE[i] + 2 * lc_[i] * DE[i]  + pq_mult_DE[i]) * DA3 +
       opq_mult_DM2[i] +
       s_g_DM2;

      // DM3
      dxdt[i + 3 * n_] =
      -(lambda_c_mu_t_vec_sum + sum_dist_g_ + la_[i]) * DM3[i] +
       (m_[i] + la_[i] * E[i] + lc_[i] * E[i] * E[i] + pq_mult_E[i]) * DA3 +
       opq_mult_DM3[i] +
       s_g_DM3;

      // E
      dxdt[i + 4 * n_] = m_[i] - (lambda_c_mu_t_vec_sum) * E[i] +
        lc_[i] * E[i] * E[i] +
        q_mult_E[i];
    }

    // DA3
    dxdt.back() = -sum_dist_g_ * DA3 + s_g_DM3;

    // DA2
    dxdt[dxdt.size() - 2] = -sum_dist_g_ * DA2 + s_g_DM2;
  }

};

struct interval4 : public interval {
  using interval::interval;

  size_t size() const noexcept {
    return 2 * n_ + 1;
  }

  void operator()(const std::vector<double>& x,
                std::vector<double>& dxdt,
                const double /* t */) const {
    auto DA1 = x.back();

    auto DM1  = vector_view_t<const double>(x.data() , n_);
    auto E   = vector_view_t<const double>(x.data() + n_, n_);

    auto q_mult_E   = q_ * E;
    auto q_mult_DM1  = q_ * DM1;

    double s_g_DM1 = calc_sum(dist_g_, DM1);

    sq_matrix pq = element_mult(p_, q_);
    sq_matrix opq = element_mult_one_minus(p_, q_);

    auto pq_mult_E    = pq * E;
    auto opq_mult_DM1 = opq * DM1;

    for (size_t i = 0; i < n_; ++i) {
      auto lambda_c_mu_t_vec_sum = lc_[i] + m_[i] + t_vec[i];

      // DM1
    //  dxdt[i] =
    //  -(lambda_c_mu_t_vec_sum + sum_dist_g_ + la_[i]) * DM1[i] +
     //  (m_[i] + la_[i] * E[i] + lc_[i] * E[i] * E[i] + p_ * q_mult_E[i]) * DA1 +
     //  (1 - p_) * q_mult_DM1[i] +
      // s_g_DM1;

      // DM1
      dxdt[i] =
          -(lambda_c_mu_t_vec_sum  + la_[i]) * DM1[i] +
          (m_[i] + la_[i] * E[i] + lc_[i] * E[i] * E[i] + pq_mult_E[i]) * DA1 +
          opq_mult_DM1[i];

      // E
      dxdt[i + n_] = m_[i] - (lambda_c_mu_t_vec_sum) * E[i] +
        lc_[i] * E[i] * E[i] +
        q_mult_E[i];
    }

    // DA3
    dxdt.back() = -sum_dist_g_ * DA1 + s_g_DM1;
  }
};

}   // namespace loglik
