//
//  Copyright (c) 2026 Thijs Janzen
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#include <cstdlib>    // std::getenv, std::atoi
#include <vector>
#include <array>
#include <chrono>
#include <string>
#include <utility>
#include <algorithm>
#include <memory>
#include "config.h"    // NOLINT [build/include_subdir]
#include <Rcpp.h>
#include "odeint.h"    // NOLINT [build/include_subdir]
#include "secsse_loglik.h"    // NOLINT [build/include_subdir]

template <typename ODE>
std::vector<double> solve_branch(std::unique_ptr<ODE> od,
                                 const std::vector<double>& states,
                                 const std::array<double, 2>& forTime,
                                 const std::string& method,
                                 double atol,
                                 double rtol) {
  auto t0 = std::min(forTime[0], forTime[1]);
  auto t1 = std::max(forTime[0], forTime[1]);

  auto states_out = std::vector<double>(states.begin(), states.end());

  auto workhorse = Integrator<ODE, odeintcpp::no_normalization>(
    std::move(od), method, atol, rtol);

  workhorse(states_out, t0, t1);

  return states_out;
}

template <typename ODE>
std::vector<std::vector<double>> solve_branch_times(std::unique_ptr<ODE> od,
                                                     const std::vector<double>& states,
                                                     const std::vector<double>& forTime,
                                                     const std::string& method,
                                                     double atol,
                                                     double rtol) {
  std::vector< std::vector<double > > states_out;
  std::vector<double> times(forTime.begin(), forTime.end());

  auto workhorse = Integrator<ODE, odeintcpp::no_normalization>(
    std::move(od), method, atol, rtol);

  std::vector<double> states_in(states.begin(), states.end());

  workhorse(states_in, times, &states_out);

  return states_out;
}

class pEC {
public:
  pEC(const Rcpp::NumericVector& brts,
      size_t missnumspec,
      double lambda_c,
      double lambda_a,
      double mu,
      double gamma,
      std::string method,
      double atol,
      double rtol) :
  lambda_c_(lambda_c),
  lambda_a_(lambda_a),
  mu_(mu),
  gamma_(gamma),
  method_(method),
  atol_(atol),
  rtol_(rtol)
  {
    t0 = brts[0];
    t1 = brts[1];
    t2 = brts[2];
    tp = 0;
    ti = std::vector<double>(brts.begin(), brts.end());
    std::sort(ti.begin(), ti.end());
    ti.pop_back(); ti.pop_back(); //remove last two entries
    double num_species = brts.size() - 1;
    rho = num_species / (num_species + missnumspec);
  }

  double calculate_likelihood(int stac) {

    std::vector<double> initial_conditions1 = {rho, 0.0, 1 - rho, 1.0};
    if (stac == 3) initial_conditions1      = {rho, 1.0, 1 - rho, 0.0};
    auto ti2 = ti;
    ti2.insert(ti2.begin(), 0);
    auto solution0 = solve_branch_times(std::make_unique<loglik::interval2_EC>(lambda_c_, lambda_a_, mu_, gamma_), initial_conditions1, ti2, method_, atol_, rtol_);

    for (size_t i = 1; i < ti2.size(); ++i) {
      std::array<double, 2> time = {ti2[i - 1], ti2[i]};
      auto solution1 = solve_branch(std::make_unique<loglik::interval2_EC>(lambda_c_, lambda_a_, mu_, gamma_), initial_conditions1, time, method_, atol_, rtol_);
      enum state {DE, DM3, E, DA3};
      initial_conditions1 = {lambda_c_ * solution0[i][DE] * solution1[DE], // DE = pars1[1] * solution0[, "DE"][idx + 1] * solution1[, "DE"][2],
                             0,                                            // DM3
                             solution0[i][E],                              // E = solution0[, "E"][idx + 1],
                             1                                             // DA3
                            };
    }



      std::vector<double> initial_conditions2;
      if (stac == 6) {
        enum state {DE, DM3, E, DA3};

        initial_conditions2 = { initial_conditions1[DE],                         // DE = initial_conditions1["DE"][[1]],
                                0,                                                // DM1 = 0,
                                initial_conditions1[DE] * solution0.back()[DA3], // DM2 = initial_conditions1["DE"][[1]] * solution0[, "DA3"][length(ti) + 1],
                                solution0.back()[DA3],                            // DM3 = solution0[, "DM3"][length(ti) + 1],
                                initial_conditions1[E],                           // E = initial_conditions1["E"][[1]],
                                0,                                                // DA2 = 0,
                                solution0.back()[DA3]                           // DA3 = solution0[, "DA3"][length(ti) + 1])
                              };
      } else {
        enum state {DE, DM3, E, DA3};

        initial_conditions2 = { initial_conditions1[DE],                         // DE = initial_conditions1["DE"][[1]],
                                initial_conditions1[DE] * solution0.back()[DA3], // DM2 = initial_conditions1["DE"][[1]] * solution0[, "DA3"][length(ti) + 1],
                                solution0.back()[DM3],                            // DM3 = solution0[, "DM3"][length(ti) + 1],
                                initial_conditions1[E],                           // E = initial_conditions1["E"][[1]],
                                solution0.back()[DA3]                           // DA3 = solution0[, "DA3"][length(ti) + 1])
                              };
      }



      std::array<double, 2> time2 = {t2, t1};

      std::vector<double> solution2 = stac == 6 ?
         solve_branch(std::make_unique<loglik::interval3_ES>(lambda_c_, lambda_a_, mu_, gamma_), initial_conditions2, time2, method_, atol_, rtol_)
         :
         solve_branch(std::make_unique<loglik::interval2_ES>(lambda_c_, lambda_a_, mu_, gamma_), initial_conditions2, time2, method_, atol_, rtol_);



      std::vector<double> initial_conditions3;
      if (stac == 6) {
        // solution2 returns (interval3_ES): DE, DM1, DM2, DM3, E, DA2, DA3
        enum ls {DE, DM1, DM2, DM3, E, DA2, DA3}; // local state
        initial_conditions3 = {solution2[ls::DA2],  // DA1 = solution2[, "DA2"][[2]],
                               solution2[ls::DM1],  // DM1 = solution2[, "DM1"][[2]],
                               solution2[ls::E]     // E   = solution2[, "E"][[2]]
                              };
      } else {
        // solution2 returns (interval2_ES): DE, DM2, DM3, E, DA3
        enum ls {DE, DM2, DM3, E, DA3}; // local state
        initial_conditions3 = {gamma_ * solution2[ls::DM2],  // DA1 = pars1[4] * solution2[, "DM2"][[2]],
                               gamma_ * solution2[ls::DM2], // DM1 = pars1[4] * solution2[, "DM2"][[2]],
                               solution2[ls::E]  // E = solution2[, "E"][[2]])
                              };
      }

      std::array<double, 2> time3 = {t1, t0};

      auto solution3 = solve_branch(std::make_unique<loglik::interval4>(lambda_c_, lambda_a_, mu_, gamma_), initial_conditions3, time3, method_, atol_, rtol_);

      // interval4 returns: DA1, DM1, E
      enum class ls {DA1, DM1, E}; // local state
      auto prob = solution3[static_cast<size_t>(ls::DA1)];
      return std::log(prob);
  }

  private:
    double t0, t1, t2, tp;
    std::vector<double> ti;
    double rho;

    const double lambda_c_;
    const double lambda_a_;
    const double mu_;
    const double gamma_;
    const std::string method_;
    const double atol_;
    const double rtol_;
};


double DAISIE_DE_logpEC_general(const Rcpp::NumericVector& brts,
                                size_t missnumspec,
                                double lambda_c,
                                double lambda_a,
                                double mu,
                                double gamma,
                                size_t stac,
                                std::string method,
                                double rtol,
                                double atol) {
    auto solver = pEC(brts,
                      missnumspec,
                      lambda_c,
                      lambda_a,
                      mu,
                      gamma,
                      method,
                      atol,
                      rtol);

    auto loglik = solver.calculate_likelihood(stac);
    return loglik;
}


//' Wrapper for the DAISIE_DE general integrator
 //'
 //' @description This is the rcpp function to do single branch DAISIE_DE calculations
 //' @name DAISIE_DE_logpEC_general_rcpp
 //' @export DAISIE_DE_logpEC_general_rcpp
 //' @return list
 RcppExport SEXP DAISIE_DE_logpEC_general_rcpp(SEXP brtsSEXP, SEXP missnumspecSEXP,
                                               SEXP lambda_cSEXP, SEXP lambda_aSEXP, SEXP muSEXP, SEXP gammaSEXP,
                                               SEXP stacSEXP,
                                               SEXP inte_methodSEXP,
                                               SEXP atolSEXP, SEXP rtolSEXP) {
   BEGIN_RCPP
   Rcpp::RObject rcpp_result_gen;

   Rcpp::traits::input_parameter< Rcpp::NumericVector >::type brts(brtsSEXP);
   Rcpp::traits::input_parameter< size_t >::type missnumspec(missnumspecSEXP);

   Rcpp::traits::input_parameter< double >::type lambda_c(lambda_cSEXP);
   Rcpp::traits::input_parameter< double >::type lambda_a(lambda_aSEXP);
   Rcpp::traits::input_parameter< double >::type mu(muSEXP);
   Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);

   Rcpp::traits::input_parameter< size_t >::type stac(stacSEXP);

   Rcpp::traits::input_parameter< std::string >::type inte_method(inte_methodSEXP);

   Rcpp::traits::input_parameter< double >::type atol(atolSEXP);
   Rcpp::traits::input_parameter< double >::type rtol(rtolSEXP);

   rcpp_result_gen = Rcpp::wrap(DAISIE_DE_logpEC_general(brts, missnumspec, lambda_c, lambda_a, mu, gamma, stac, inte_method, atol, rtol));
   return rcpp_result_gen;
   END_RCPP
 }
