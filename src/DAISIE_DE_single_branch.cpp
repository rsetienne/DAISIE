//
//  Copyright (c) 2025 Thijs Janzen
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#include <cstdlib>    // std::getenv, std::atoi
#include <vector>
#include <chrono>
#include <string>
#include <utility>
#include <algorithm>
#include <memory>
#include "config.h"    // NOLINT [build/include_subdir]
#include <Rcpp.h>

#include <complex>

#include "DAISIE_DE_rhs.h"   // NOLINT [build/include_subdir]
#include "DAISIE_DE_odeint.h"    // NOLINT [build/include_subdir]

template <typename ODE, typename datatype>
std::vector<datatype> solve_branch(std::unique_ptr<ODE> od,
                                 const std::vector<datatype>& states,
                                 const std::array<double, 2>& forTime,
                                 const std::string& method,
                                 double atol,
                                 double rtol) {
  auto t0 = std::min(forTime[0], forTime[1]);
  auto t1 = std::max(forTime[0], forTime[1]);

  auto states_out = std::vector<datatype>(states.begin(), states.end());

  auto workhorse = Integrator<ODE, odeintcpp::no_normalization, datatype>(std::move(od), method, atol, rtol);

  workhorse(states_out, t0, t1);

  return states_out;
}

template <typename ODE, typename datatype>
std::vector<std::vector<datatype>> solve_branch_times(std::unique_ptr<ODE> od,
                                                    const std::vector<datatype>& states,
                                                    const std::vector<double>& forTime,
                                                    const std::string& method,
                                                    double atol,
                                                    double rtol) {
  std::vector< std::vector< datatype > > states_out;
  std::vector<double> times(forTime.begin(), forTime.end());

  auto workhorse = Integrator<ODE, odeintcpp::no_normalization, datatype>(std::move(od), method, atol, rtol);

  std::vector<datatype> states_in(states.begin(), states.end());

  workhorse(states_in, times, &states_out);

  return states_out;
}

template <typename ODE, typename datatype>
Rcpp::List calc_ll_single_branch(std::unique_ptr<ODE> od,
                                 const std::vector<datatype>& states,
                                 const std::vector<double>& forTime,
                                 const std::string& method,
                                 double atol,
                                 double rtol) {
    if (forTime.size() == 2) {
      auto states_out = solve_branch(std::move(od), states, {forTime[0], forTime[1]}, method, atol, rtol);
      return Rcpp::List::create(Rcpp::Named("states") = states_out);
    } else if (forTime.size() < 2) {
      throw std::invalid_argument("forTime should have at least 2 entries");
    }

    // and if forTime is a vector:
    auto states_out = solve_branch_times(std::move(od), states, std::vector<double>(forTime.begin(), forTime.end()), method, atol, rtol);

    return Rcpp::List::create(Rcpp::Named("states") = states_out);

}

enum class string_code {
  interval2_NE,
  interval2_ES,
  interval2_EC,
  interval3_ES,
  interval3_NE,
  interval4
};

string_code hash_string(const std::string& s) {
  if (s == "interval2_NE") return string_code::interval2_NE;
  if (s == "interval2_ES") return string_code::interval2_ES;
  if (s == "interval2_EC") return string_code::interval2_EC;
  if (s == "interval3_ES") return string_code::interval3_ES;
  if (s == "interval3_NE") return string_code::interval3_NE;
  if (s == "interval4")    return string_code::interval4;

  return string_code::interval4;
}


template <typename datatype>
Rcpp::List DAISIE_DE_cpp_solve_local(const double& lambda_c,
                                     const double& lambda_a,
                                     const double& mu_E,
                                     const double& mu_NE,
                                     const double& gamma,
                                     const std::string& chosen_interval,
                                     const std::string& inte_method,
                                     const std::vector<datatype>& init_states,
                                     const std::vector<double>& time,
                                     double atol,
                                     double rtol) {

  switch( hash_string(chosen_interval)) {
    case string_code::interval2_NE:
      return calc_ll_single_branch(std::make_unique<loglik::interval2_NE<datatype>>(lambda_c, lambda_a, mu_E, mu_NE, gamma), init_states, time, inte_method, atol, rtol);
    case string_code::interval2_ES:
      return calc_ll_single_branch(std::make_unique<loglik::interval2_ES<datatype>>(lambda_c, lambda_a, mu_E, mu_NE, gamma), init_states, time, inte_method, atol, rtol);
    case string_code::interval2_EC:
      return calc_ll_single_branch(std::make_unique<loglik::interval2_EC<datatype>>(lambda_c, lambda_a, mu_E, mu_NE, gamma), init_states, time, inte_method, atol, rtol);
    case string_code::interval3_ES:
      return calc_ll_single_branch(std::make_unique<loglik::interval3_ES<datatype>>(lambda_c, lambda_a, mu_E, mu_NE, gamma), init_states, time, inte_method, atol, rtol);
    case string_code::interval3_NE:
      return calc_ll_single_branch(std::make_unique<loglik::interval3_NE<datatype>>(lambda_c, lambda_a, mu_E, mu_NE, gamma), init_states, time, inte_method, atol, rtol);
    case string_code::interval4:
      return calc_ll_single_branch(std::make_unique<loglik::interval4<datatype>   >(lambda_c, lambda_a, mu_E, mu_NE, gamma), init_states, time, inte_method, atol, rtol);
  }
  return NA_REAL;
}

std::vector<std::complex<double>> as_std_vector(const Rcpp::ComplexVector& x) {
    std::vector<std::complex<double>> result;
    result.reserve(x.size());

    for (size_t i = 0; i < x.size(); ++i) {
        result.emplace_back(x[i].r, x[i].i);
    }

    return result;
}

Rcpp::List DAISIE_DE_cpp_solve_local_double(const double& lambda_c,
                                     const double& lambda_a,
                                     const double& mu_E,
                                     const double& mu_NE,
                                     const double& gamma,
                                     const std::string& chosen_interval,
                                     const std::string& inte_method,
                                     const Rcpp::NumericVector& init_states,
                                     const Rcpp::NumericVector& time,
                                     double atol,
                                     double rtol) {
  auto init_states_vec = Rcpp::as<std::vector<double>>(init_states);
  auto time_vec = Rcpp::as<std::vector<double>>(time);

  return DAISIE_DE_cpp_solve_local<double>(lambda_c, lambda_a, mu_E, mu_NE, gamma, chosen_interval, inte_method, init_states_vec, time_vec, atol, rtol);
}

Rcpp::List DAISIE_DE_cpp_solve_local_complex(const double& lambda_c,
                                     const double& lambda_a,
                                     const double& mu_E,
                                     const double& mu_NE,
                                     const double& gamma,
                                     const std::string& chosen_interval,
                                     const std::string& inte_method,
                                     const Rcpp::ComplexVector& init_states,
                                     const Rcpp::NumericVector& time,
                                     double atol,
                                     double rtol) {
  auto init_states_vec = as_std_vector(init_states);
  auto time_vec = Rcpp::as<std::vector<double>>(time);

  return DAISIE_DE_cpp_solve_local<loglik::complex>(lambda_c, lambda_a, mu_E, mu_NE, gamma, chosen_interval, inte_method, init_states_vec, time_vec, atol, rtol);
}

RcppExport SEXP DAISIE_DE_cpp_solve(SEXP lambda_cSEXP, SEXP lambda_aSEXP, SEXP mu_ESEXP, SEXP mu_NESEXP, SEXP gammaSEXP,
                                    SEXP chosen_intervalSEXP, SEXP inte_methodSEXP,
                                    SEXP init_statesSEXP, SEXP timeSEXP,
                                    SEXP atolSEXP, SEXP rtolSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< double >::type lambda_c(lambda_cSEXP);
  Rcpp::traits::input_parameter< double >::type lambda_a(lambda_aSEXP);
  Rcpp::traits::input_parameter< double >::type mu_E(mu_ESEXP);
  Rcpp::traits::input_parameter< double >::type mu_NE(mu_NESEXP);
  Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);

  Rcpp::traits::input_parameter< std::string >::type chosen_interval(chosen_intervalSEXP);
  Rcpp::traits::input_parameter< std::string >::type inte_method(inte_methodSEXP);

  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type init_states(init_statesSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type time(timeSEXP);

  Rcpp::traits::input_parameter< double >::type atol(atolSEXP);
  Rcpp::traits::input_parameter< double >::type rtol(rtolSEXP);

  rcpp_result_gen = Rcpp::wrap(DAISIE_DE_cpp_solve_local_double(lambda_c, lambda_a, mu_E, mu_NE, gamma, chosen_interval, inte_method, init_states, time, atol, rtol));
  return rcpp_result_gen;
  END_RCPP
}


RcppExport SEXP DAISIE_DE_cpp_solve_complex(SEXP lambda_cSEXP, SEXP lambda_aSEXP, SEXP mu_ESEXP, SEXP mu_NESEXP, SEXP gammaSEXP,
                                    SEXP chosen_intervalSEXP, SEXP inte_methodSEXP,
                                    SEXP init_statesSEXP, SEXP timeSEXP,
                                    SEXP atolSEXP, SEXP rtolSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< double >::type lambda_c(lambda_cSEXP);
  Rcpp::traits::input_parameter< double >::type lambda_a(lambda_aSEXP);
  Rcpp::traits::input_parameter< double >::type mu_E(mu_ESEXP);
  Rcpp::traits::input_parameter< double >::type mu_NE(mu_NESEXP);
  Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);

  Rcpp::traits::input_parameter< std::string >::type chosen_interval(chosen_intervalSEXP);
  Rcpp::traits::input_parameter< std::string >::type inte_method(inte_methodSEXP);

  Rcpp::traits::input_parameter< Rcpp::ComplexVector >::type init_states(init_statesSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type time(timeSEXP);

  Rcpp::traits::input_parameter< double >::type atol(atolSEXP);
  Rcpp::traits::input_parameter< double >::type rtol(rtolSEXP);

  rcpp_result_gen = Rcpp::wrap(DAISIE_DE_cpp_solve_local_complex(lambda_c, lambda_a, mu_E, mu_NE, gamma, chosen_interval, inte_method, init_states, time, atol, rtol));
  return rcpp_result_gen;
  END_RCPP
}
