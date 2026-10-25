//
//  Copyright (c) 2026 Thijs Janzen
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
#include <complex>
#include "config.h"    // NOLINT [build/include_subdir]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "TRAISIE_loglik.h"    // NOLINT [build/include_subdir]

// datatype is double (real numbers) or std::complex<double> (complex numbers)
template <typename ODE, typename datatype>
Rcpp::List TRAISIE_calc_ll_single_branch(std::unique_ptr<ODE> od,
                                         const std::vector<datatype>& states,
                                         const Rcpp::NumericVector& forTime,
                                         const std::string& method,
                                         double atol,
                                         double rtol) {
  try {
    auto t0 = std::min(forTime[0], forTime[1]);
    auto t1 = std::max(forTime[0], forTime[1]);

    auto T0 = std::chrono::high_resolution_clock::now();

    auto states_out = std::vector<datatype>(states.begin(), states.end());

    auto workhorse = Integrator<ODE, datatype>(
      std::move(od), method, atol, rtol);

    workhorse(states_out, t0, t1);

    auto T1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> DT = (T1 - T0);


    return Rcpp::List::create(Rcpp::Named("states") = states_out,
                              Rcpp::Named("duration") = DT.count());
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (const char* msg) {
    Rcpp::Rcout << msg << std::endl;
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}

template <typename datatype>
Rcpp::List TRAISIE_branch_rcpp(const Rcpp::NumericVector& lambda_cs,
                               const Rcpp::NumericVector& lambda_as,
                               const Rcpp::NumericVector& mus,
                               const Rcpp::NumericVector& gammas,
                               const Rcpp::NumericMatrix& qs,
                               const Rcpp::NumericMatrix& p,
                               const Rcpp::NumericVector& tma,
                               const std::string& chosen_interval,
                               const std::string& inte_method,
                               const std::vector<datatype>& init_states,
                               const Rcpp::NumericVector& time,
                               double atol,
                               double rtol) {
  auto num_unique_states = lambda_cs.size();

  if (chosen_interval == "TRAISIE_interval2") {
    return TRAISIE_calc_ll_single_branch(
      std::make_unique<TRAISIE::interval2<datatype>>(lambda_cs,
                                           lambda_as,
                                           mus,
                                           gammas,
                                           qs,
                                           p,
                                           tma,
                                           num_unique_states),
                                           init_states,
                                           time,
                                           inte_method,
                                           atol,
                                           rtol);
  }
  if (chosen_interval == "TRAISIE_interval3") {
    return TRAISIE_calc_ll_single_branch(
      std::make_unique<TRAISIE::interval3<datatype>>(lambda_cs,
                                           lambda_as,
                                           mus,
                                           gammas,
                                           qs,
                                           p,
                                           tma,
                                           num_unique_states),
                                           init_states,
                                           time,
                                           inte_method,
                                           atol,
                                           rtol);
  }

  if (chosen_interval == "TRAISIE_interval4") {
    return TRAISIE_calc_ll_single_branch(
      std::make_unique<TRAISIE::interval4<datatype>>(lambda_cs,
                                           lambda_as,
                                           mus,
                                           gammas,
                                           qs,
                                           p,
                                           tma,
                                           num_unique_states),
                                           init_states,
                                           time,
                                           inte_method,
                                           atol,
                                           rtol);
  }

  return NA_REAL;
}



RcppExport SEXP TRAISIE_branch_cpp(SEXP lambda_csSEXP, SEXP lambda_asSEXP, SEXP musSEXP, SEXP gammasSEXP, SEXP qsSEXP,
                                   SEXP pSEXP, SEXP tmaSEXP,
                                   SEXP chosen_intervalSEXP, SEXP inte_methodSEXP,
                                   SEXP init_statesSEXP, SEXP timeSEXP,
                                   SEXP atolSEXP, SEXP rtolSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda_cs(lambda_csSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda_as(lambda_asSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mus(musSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type gammas(gammasSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type qs(qsSEXP);

  Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type p(pSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tma(tmaSEXP);

  Rcpp::traits::input_parameter< std::string >::type chosen_interval(chosen_intervalSEXP);
  Rcpp::traits::input_parameter< std::string >::type inte_method(inte_methodSEXP);

  Rcpp::traits::input_parameter< std::vector<double> >::type init_states(init_statesSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type time(timeSEXP);

  Rcpp::traits::input_parameter< double >::type atol(atolSEXP);
  Rcpp::traits::input_parameter< double >::type rtol(rtolSEXP);

  rcpp_result_gen = Rcpp::wrap(TRAISIE_branch_rcpp<double>(lambda_cs, lambda_as, mus, gammas, qs, p, tma, chosen_interval, inte_method, init_states, time, atol, rtol));
  return rcpp_result_gen;
  END_RCPP
}


// copy an R complex vector into a C++ vector of complex numbers
inline std::vector<std::complex<double>> TRAISIE_as_complex_vector(const Rcpp::ComplexVector& x) {
  std::vector<std::complex<double>> result;
  result.reserve(x.size());
  for (R_xlen_t i = 0; i < x.size(); ++i) {
    result.emplace_back(x[i].r, x[i].i);
  }
  return result;
}

// same as TRAISIE_branch_cpp, but the start values are complex numbers
// (used when missing species are given, see TRAISIE_n)
RcppExport SEXP TRAISIE_branch_cpp_complex(SEXP lambda_csSEXP, SEXP lambda_asSEXP, SEXP musSEXP, SEXP gammasSEXP, SEXP qsSEXP,
                                           SEXP pSEXP, SEXP tmaSEXP,
                                           SEXP chosen_intervalSEXP, SEXP inte_methodSEXP,
                                           SEXP init_statesSEXP, SEXP timeSEXP,
                                           SEXP atolSEXP, SEXP rtolSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda_cs(lambda_csSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda_as(lambda_asSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mus(musSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type gammas(gammasSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type qs(qsSEXP);

  Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type p(pSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tma(tmaSEXP);

  Rcpp::traits::input_parameter< std::string >::type chosen_interval(chosen_intervalSEXP);
  Rcpp::traits::input_parameter< std::string >::type inte_method(inte_methodSEXP);

  Rcpp::traits::input_parameter< Rcpp::ComplexVector >::type init_states(init_statesSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type time(timeSEXP);

  Rcpp::traits::input_parameter< double >::type atol(atolSEXP);
  Rcpp::traits::input_parameter< double >::type rtol(rtolSEXP);

  rcpp_result_gen = Rcpp::wrap(TRAISIE_branch_rcpp<std::complex<double>>(lambda_cs, lambda_as, mus, gammas, qs, p, tma, chosen_interval, inte_method, TRAISIE_as_complex_vector(init_states), time, atol, rtol));
  return rcpp_result_gen;
  END_RCPP
}
