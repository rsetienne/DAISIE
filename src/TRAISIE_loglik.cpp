// Copyright 2023 Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)


#include <cstdlib>    // std::getenv, std::atoi
#include <vector>
#include <chrono>
#include <utility>
#include <algorithm>
#include <string>
#include <complex>

#include "config.h"    // NOLINT [build/include_subdir]
#include <Rcpp.h>
#include <RcppParallel.h>
#include "TRAISIE_loglik.h"    // NOLINT [build/include_subdir]
#include "tbb_stub.h"

// copy each row of the R matrix into a C++ vector (real numbers)
inline std::vector<std::vector<double>> TRAISIE_matrix_rows(const Rcpp::NumericMatrix& states) {
  std::vector<std::vector<double>> tstates{};
  for (int i = 0; i < states.nrow(); ++i) {
    tstates.emplace_back(states.row(i).begin(), states.row(i).end());
  }
  return tstates;
}

// copy each row of the R matrix into a C++ vector (complex numbers)
inline std::vector<std::vector<std::complex<double>>> TRAISIE_matrix_rows(const Rcpp::ComplexMatrix& states) {
  std::vector<std::vector<std::complex<double>>> tstates(states.nrow());
  for (int i = 0; i < states.nrow(); ++i) {
    for (int j = 0; j < states.ncol(); ++j) {
      tstates[i].emplace_back(states(i, j).r, states(i, j).i);
    }
  }
  return tstates;
}

// copy the C++ rows back into the R matrix (real numbers)
inline void TRAISIE_copy_rows_back(const std::vector<std::vector<double>>& tstates,
                                   Rcpp::NumericMatrix& states_out) {
  for (int i = 0; i < states_out.nrow(); ++i) {
    std::copy(std::begin(tstates[i]), std::end(tstates[i]),
              states_out.row(i).begin());
  }
}

// copy the C++ rows back into the R matrix (complex numbers)
inline void TRAISIE_copy_rows_back(const std::vector<std::vector<std::complex<double>>>& tstates,
                                   Rcpp::ComplexMatrix& states_out) {
  for (int i = 0; i < states_out.nrow(); ++i) {
    for (size_t j = 0; j < tstates[i].size(); ++j) {
      Rcomplex value;
      value.r = tstates[i][j].real();
      value.i = tstates[i][j].imag();
      states_out(i, j) = value;
    }
  }
}

// datatype is double (real numbers) or std::complex<double> (complex numbers)
// RMATRIX is Rcpp::NumericMatrix or Rcpp::ComplexMatrix
template <typename ODE, typename datatype, typename RMATRIX>
Rcpp::List TRAISIE_calc_ll(std::unique_ptr<ODE> od,
                           const Rcpp::IntegerVector& ances,
                           const RMATRIX& states,
                           const Rcpp::NumericMatrix& forTime,
                           const std::string& method,
                           double atol,
                           double rtol,
                           bool see_states,
                           bool use_normalization) {
  auto num_threads = get_rcpp_num_threads();

  auto T0 = std::chrono::high_resolution_clock::now();
  std::vector<std::vector<datatype>> tstates = TRAISIE_matrix_rows(states);
  const auto phy_edge = make_phy_edge_vector(TRAISIE::rmatrix<const double>(forTime));

  auto inodes = find_inte_nodes(phy_edge,
                                TRAISIE::rvector<const int>(ances),
                                tstates, num_threads);

  calc_ll_res<datatype> ll_res = calc_ll(TreeIntegrator<ODE, datatype>(std::move(od),
                                                                      method,
                                                                      atol,
                                                                      rtol),
                                                                      inodes,
                                                                      tstates,
                                                                      num_threads);

  auto T1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> DT = (T1 - T0);
  RMATRIX states_out;
  if (see_states) {
    // R side expect full states back.
    states_out = RMATRIX(states.nrow(), states.ncol());
    TRAISIE_copy_rows_back(tstates, states_out);
  }
  return Rcpp::List::create(Rcpp::Named("loglik") = ll_res.loglik,
                            Rcpp::Named("node_M") = ll_res.node_M,
                            Rcpp::Named("merge_branch") = ll_res.merge_branch,
                            Rcpp::Named("states") = states_out,
                            Rcpp::Named("duration") = DT.count());
}

template <typename datatype, typename RMATRIX>
Rcpp::List TRAISIE_calc_ll_cpp_local(const Rcpp::IntegerVector& ances,
                                     const RMATRIX& states,
                                     const Rcpp::NumericMatrix& forTime,
                                     const Rcpp::NumericVector& lambda_cs,
                                     const Rcpp::NumericVector& lambda_as,
                                     const Rcpp::NumericVector& mus,
                                     const Rcpp::NumericVector& gammas,
                                     const Rcpp::NumericMatrix& qs,
                                     const Rcpp::NumericMatrix& p,
                                     const Rcpp::NumericVector& trait_mainland_ancestor,
                                     const std::string& method,
                                     double atol,
                                     double rtol,
                                     bool see_states,
                                     bool use_normalization) {
  try {
    size_t num_unique_states = (states.ncol() - 1) / 3;

    return TRAISIE_calc_ll<TRAISIE::interval1<datatype>, datatype, RMATRIX>(std::make_unique<TRAISIE::interval1<datatype>>(lambda_cs,
                                                                lambda_as,
                                                                mus,
                                                                gammas,
                                                                qs,
                                                                p,
                                                                trait_mainland_ancestor,
                                                                num_unique_states),
                                                                ances,
                                                                states,
                                                                forTime,
                                                                method,
                                                                atol,
                                                                rtol,
                                                                see_states,
                                                                use_normalization);
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (const char* msg) {
    Rcpp::Rcout << msg << std::endl;
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}

RcppExport SEXP TRAISIE_calc_ll_cpp(SEXP ancesSEXP, SEXP statesSEXP, SEXP forTimeSEXP,
                                    SEXP lambda_csSEXP, SEXP lambda_asSEXP, SEXP musSEXP,
                                    SEXP gammasSEXP, SEXP qsSEXP,
                                    SEXP pSEXP, SEXP tmaSEXP,
                                    SEXP methodSEXP,
                                    SEXP atolSEXP, SEXP rtolSEXP,
                                    SEXP see_statesSEXP, SEXP use_normalizationSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ances(ancesSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type states(statesSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type forTime(forTimeSEXP);

  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda_cs(lambda_csSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda_as(lambda_asSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mus(musSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type gammas(gammasSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type qs(qsSEXP);

  Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type p(pSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tma(tmaSEXP);

  Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);

  Rcpp::traits::input_parameter< double >::type atol(atolSEXP);
  Rcpp::traits::input_parameter< double >::type rtol(rtolSEXP);

  Rcpp::traits::input_parameter< bool >::type see_states(see_statesSEXP);
  Rcpp::traits::input_parameter< bool >::type use_normalization(use_normalizationSEXP);

  rcpp_result_gen = Rcpp::wrap(TRAISIE_calc_ll_cpp_local<double, Rcpp::NumericMatrix>(ances, states, forTime,
                                                         lambda_cs, lambda_as, mus, gammas, qs,
                                                         p, tma,
                                                         method, atol, rtol,
                                                         see_states, use_normalization));
  return rcpp_result_gen;
  END_RCPP
}


// same as TRAISIE_calc_ll_cpp, but the tip states are complex numbers
// (used when missing species are given, see TRAISIE_n)
RcppExport SEXP TRAISIE_calc_ll_cpp_complex(SEXP ancesSEXP, SEXP statesSEXP, SEXP forTimeSEXP,
                                    SEXP lambda_csSEXP, SEXP lambda_asSEXP, SEXP musSEXP,
                                    SEXP gammasSEXP, SEXP qsSEXP,
                                    SEXP pSEXP, SEXP tmaSEXP,
                                    SEXP methodSEXP,
                                    SEXP atolSEXP, SEXP rtolSEXP,
                                    SEXP see_statesSEXP, SEXP use_normalizationSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ances(ancesSEXP);
  Rcpp::traits::input_parameter< Rcpp::ComplexMatrix >::type states(statesSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type forTime(forTimeSEXP);

  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda_cs(lambda_csSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type lambda_as(lambda_asSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mus(musSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type gammas(gammasSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type qs(qsSEXP);

  Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type p(pSEXP);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type tma(tmaSEXP);

  Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);

  Rcpp::traits::input_parameter< double >::type atol(atolSEXP);
  Rcpp::traits::input_parameter< double >::type rtol(rtolSEXP);

  Rcpp::traits::input_parameter< bool >::type see_states(see_statesSEXP);
  Rcpp::traits::input_parameter< bool >::type use_normalization(use_normalizationSEXP);

  rcpp_result_gen = Rcpp::wrap(TRAISIE_calc_ll_cpp_local<std::complex<double>, Rcpp::ComplexMatrix>(ances, states, forTime,
                                                         lambda_cs, lambda_as, mus, gammas, qs,
                                                         p, tma,
                                                         method, atol, rtol,
                                                         see_states, use_normalization));
  return rcpp_result_gen;
  END_RCPP
}
