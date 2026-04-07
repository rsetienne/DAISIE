//
//  Copyright (c) 2026 Thijs Janzen
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#include "DAISIE_DE_util.h"    // NOLINT [build/include_subdir]

double DAISIE_DE_general_cpp(const Rcpp::NumericVector& brts,
                             size_t missnumspec,
                             double lambda_c,
                             double lambda_a,
                             double mu,
                             double gamma,
                             size_t stac,
                             std::string method,
                             double rtol,
                             double atol,
                             std::string type) {
  solver solution(brts, missnumspec, lambda_c, lambda_a, mu, gamma, method, atol, rtol);
  return solution.calculate_likelihood(type, stac);
}

//' Wrapper for the DAISIE_DE general integrator
 //'
 //' @description This is the rcpp function to do single branch DAISIE_DE calculations
 //' @name DAISIE_DE_logpEC_general_rcpp
 //' @export DAISIE_DE_logpEC_general_rcpp
 //' @return list
 RcppExport SEXP DAISIE_DE_general_cpp(SEXP brtsSEXP, SEXP missnumspecSEXP,
                                               SEXP lambda_cSEXP, SEXP lambda_aSEXP, SEXP muSEXP, SEXP gammaSEXP,
                                               SEXP stacSEXP,
                                               SEXP inte_methodSEXP,
                                               SEXP atolSEXP, SEXP rtolSEXP,
                                               SEXP typeSEXP) {
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

   Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);

   rcpp_result_gen = Rcpp::wrap(DAISIE_DE_general_cpp(brts, missnumspec, lambda_c, lambda_a, mu, gamma, stac, inte_method, atol, rtol, type));
   return rcpp_result_gen;
   END_RCPP
 }
