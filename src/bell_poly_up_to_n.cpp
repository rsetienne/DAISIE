#include <span>
#include "config.h"
#include "DAISIE_types.h"

using namespace Rcpp;


namespace detail {

  inline double bell_polynomials_up_to_n(int n, std::span<double> g_derives) {
    return 0.0;
  }
}
/*
bell_polynomials_up_to_n <- function(n, g_derivs) {
    B <- numeric(n + 1)
    B[1] <- 1  # B_0

    for (m in 1:n) {
      tmp <- numeric(0)
      for (k in 1:m) {
        tmp <- c(tmp, choose(m - 1, k - 1) * B[m - k + 1] * g_derivs[k])
      }
      o <- order(abs(tmp))
      B[m + 1] <- sum(tmp[o])
    }
    return(B)  # B[1] = B_0, ..., B[n+1] = B_n
  }
*/

RcppExport SEXP bell_polynomials_up_to_n(SEXP rn, SEXP rg_derives) {
  BEGIN_RCPP
  auto n = as<int>(rn);
  auto g_derives = as<std::vector<double>>(rg_derives);
  return wrap(detail::bell_polynomials_up_to_n(n, g_derives));
  END_RCPP
}
