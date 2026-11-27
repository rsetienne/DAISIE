//
//  Solver for the diversity-dependent DAISIE_DE equations.
//
//  In the diversity-dependent version every DAISIE_DE variable carries an
//  index n, the number of hidden species (alive, not in the data). The
//  equations are then linear: dy/dt = Q y, where the matrix Q is built in R
//  (R/DAISIE_DE_DD_intervals.R) and passed here as a list of its non-zero
//  entries (row, column, value). This file only integrates dy/dt = Q y with
//  the same odeint steppers as the other DAISIE_DE functions.
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <string>
#include <memory>
#include "config.h"    // NOLINT [build/include_subdir]
#include <Rcpp.h>
#include "DAISIE_DE_odeint.h"    // NOLINT [build/include_subdir]

namespace DAISIEDEDD {

// dy/dt = Q y, with Q stored row by row (only the non-zero entries)
struct linear_system {
  std::vector<size_t> row_start;   // entries of row r: row_start[r] .. row_start[r + 1] - 1
  std::vector<size_t> col;
  std::vector<double> val;

  linear_system(const std::vector<int>& i,
                const std::vector<int>& j,
                const std::vector<double>& x,
                size_t n_states) {
    row_start.assign(n_states + 1, 0);
    for (auto r : i) row_start[r + 1]++;
    for (size_t r = 0; r < n_states; ++r) row_start[r + 1] += row_start[r];
    col.resize(x.size());
    val.resize(x.size());
    auto next = row_start;
    for (size_t e = 0; e < x.size(); ++e) {
      auto pos = next[i[e]]++;
      col[pos] = j[e];
      val[pos] = x[e];
    }
  }

  size_t size() const noexcept { return row_start.size() - 1; }

  void operator()(const std::vector<double>& y,
                  std::vector<double>& dydt,
                  const double /* t */) const {
    for (size_t r = 0; r + 1 < row_start.size(); ++r) {
      double s = 0.0;
      for (size_t e = row_start[r]; e < row_start[r + 1]; ++e) {
        s += val[e] * y[col[e]];
      }
      dydt[r] = s;
    }
  }
};

}   // namespace DAISIEDEDD

// i, j: row and column (starting at 0) of the non-zero entries of Q, x: values
RcppExport SEXP DAISIE_DE_DD_cpp_solve(SEXP iSEXP, SEXP jSEXP, SEXP xSEXP,
                                       SEXP init_statesSEXP, SEXP timeSEXP,
                                       SEXP inte_methodSEXP,
                                       SEXP atolSEXP, SEXP rtolSEXP) {
  BEGIN_RCPP
  auto i = Rcpp::as<std::vector<int>>(iSEXP);
  auto j = Rcpp::as<std::vector<int>>(jSEXP);
  auto x = Rcpp::as<std::vector<double>>(xSEXP);
  auto states = Rcpp::as<std::vector<double>>(init_statesSEXP);
  auto time = Rcpp::as<std::vector<double>>(timeSEXP);
  auto inte_method = Rcpp::as<std::string>(inte_methodSEXP);
  auto atol = Rcpp::as<double>(atolSEXP);
  auto rtol = Rcpp::as<double>(rtolSEXP);

  auto od = std::make_unique<DAISIEDEDD::linear_system>(i, j, x, states.size());
  auto workhorse = Integrator<DAISIEDEDD::linear_system, double>(std::move(od),
                                                                 inte_method,
                                                                 atol, rtol);
  auto t0 = std::min(time[0], time[1]);
  auto t1 = std::max(time[0], time[1]);
  workhorse(states, t0, t1);
  return Rcpp::wrap(states);
  END_RCPP
}
