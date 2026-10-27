//
//  C++ (odeint) integration of the likelihood equations of the archipelago
//  model DAISIE_A for infinite carrying capacities.
//
//  As for DAISIE_DE, the tree traversal stays in R and the integration of
//  the equations along branches is done here. The rate structures of a
//  given set of existing islands are passed once and kept in an external
//  pointer; the power-series algebra (for missing species) is passed with
//  every integration.
//
//  State: a matrix with rows E (S), A (1), F (S) and optionally DE (S) and
//  DM (S), and one column per power-series coefficient, stored column-major
//  as in R: element (row, k) at row + nrow * k.
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <string>
#include <memory>
#include <algorithm>
#include <stdexcept>
#include "config.h"    // NOLINT [build/include_subdir]
#include <Rcpp.h>
#include "DAISIE_DE_odeint.h"    // NOLINT [build/include_subdir]

namespace DAISIEA {

struct lik_ode {
  int S = 0;          // number of range states
  int P = 1;          // number of power-series coefficients

  // sparse transition matrices (0-based): d[i] += v * x[j]
  std::vector<int> qe_i, qe_j, qm_i, qm_j;
  std::vector<double> qe_v, qm_v;

  std::vector<double> omega_E, omega_F, omega_M, g;
  std::vector<double> col;      // colonisation rate per island
  std::vector<int> single;      // range index of every single-island range
  double sumcol = 0.0;
  double laa = 0.0;

  // cladogenetic events of endemics: parent, daughters a and b, rate
  std::vector<int> tri_p, tri_a, tri_b;
  std::vector<double> tri_r;
  // events in which the non-endemic survives: parent, parent afterwards,
  // new endemic, rate
  std::vector<int> bud_p, bud_y, bud_x;
  std::vector<double> bud_r;
  // events in which the non-endemic is replaced by two endemics
  std::vector<int> end_p, end_a, end_b;
  std::vector<double> end_r;

  // truncated power-series product: out[pc] += u[pa] * v[pb]
  std::vector<int> pa, pb, pc;

  // work buffers (S x P, contiguous per range)
  mutable std::vector<double> endE, endD;

  // u and v strided (stride su, sv), output strided (so)
  inline void conv(const double* u, int su, const double* v, int sv,
                   double r, double* out, int so) const {
    const size_t np = pa.size();
    for (size_t q = 0; q < np; ++q) {
      out[pc[q] * so] += r * u[pa[q] * su] * v[pb[q] * sv];
    }
  }

  void operator()(const std::vector<double>& x, std::vector<double>& dx,
                  double /* t */) const {
    const int nr = static_cast<int>(x.size()) / P;
    const bool has_D = nr > 2 * S + 1;
    const int iA = S;
    const int oF = S + 1;
    const int oDE = 2 * S + 1;
    const int oDM = 3 * S + 1;
    const double* X = x.data();
    double* D = dx.data();

    // loss terms
    for (int k = 0; k < P; ++k) {
      const int c = nr * k;
      for (int s = 0; s < S; ++s) {
        D[c + s] = -omega_E[s] * X[c + s];
        D[c + oF + s] = -omega_F[s] * X[c + oF + s];
        if (has_D) {
          D[c + oDE + s] = -omega_E[s] * X[c + oDE + s];
          D[c + oDM + s] = -omega_M[s] * X[c + oDM + s];
        }
      }
      // mainland species absent
      double dA = -sumcol * X[c + iA];
      for (size_t i = 0; i < col.size(); ++i) {
        dA += col[i] * X[c + oF + single[i]];
      }
      D[c + iA] = dA;
    }
    // range changes (dispersal, local extinction, recolonisation)
    for (size_t q = 0; q < qe_v.size(); ++q) {
      const int i = qe_i[q], j = qe_j[q];
      const double v = qe_v[q];
      for (int k = 0; k < P; ++k) {
        const int c = nr * k;
        D[c + i] += v * X[c + j];
        if (has_D) D[c + oDE + i] += v * X[c + oDE + j];
      }
    }
    for (size_t q = 0; q < qm_v.size(); ++q) {
      const int i = qm_i[q], j = qm_j[q];
      const double v = qm_v[q];
      for (int k = 0; k < P; ++k) {
        const int c = nr * k;
        D[c + oF + i] += v * X[c + oF + j];
        if (has_D) D[c + oDM + i] += v * X[c + oDM + j];
      }
    }
    // global extinction
    for (int s = 0; s < S; ++s) {
      D[s] += g[s];   // constant series: coefficient 0 only
      for (int k = 0; k < P; ++k) {
        const int c = nr * k;
        D[c + oF + s] += g[s] * X[c + iA];
      }
    }
    // cladogenesis of endemics
    for (size_t e = 0; e < tri_r.size(); ++e) {
      const int p = tri_p[e], a = tri_a[e], b = tri_b[e];
      const double r = tri_r[e];
      conv(X + a, nr, X + b, nr, r, D + p, nr);
      if (has_D) {
        conv(X + oDE + a, nr, X + b, nr, r, D + oDE + p, nr);
        conv(X + a, nr, X + oDE + b, nr, r, D + oDE + p, nr);
      }
    }
    // events in which the non-endemic survives (budding, island-level
    // anagenesis)
    for (size_t e = 0; e < bud_r.size(); ++e) {
      const int p = bud_p[e], y = bud_y[e], xx = bud_x[e];
      const double r = bud_r[e];
      conv(X + oF + y, nr, X + xx, nr, r, D + oF + p, nr);
      if (has_D) {
        conv(X + oDM + y, nr, X + xx, nr, r, D + oDM + p, nr);
        conv(X + oF + y, nr, X + oDE + xx, nr, r, D + oDM + p, nr);
      }
    }
    // events after which the mainland species is absent (anagenesis,
    // single-island split, vicariance of the non-endemic)
    endE.assign(static_cast<size_t>(S) * P, 0.0);
    for (int s = 0; s < S; ++s) {
      for (int k = 0; k < P; ++k) endE[s * P + k] = laa * X[nr * k + s];
    }
    for (size_t e = 0; e < end_r.size(); ++e) {
      conv(X + end_a[e], nr, X + end_b[e], nr, end_r[e],
           endE.data() + end_p[e] * P, 1);
    }
    for (int s = 0; s < S; ++s) {
      conv(endE.data() + s * P, 1, X + iA, nr, 1.0, D + oF + s, nr);
    }
    if (has_D) {
      endD.assign(static_cast<size_t>(S) * P, 0.0);
      for (int s = 0; s < S; ++s) {
        for (int k = 0; k < P; ++k) endD[s * P + k] = laa * X[nr * k + oDE + s];
      }
      for (size_t e = 0; e < end_r.size(); ++e) {
        conv(X + oDE + end_a[e], nr, X + end_b[e], nr, end_r[e],
             endD.data() + end_p[e] * P, 1);
        conv(X + end_a[e], nr, X + oDE + end_b[e], nr, end_r[e],
             endD.data() + end_p[e] * P, 1);
      }
      for (int s = 0; s < S; ++s) {
        conv(endD.data() + s * P, 1, X + iA, nr, 1.0, D + oDM + s, nr);
      }
    }
  }
};

template <typename T>
std::vector<T> get_vec(const Rcpp::List& l, const char* name) {
  return Rcpp::as<std::vector<T>>(l[name]);
}

}   // namespace DAISIEA

// Create the rate structures of one set of existing islands. All indices
// in 'st' are 0-based.
RcppExport SEXP DAISIE_A_cpp_make(SEXP stSEXP) {
  BEGIN_RCPP
  Rcpp::List st(stSEXP);
  auto od = new DAISIEA::lik_ode();
  od->S = Rcpp::as<int>(st["S"]);
  od->qe_i = DAISIEA::get_vec<int>(st, "qe_i");
  od->qe_j = DAISIEA::get_vec<int>(st, "qe_j");
  od->qe_v = DAISIEA::get_vec<double>(st, "qe_v");
  od->qm_i = DAISIEA::get_vec<int>(st, "qm_i");
  od->qm_j = DAISIEA::get_vec<int>(st, "qm_j");
  od->qm_v = DAISIEA::get_vec<double>(st, "qm_v");
  od->omega_E = DAISIEA::get_vec<double>(st, "omega_E");
  od->omega_F = DAISIEA::get_vec<double>(st, "omega_F");
  od->omega_M = DAISIEA::get_vec<double>(st, "omega_M");
  od->g = DAISIEA::get_vec<double>(st, "g");
  od->col = DAISIEA::get_vec<double>(st, "col");
  od->single = DAISIEA::get_vec<int>(st, "single");
  od->sumcol = 0.0;
  for (double c : od->col) od->sumcol += c;
  od->laa = Rcpp::as<double>(st["laa"]);
  od->tri_p = DAISIEA::get_vec<int>(st, "tri_p");
  od->tri_a = DAISIEA::get_vec<int>(st, "tri_a");
  od->tri_b = DAISIEA::get_vec<int>(st, "tri_b");
  od->tri_r = DAISIEA::get_vec<double>(st, "tri_r");
  od->bud_p = DAISIEA::get_vec<int>(st, "bud_from");
  od->bud_y = DAISIEA::get_vec<int>(st, "bud_p");
  od->bud_x = DAISIEA::get_vec<int>(st, "bud_x");
  od->bud_r = DAISIEA::get_vec<double>(st, "bud_r");
  od->end_p = DAISIEA::get_vec<int>(st, "end_p");
  od->end_a = DAISIEA::get_vec<int>(st, "end_a");
  od->end_b = DAISIEA::get_vec<int>(st, "end_b");
  od->end_r = DAISIEA::get_vec<double>(st, "end_r");
  Rcpp::XPtr<DAISIEA::lik_ode> ptr(od, true);
  return ptr;
  END_RCPP
}

// Is the external pointer empty (e.g. after saving and reloading)?
RcppExport SEXP DAISIE_A_cpp_ptr_is_null(SEXP ptrSEXP) {
  BEGIN_RCPP
  return Rcpp::wrap(TYPEOF(ptrSEXP) != EXTPTRSXP ||
                    R_ExternalPtrAddr(ptrSEXP) == nullptr);
  END_RCPP
}

// Integrate the state y (vector of a nrow x P matrix) from age t0 to t1.
RcppExport SEXP DAISIE_A_cpp_integrate(SEXP ptrSEXP, SEXP ySEXP,
                                       SEXP paSEXP, SEXP pbSEXP, SEXP pcSEXP,
                                       SEXP PSEXP, SEXP t0SEXP, SEXP t1SEXP,
                                       SEXP methodSEXP, SEXP atolSEXP,
                                       SEXP rtolSEXP) {
  BEGIN_RCPP
  Rcpp::XPtr<DAISIEA::lik_ode> od(ptrSEXP);
  od->P = Rcpp::as<int>(PSEXP);
  od->pa = Rcpp::as<std::vector<int>>(paSEXP);
  od->pb = Rcpp::as<std::vector<int>>(pbSEXP);
  od->pc = Rcpp::as<std::vector<int>>(pcSEXP);
  auto y = Rcpp::as<std::vector<double>>(ySEXP);
  if (y.size() % od->P != 0) {
    throw std::invalid_argument("state length is not a multiple of P");
  }
  const double t0 = Rcpp::as<double>(t0SEXP);
  const double t1 = Rcpp::as<double>(t1SEXP);
  const std::string method = Rcpp::as<std::string>(methodSEXP);
  const double atol = Rcpp::as<double>(atolSEXP);
  const double rtol = Rcpp::as<double>(rtolSEXP);
  if (t1 > t0) {
    odeintcpp::integrate(method, od.get(), &y, t0, t1,
                         SECSSE_DEFAULT_DTF * (t1 - t0), atol, rtol);
  }
  return Rcpp::wrap(y);
  END_RCPP
}
