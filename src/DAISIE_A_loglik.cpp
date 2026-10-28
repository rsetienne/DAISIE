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
#include <cmath>
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

// Finite carrying capacities (mean-field approximation). The rates are split
// into components: 0 (not diversity dependent), 1..n (cladogenesis on island
// i) and n+1..2n (colonisation of and dispersal into island i). The
// right-hand side is the sum of the right-hand sides of the components, each
// multiplied by a weight that depends on the block of equations (E, A, F or
// D = DE and DM): with nu_i = l_i + h_i(t) + p_i(t) (reconstructed lineages,
// unobserved species and, for a mainland species without clades, the
// presence of its non-endemic lineage), m_i = l_i + h_i(t) and
// o_i = nu_i - own_i (own_i: contribution of the branch being integrated to
// l_i),
//   cladogenesis on i: E: f(nu_i + 1), A: f(nu_i), F: f(m_i + 1), D: f(o_i + 1)
//   dispersal into i:  E: f(nu_i),     A: f(m_i),  F: f(m_i),     D: f(o_i)
// with f(x) = max(0, 1 - x / K_i): a lineage counts itself when it speciates
// on its own island and not when it disperses to another island. The
// recolonisation of the single island of the non-endemic lineage (which
// resets its colonisation time; a loss term of DM) is weighted with
// f(o_i + 1) instead of f(o_i), because that lineage is on island i.
// h and p are linearly interpolated on a grid of ages.
struct lik_ode_dd {
  std::vector<lik_ode*> comps;
  int n = 0;
  int G = 0;
  int S = 0;
  std::vector<double> l, K, gt, gh, gp;   // gh, gp: value at gt[g] in [i * G + g]
  std::vector<double> own;                // contribution of the branch to l
  mutable std::vector<double> tmp, w;     // w: 4 weights (E, A, F, D) per component
  mutable std::vector<double> wr;         // weight of the recolonisation reset

  double interp(const std::vector<double>& v, int i, int g, double t) const {
    if (G == 1) return v[i];
    const double t0 = gt[g], t1 = gt[g + 1];
    const double h0 = v[i * G + g], h1 = v[i * G + g + 1];
    double a = (t1 > t0) ? (t - t0) / (t1 - t0) : 1.0;
    a = std::min(1.0, std::max(0.0, a));
    return h0 + a * (h1 - h0);
  }

  // the grid may contain an age twice (a jump, e.g. at island submergence):
  // the interval with gt[g] < t <= gt[g + 1] is used, as with
  // findInterval(t, gt, left.open = TRUE) in R
  void weights(double t) const {
    const int nc = 2 * n + 1;
    w.assign(4 * nc, 1.0);
    wr.assign(n, 1.0);
    int g = 0;
    if (G > 1) {
      const int below = static_cast<int>(
        std::lower_bound(gt.begin(), gt.end(), t) - gt.begin());
      g = std::min(std::max(below - 1, 0), G - 2);
    }
    auto f = [](double x, double k) { return std::max(0.0, 1.0 - x / k); };
    for (int i = 0; i < n; ++i) {
      if (!std::isfinite(K[i])) continue;
      const double h = interp(gh, i, g, t);
      const double p = interp(gp, i, g, t);
      const double nu = l[i] + h + p;
      const double m = l[i] + h;
      const double o = nu - (own.empty() ? 0.0 : own[i]);
      double* wc = w.data() + 4 * (1 + i);       // cladogenesis on i
      wc[0] = f(nu + 1.0, K[i]);
      wc[1] = f(nu, K[i]);
      wc[2] = f(m + 1.0, K[i]);
      wc[3] = f(o + 1.0, K[i]);
      double* wd = w.data() + 4 * (1 + n + i);   // dispersal into i
      wd[0] = f(nu, K[i]);
      wd[1] = f(m, K[i]);
      wd[2] = f(m, K[i]);
      wd[3] = f(o, K[i]);
      wr[i] = f(o + 1.0, K[i]);
    }
  }

  void operator()(const std::vector<double>& x, std::vector<double>& dx,
                  double t) const {
    (*comps[0])(x, dx, t);
    weights(t);
    tmp.resize(x.size());
    const int nr = static_cast<int>(x.size()) / comps[0]->P;
    const int nc = static_cast<int>(comps.size());
    for (int c = 1; c < nc; ++c) {
      const double* wc = w.data() + 4 * c;
      if (wc[0] <= 0.0 && wc[1] <= 0.0 && wc[2] <= 0.0 && wc[3] <= 0.0) {
        continue;
      }
      (*comps[c])(x, tmp, t);
      for (size_t q = 0; q < x.size(); ++q) {
        const int r = static_cast<int>(q % nr);
        const int b = r < S ? 0 : (r == S ? 1 : (r <= 2 * S ? 2 : 3));
        dx[q] += wc[b] * tmp[q];
      }
    }
    // recolonisation reset (loss term of DM for single-island ranges)
    if (nr > 2 * S + 1) {
      const int P = comps[0]->P;
      for (int i = 0; i < n; ++i) {
        const lik_ode& st = *comps[1 + n + i];
        const double dw = wr[i] - w[4 * (1 + n + i) + 3];
        if (dw == 0.0 || st.col.empty() || st.col[i] == 0.0) continue;
        const int row = 3 * S + 1 + st.single[i];
        for (int k = 0; k < P; ++k) {
          dx[row + nr * k] -= dw * st.col[i] * x[row + nr * k];
        }
      }
    }
  }
};

// Forward-in-time equations of the mean-field approximation, integrated
// over one grid interval at a time (age from 'hi' down to 'lo'; u = -age).
//   mode 0: G_R = h_R / E_R, the expected numbers h_R of unobserved species
//           with range R divided by E_R (this removes the singular rates of
//           the conditioned branching process):
//           dG_R/du = -omega_R G_R + sum_R' Qe[R', R] G_R'
//                     + sum_events rho E_b G_p (daughter a = R) + ...
//                     + sources / E_R;
//   mode 1: U, the solution of the adjoint of the equations for (A, F),
//           giving the posterior state of the non-endemic lineage of a
//           mainland species without clades.
struct dd_forward_ode {
  const lik_ode_dd* dd = nullptr;
  int mode = 0;
  int S = 0;
  double lo = 0.0, hi = 0.0;
  std::vector<double> Elo, Ehi;   // E at lo and hi
  std::vector<double> Qt;         // reconstructed lineages per range
  bool hidden_col = false;
  std::vector<double> slo, shi;   // sources (divided by E) at lo and hi
  mutable std::vector<double> E;

  void operator()(const std::vector<double>& x, std::vector<double>& dx,
                  double u) const {
    const double t = -u;
    double a = (hi > lo) ? (t - lo) / (hi - lo) : 0.0;
    a = std::min(1.0, std::max(0.0, a));
    E.resize(S);
    for (int s = 0; s < S; ++s) E[s] = Elo[s] + a * (Ehi[s] - Elo[s]);
    dd->weights(t);
    std::fill(dx.begin(), dx.end(), 0.0);
    const int nc = static_cast<int>(dd->comps.size());
    for (int c = 0; c < nc; ++c) {
      const lik_ode& st = *dd->comps[c];
      const double* wc = dd->w.data() + 4 * c;
      const double wE = wc[0], wA = wc[1], wF = wc[2], wD = wc[3];
      if (mode == 0) {
        if (wE > 0.0) {
          for (int s = 0; s < S; ++s) dx[s] -= wE * st.omega_E[s] * x[s];
          for (size_t q = 0; q < st.qe_v.size(); ++q) {
            dx[st.qe_j[q]] += wE * st.qe_v[q] * x[st.qe_i[q]];
          }
          for (size_t e = 0; e < st.tri_r.size(); ++e) {
            const double r = wE * st.tri_r[e] * x[st.tri_p[e]];
            dx[st.tri_a[e]] += r * E[st.tri_b[e]];
            dx[st.tri_b[e]] += r * E[st.tri_a[e]];
          }
        }
        if (wD > 0.0) {
          for (size_t e = 0; e < st.tri_r.size(); ++e) {
            const double r = wD * st.tri_r[e] * Qt[st.tri_p[e]];
            dx[st.tri_a[e]] += r;
            dx[st.tri_b[e]] += r;
          }
        }
        if (hidden_col && wA > 0.0) {
          for (size_t i = 0; i < st.col.size(); ++i) {
            dx[st.single[i]] += wA * st.col[i];
          }
        }
      } else {
        // x = (U_A, U_F[0..S-1])
        if (wA > 0.0) {
          for (size_t i = 0; i < st.col.size(); ++i) {
            dx[0] -= wA * st.col[i] * x[0];
            dx[1 + st.single[i]] += wA * st.col[i] * x[0];
          }
        }
        if (wF > 0.0) {
          for (int s = 0; s < S; ++s) {
            dx[1 + s] -= wF * st.omega_F[s] * x[1 + s];
            dx[0] += wF * (st.g[s] + st.laa * E[s]) * x[1 + s];
          }
          for (size_t q = 0; q < st.qm_v.size(); ++q) {
            dx[1 + st.qm_j[q]] += wF * st.qm_v[q] * x[1 + st.qm_i[q]];
          }
          for (size_t e = 0; e < st.end_r.size(); ++e) {
            dx[0] += wF * st.end_r[e] * E[st.end_a[e]] * E[st.end_b[e]] *
              x[1 + st.end_p[e]];
          }
          for (size_t e = 0; e < st.bud_r.size(); ++e) {
            dx[1 + st.bud_y[e]] += wF * st.bud_r[e] * E[st.bud_x[e]] *
              x[1 + st.bud_p[e]];
          }
        }
      }
    }
    if (mode == 0 && !slo.empty()) {
      for (int s = 0; s < S; ++s) dx[s] += slo[s] + a * (shi[s] - slo[s]);
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

// As DAISIE_A_cpp_integrate, for finite carrying capacities: 'ptrs' is a list
// of 2n + 1 external pointers (the rate components), 'l' the contribution of
// the reconstructed lineages to nu (per island), 'K' the carrying
// capacities, and 'gt' (increasing ages), 'gh' and 'gp' (n x length(gt),
// column-major as in R) the grids of the other contributions.
RcppExport SEXP DAISIE_A_cpp_integrate_dd(SEXP ptrsSEXP, SEXP ySEXP,
                                          SEXP paSEXP, SEXP pbSEXP,
                                          SEXP pcSEXP, SEXP PSEXP,
                                          SEXP t0SEXP, SEXP t1SEXP,
                                          SEXP methodSEXP, SEXP atolSEXP,
                                          SEXP rtolSEXP, SEXP lSEXP,
                                          SEXP KSEXP, SEXP gtSEXP,
                                          SEXP ghSEXP, SEXP gpSEXP,
                                          SEXP ownSEXP) {
  BEGIN_RCPP
  Rcpp::List ptrs(ptrsSEXP);
  DAISIEA::lik_ode_dd od;
  const int P = Rcpp::as<int>(PSEXP);
  od.own = Rcpp::as<std::vector<double>>(ownSEXP);
  const auto pa = Rcpp::as<std::vector<int>>(paSEXP);
  const auto pb = Rcpp::as<std::vector<int>>(pbSEXP);
  const auto pc = Rcpp::as<std::vector<int>>(pcSEXP);
  for (R_xlen_t c = 0; c < ptrs.size(); ++c) {
    Rcpp::XPtr<DAISIEA::lik_ode> p(Rcpp::as<SEXP>(ptrs[c]));
    p->P = P;
    p->pa = pa;
    p->pb = pb;
    p->pc = pc;
    od.comps.push_back(p.get());
  }
  od.S = od.comps[0]->S;
  od.n = (static_cast<int>(ptrs.size()) - 1) / 2;
  od.l = Rcpp::as<std::vector<double>>(lSEXP);
  od.K = Rcpp::as<std::vector<double>>(KSEXP);
  od.gt = Rcpp::as<std::vector<double>>(gtSEXP);
  const auto ghm = Rcpp::as<std::vector<double>>(ghSEXP);
  const auto gpm = Rcpp::as<std::vector<double>>(gpSEXP);
  od.G = static_cast<int>(od.gt.size());
  if (2 * od.n + 1 != static_cast<int>(ptrs.size()) ||
      static_cast<int>(od.l.size()) != od.n ||
      static_cast<int>(od.K.size()) != od.n ||
      static_cast<int>(ghm.size()) != od.n * od.G ||
      static_cast<int>(gpm.size()) != od.n * od.G) {
    throw std::invalid_argument("inconsistent diversity-dependence input");
  }
  // R stores the n x G matrices column-major: element (i, g) at i + n * g
  od.gh.assign(static_cast<size_t>(od.n) * od.G, 0.0);
  od.gp.assign(static_cast<size_t>(od.n) * od.G, 0.0);
  for (int i = 0; i < od.n; ++i) {
    for (int g = 0; g < od.G; ++g) {
      od.gh[i * od.G + g] = ghm[i + od.n * g];
      od.gp[i * od.G + g] = gpm[i + od.n * g];
    }
  }
  auto y = Rcpp::as<std::vector<double>>(ySEXP);
  if (y.size() % P != 0) {
    throw std::invalid_argument("state length is not a multiple of P");
  }
  const double t0 = Rcpp::as<double>(t0SEXP);
  const double t1 = Rcpp::as<double>(t1SEXP);
  const std::string method = Rcpp::as<std::string>(methodSEXP);
  const double atol = Rcpp::as<double>(atolSEXP);
  const double rtol = Rcpp::as<double>(rtolSEXP);
  if (t1 > t0) {
    odeintcpp::integrate(method, &od, &y, t0, t1,
                         SECSSE_DEFAULT_DTF * (t1 - t0), atol, rtol);
  }
  return Rcpp::wrap(y);
  END_RCPP
}

// Forward-in-time pass of the mean-field approximation over grid ages
// gp[0] < ... < gp[m-1] (no submergence in between): starting from x0 at
// age gp[m-1], returns the state at every grid age (S' x m, column-major).
// E: S x m; lseg: n x (m-1) (reconstructed lineages per interval); Qt:
// S x (m-1); hcol: m-1 flags (earlier colonisations without sampled
// descendants); src: S x m or empty (sources divided by E).
RcppExport SEXP DAISIE_A_cpp_dd_forward(SEXP ptrsSEXP, SEXP modeSEXP,
                                        SEXP gpSEXP, SEXP ESEXP,
                                        SEXP lsegSEXP, SEXP KSEXP,
                                        SEXP gtSEXP, SEXP ghSEXP,
                                        SEXP gppSEXP, SEXP QtSEXP,
                                        SEXP hcolSEXP, SEXP srcSEXP,
                                        SEXP x0SEXP, SEXP methodSEXP,
                                        SEXP atolSEXP, SEXP rtolSEXP) {
  BEGIN_RCPP
  Rcpp::List ptrs(ptrsSEXP);
  DAISIEA::lik_ode_dd od;
  for (R_xlen_t c = 0; c < ptrs.size(); ++c) {
    Rcpp::XPtr<DAISIEA::lik_ode> p(Rcpp::as<SEXP>(ptrs[c]));
    od.comps.push_back(p.get());
  }
  od.S = od.comps[0]->S;
  od.n = (static_cast<int>(ptrs.size()) - 1) / 2;
  od.K = Rcpp::as<std::vector<double>>(KSEXP);
  od.gt = Rcpp::as<std::vector<double>>(gtSEXP);
  od.G = static_cast<int>(od.gt.size());
  const auto ghm = Rcpp::as<std::vector<double>>(ghSEXP);
  const auto gpm = Rcpp::as<std::vector<double>>(gppSEXP);
  od.gh.assign(static_cast<size_t>(od.n) * od.G, 0.0);
  od.gp.assign(static_cast<size_t>(od.n) * od.G, 0.0);
  for (int i = 0; i < od.n; ++i) {
    for (int g = 0; g < od.G; ++g) {
      od.gh[i * od.G + g] = ghm[i + od.n * g];
      od.gp[i * od.G + g] = gpm[i + od.n * g];
    }
  }
  // unobserved daughters of reconstructed lineages: the cladogenesis of a
  // reconstructed lineage on island i is weighted with f(nu_i), i.e. its
  // contribution to l_i is taken to be its presence on i
  od.own.assign(od.n, 1.0);
  const int mode = Rcpp::as<int>(modeSEXP);
  const auto gp = Rcpp::as<std::vector<double>>(gpSEXP);
  const auto Em = Rcpp::as<std::vector<double>>(ESEXP);
  const auto lseg = Rcpp::as<std::vector<double>>(lsegSEXP);
  const auto Qtm = Rcpp::as<std::vector<double>>(QtSEXP);
  const auto hcol = Rcpp::as<std::vector<int>>(hcolSEXP);
  const auto src = Rcpp::as<std::vector<double>>(srcSEXP);
  auto x = Rcpp::as<std::vector<double>>(x0SEXP);
  const std::string method = Rcpp::as<std::string>(methodSEXP);
  const double atol = Rcpp::as<double>(atolSEXP);
  const double rtol = Rcpp::as<double>(rtolSEXP);
  const int m = static_cast<int>(gp.size());
  const int S = od.S;
  const int dim = static_cast<int>(x.size());
  Rcpp::NumericMatrix out(dim, m);
  for (int k = 0; k < dim; ++k) out(k, m - 1) = x[k];
  DAISIEA::dd_forward_ode fo;
  fo.dd = &od;
  fo.mode = mode;
  fo.S = S;
  od.l.assign(od.n, 0.0);
  for (int g = m - 1; g > 0; --g) {
    fo.lo = gp[g - 1];
    fo.hi = gp[g];
    fo.Elo.assign(Em.begin() + S * (g - 1), Em.begin() + S * g);
    fo.Ehi.assign(Em.begin() + S * g, Em.begin() + S * (g + 1));
    fo.Qt.assign(Qtm.begin() + S * (g - 1), Qtm.begin() + S * g);
    fo.hidden_col = hcol[g - 1] != 0;
    for (int i = 0; i < od.n; ++i) od.l[i] = lseg[i + od.n * (g - 1)];
    if (src.empty()) {
      fo.slo.clear();
      fo.shi.clear();
    } else {
      fo.slo.assign(src.begin() + S * (g - 1), src.begin() + S * g);
      fo.shi.assign(src.begin() + S * g, src.begin() + S * (g + 1));
    }
    if (fo.hi > fo.lo) {
      const double u0 = -fo.hi, u1 = -fo.lo;
      odeintcpp::integrate(method, &fo, &x, u0, u1,
                           SECSSE_DEFAULT_DTF * (u1 - u0), atol, rtol);
    }
    for (int k = 0; k < dim; ++k) out(k, g - 1) = x[k];
  }
  return out;
  END_RCPP
}

// Integrates the likelihood equations (finite carrying capacities) through
// the ages gp[0] < ... < gp[m-1] (no breaks in between other than changes of
// the reconstructed lineages, given per interval in 'lseg', n x (m-1)), from
// the state y at gp[0]; returns the state at every age (nrow x m).
RcppExport SEXP DAISIE_A_cpp_integrate_dd_grid(SEXP ptrsSEXP, SEXP ySEXP,
                                               SEXP PSEXP, SEXP gpSEXP,
                                               SEXP lsegSEXP,
                                               SEXP methodSEXP,
                                               SEXP atolSEXP, SEXP rtolSEXP,
                                               SEXP KSEXP, SEXP gtSEXP,
                                               SEXP ghSEXP, SEXP gppSEXP) {
  BEGIN_RCPP
  Rcpp::List ptrs(ptrsSEXP);
  DAISIEA::lik_ode_dd od;
  const int P = Rcpp::as<int>(PSEXP);
  for (R_xlen_t c = 0; c < ptrs.size(); ++c) {
    Rcpp::XPtr<DAISIEA::lik_ode> p(Rcpp::as<SEXP>(ptrs[c]));
    p->P = P;
    if (P == 1) {
      p->pa.assign(1, 0);
      p->pb.assign(1, 0);
      p->pc.assign(1, 0);
    }
    od.comps.push_back(p.get());
  }
  if (P != 1) throw std::invalid_argument("only one coefficient supported");
  od.S = od.comps[0]->S;
  od.n = (static_cast<int>(ptrs.size()) - 1) / 2;
  od.K = Rcpp::as<std::vector<double>>(KSEXP);
  od.gt = Rcpp::as<std::vector<double>>(gtSEXP);
  od.G = static_cast<int>(od.gt.size());
  const auto ghm = Rcpp::as<std::vector<double>>(ghSEXP);
  const auto gpm = Rcpp::as<std::vector<double>>(gppSEXP);
  od.gh.assign(static_cast<size_t>(od.n) * od.G, 0.0);
  od.gp.assign(static_cast<size_t>(od.n) * od.G, 0.0);
  for (int i = 0; i < od.n; ++i) {
    for (int g = 0; g < od.G; ++g) {
      od.gh[i * od.G + g] = ghm[i + od.n * g];
      od.gp[i * od.G + g] = gpm[i + od.n * g];
    }
  }
  const auto gp = Rcpp::as<std::vector<double>>(gpSEXP);
  const auto lseg = Rcpp::as<std::vector<double>>(lsegSEXP);
  auto y = Rcpp::as<std::vector<double>>(ySEXP);
  const std::string method = Rcpp::as<std::string>(methodSEXP);
  const double atol = Rcpp::as<double>(atolSEXP);
  const double rtol = Rcpp::as<double>(rtolSEXP);
  const int m = static_cast<int>(gp.size());
  const int dim = static_cast<int>(y.size());
  Rcpp::NumericMatrix out(dim, m);
  for (int k = 0; k < dim; ++k) out(k, 0) = y[k];
  od.l.assign(od.n, 0.0);
  for (int g = 1; g < m; ++g) {
    for (int i = 0; i < od.n; ++i) od.l[i] = lseg[i + od.n * (g - 1)];
    const double t0 = gp[g - 1], t1 = gp[g];
    if (t1 > t0) {
      odeintcpp::integrate(method, &od, &y, t0, t1,
                           SECSSE_DEFAULT_DTF * (t1 - t0), atol, rtol);
    }
    for (int k = 0; k < dim; ++k) out(k, g) = y[k];
  }
  return out;
  END_RCPP
}
