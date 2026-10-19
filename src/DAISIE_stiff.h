//
//  Stiff integration with a sparse-Jacobian Rosenbrock method.
//
//  boost::odeint::rosenbrock4 factorizes a dense ublas Jacobian in every step,
//  which costs O(n^3) and makes it 20-70x slower than lsodes for the DAISIE
//  systems (n = 2 * lx + 1 and more). The DAISIE systems have a very sparse
//  Jacobian (a handful of non-zeros per row). This header implements the same
//  Rosenbrock method (same coefficients, same error estimate and the same step
//  size controller as boost::odeint::rosenbrock4 and rosenbrock4_controller),
//  but
//    * the sparsity pattern of the Jacobian is detected once by probing,
//    * the unknowns are reordered with reverse Cuthill-McKee,
//    * the Jacobian is evaluated by grouped probing (one right-hand side
//      evaluation per band diagonal instead of one per column) and
//    * the linear systems are solved with a banded LU factorization with
//      partial pivoting.
//
//  The Rosenbrock coefficients below are those of
//  boost/numeric/odeint/stepper/rosenbrock4.hpp
//  (Copyright 2011-2012 Karsten Ahnert, Copyright 2011-2012 Mario Mulansky,
//  Boost Software License, Version 1.0).
//
//  Systems are described by a 'system' type providing
//    void eval(const state& x, state& dx, double t) const   (non counting)
//    void operator()(const state& x, state& dx, double t) const (may count)
//  and the compile time flags
//    static constexpr bool stiff_linear;   dx/dt = A(t) x
//    static constexpr bool stiff_time_dependent;  A or f depends on t
//

#ifndef DAISIE_STIFF_H_INCLUDED
#define DAISIE_STIFF_H_INCLUDED

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <string>
#include <vector>


namespace daisie_stiff {

  // ---------------------------------------------------------------------------
  // Banded LU with partial pivoting (LAPACK dgbtf2 / dgbtrs, LAPACK band layout)
  // ---------------------------------------------------------------------------
  class band_lu
  {
  public:
    void init(int n, int kl, int ku)
    {
      n_ = n; kl_ = kl; ku_ = ku;
      ldab_ = 2 * kl + ku + 1;
      ab_.assign(static_cast<size_t>(ldab_) * n, 0.0);
      piv_.assign(n, 0);
    }

    int n() const { return n_; }

    // element (i, j) of the matrix; valid within the band before factorization
    double& at(int i, int j)
    {
      return ab_[static_cast<size_t>(kl_ + ku_ + i - j) + static_cast<size_t>(j) * ldab_];
    }

    void zero() { std::fill(ab_.begin(), ab_.end(), 0.0); }

    void factorize()
    {
      int ju = 0;
      for (int j = 0; j < n_; ++j) {
        const int km = std::min(kl_, n_ - 1 - j);
        int p = 0;
        double amax = std::fabs(at(j, j));
        for (int i = 1; i <= km; ++i) {
          const double v = std::fabs(at(j + i, j));
          if (v > amax) { amax = v; p = i; }
        }
        piv_[j] = j + p;
        if (!(amax > 0.0)) {
          throw std::runtime_error("daisie_stiff: singular matrix in banded LU");
        }
        ju = std::max(ju, std::min(j + ku_ + p, n_ - 1));
        if (p != 0) {
          for (int c = j; c <= ju; ++c) std::swap(at(j, c), at(j + p, c));
        }
        if (km > 0) {
          const double inv = 1.0 / at(j, j);
          for (int i = 1; i <= km; ++i) at(j + i, j) *= inv;
          for (int c = j + 1; c <= ju; ++c) {
            const double ajc = at(j, c);
            if (ajc != 0.0) {
              for (int i = 1; i <= km; ++i) at(j + i, c) -= at(j + i, j) * ajc;
            }
          }
        }
      }
    }

    // solve in place
    template <typename Vec>
    void solve(Vec& b)
    {
      for (int j = 0; j < n_; ++j) {
        const int l = piv_[j];
        if (l != j) std::swap(b[l], b[j]);
        const int km = std::min(kl_, n_ - 1 - j);
        const double bj = b[j];
        if (bj != 0.0) {
          for (int i = 1; i <= km; ++i) b[j + i] -= at(j + i, j) * bj;
        }
      }
      const int kv = kl_ + ku_;
      for (int j = n_ - 1; j >= 0; --j) {
        b[j] /= at(j, j);
        const double bj = b[j];
        if (bj != 0.0) {
          for (int i = std::max(0, j - kv); i < j; ++i) b[i] -= at(i, j) * bj;
        }
      }
    }

  private:
    int n_ = 0, kl_ = 0, ku_ = 0, ldab_ = 0;
    std::vector<double> ab_;
    std::vector<int> piv_;
  };


  // ---------------------------------------------------------------------------
  // reverse Cuthill-McKee ordering of the symmetrized pattern.
  // Returns perm with perm[new] = old.
  // ---------------------------------------------------------------------------
  inline std::vector<int> rcm_order(int n, const std::vector<std::pair<int, int>>& nz)
  {
    std::vector<std::vector<int>> adj(n);
    for (const auto& e : nz) {
      if (e.first != e.second) {
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
      }
    }
    for (auto& a : adj) {
      std::sort(a.begin(), a.end());
      a.erase(std::unique(a.begin(), a.end()), a.end());
    }
    std::vector<int> order;
    order.reserve(n);
    std::vector<char> seen(n, 0);
    std::vector<int> byDeg(n);
    std::iota(byDeg.begin(), byDeg.end(), 0);
    std::stable_sort(byDeg.begin(), byDeg.end(),
                     [&](int a, int b) { return adj[a].size() < adj[b].size(); });
    for (int s : byDeg) {
      if (seen[s]) continue;
      size_t head = order.size();
      order.push_back(s);
      seen[s] = 1;
      while (head < order.size()) {
        const int v = order[head++];
        std::vector<int> nb;
        for (int w : adj[v]) if (!seen[w]) { seen[w] = 1; nb.push_back(w); }
        std::stable_sort(nb.begin(), nb.end(),
                         [&](int a, int b) { return adj[a].size() < adj[b].size(); });
        order.insert(order.end(), nb.begin(), nb.end());
      }
    }
    std::reverse(order.begin(), order.end());
    return order;
  }


  // ---------------------------------------------------------------------------
  // Jacobian of a system in banded form in the reordered unknowns.
  //   linear systems:     J(t) = columns of f(e_j, t) (grouped probing)
  //   nonlinear systems:  forward finite differences (grouped probing)
  // ---------------------------------------------------------------------------
  template <typename State>
  class sparse_jacobian
  {
  public:
    // dense_limit: systems with n <= dense_limit use a full (dense) pattern
    // which is the only safe choice for nonlinear systems whose pattern
    // depends on the state.
    sparse_jacobian(int n, bool linear, bool time_dependent, int dense_limit = 0) :
      n_(n), linear_(linear), time_dependent_(time_dependent),
      dense_(!linear || n <= dense_limit)
    {
      perm_.resize(n);
      iperm_.resize(n);
      std::iota(perm_.begin(), perm_.end(), 0);
      std::iota(iperm_.begin(), iperm_.end(), 0);
      if (dense_) { kl_ = n - 1; ku_ = n - 1; have_pattern_ = true; }
    }

    int n() const { return n_; }
    int kl() const { return kl_; }
    int ku() const { return ku_; }
    const std::vector<int>& perm() const { return perm_; }

    // evaluates J at (x, t) into the band (J(pi, pj) at band_(ku + pi - pj, pj))
    template <typename Sys>
    void evaluate(Sys& sys, const State& x, double t)
    {
      if (!have_pattern_ || (!linear_ && !dense_)) {
        discover(sys, x, t);
        return;
      }
      // linear time independent systems: the Jacobian is constant
      if (linear_ && !time_dependent_ && evaluated_) return;
      evaluated_ = true;
      grouped_probe(sys, x, t);
      if (linear_ && !dense_ && time_dependent_) {
        // safeguard: the pattern was found at another time
        if (!consistent(sys, t)) discover(sys, x, t);
      }
    }

    // M = scale * I - J   (in LAPACK layout of 'lu'), lu already init()'ed
    void build_matrix(band_lu& lu, double scale) const
    {
      lu.zero();
      for (int pj = 0; pj < n_; ++pj) {
        const int lo = std::max(0, pj - ku_);
        const int hi = std::min(n_ - 1, pj + kl_);
        for (int pi = lo; pi <= hi; ++pi) {
          lu.at(pi, pj) = -band(pi, pj);
        }
        lu.at(pj, pj) += scale;
      }
    }

  private:
    double& band(int pi, int pj)
    {
      return jb_[static_cast<size_t>(ku_ + pi - pj) + static_cast<size_t>(pj) * (kl_ + ku_ + 1)];
    }
    double band(int pi, int pj) const
    {
      return jb_[static_cast<size_t>(ku_ + pi - pj) + static_cast<size_t>(pj) * (kl_ + ku_ + 1)];
    }

    // full column-by-column probing. Establishes the pattern (linear systems
    // with n > dense_limit) and the Jacobian.
    template <typename Sys>
    void discover(Sys& sys, const State& x, double t)
    {
      std::vector<std::pair<int, int>> nz;      // (row, col) in original numbering
      std::vector<double> val;
      State e(n_), dx(n_), f0(n_), xp(n_);
      if (!linear_) sys.eval(x, f0, t);
      for (int j = 0; j < n_; ++j) {
        if (linear_) {
          for (int i = 0; i < n_; ++i) e[i] = 0.0;
          e[j] = 1.0;
          sys.eval(e, dx, t);
          for (int i = 0; i < n_; ++i) {
            if (dx[i] != 0.0) { nz.emplace_back(i, j); val.push_back(dx[i]); }
          }
        } else {
          for (int i = 0; i < n_; ++i) xp[i] = x[i];
          const double h = fd_step(x[j]);
          xp[j] += h;
          sys.eval(xp, dx, t);
          for (int i = 0; i < n_; ++i) {
            nz.emplace_back(i, j); val.push_back((dx[i] - f0[i]) / h);
          }
        }
      }
      if (dense_) {
        kl_ = n_ - 1; ku_ = n_ - 1;
        std::iota(perm_.begin(), perm_.end(), 0);
        std::iota(iperm_.begin(), iperm_.end(), 0);
      } else {
        perm_ = rcm_order(n_, nz);
        for (int p = 0; p < n_; ++p) iperm_[perm_[p]] = p;
        kl_ = 0; ku_ = 0;
        for (const auto& e2 : nz) {
          const int pr = iperm_[e2.first], pc = iperm_[e2.second];
          kl_ = std::max(kl_, pr - pc);
          ku_ = std::max(ku_, pc - pr);
        }
      }
      jb_.assign(static_cast<size_t>(kl_ + ku_ + 1) * n_, 0.0);
      for (size_t k = 0; k < nz.size(); ++k) {
        band(iperm_[nz[k].first], iperm_[nz[k].second]) = val[k];
      }
      have_pattern_ = true;
      evaluated_ = true;
      if (dense_) return;
      // the pattern of linear systems is state independent: remember x-free
      // evaluation is complete, nothing else to do
    }

    static double fd_step(double xj)
    {
      return 1.4901161193847656e-08 * std::max(1.0, std::fabs(xj));
    }

    // grouped probing using the established pattern. One right-hand side
    // evaluation per diagonal of the band.
    template <typename Sys>
    void grouped_probe(Sys& sys, const State& x, double t)
    {
      jb_.assign(static_cast<size_t>(kl_ + ku_ + 1) * n_, 0.0);
      const int m = kl_ + ku_ + 1;
      State v(n_), dx(n_), f0(n_);
      if (!linear_) sys.eval(x, f0, t);
      std::vector<double> h(n_, 0.0);
      for (int g = 0; g < std::min(m, n_); ++g) {
        if (linear_) {
          for (int i = 0; i < n_; ++i) v[i] = 0.0;
          for (int pc = g; pc < n_; pc += m) v[perm_[pc]] = 1.0;
        } else {
          for (int i = 0; i < n_; ++i) v[i] = x[i];
          for (int pc = g; pc < n_; pc += m) {
            h[pc] = fd_step(x[perm_[pc]]);
            v[perm_[pc]] += h[pc];
          }
        }
        sys.eval(v, dx, t);
        for (int pc = g; pc < n_; pc += m) {
          const int lo = std::max(0, pc - ku_);
          const int hi = std::min(n_ - 1, pc + kl_);
          for (int pr = lo; pr <= hi; ++pr) {
            const int r = perm_[pr];
            band(pr, pc) = linear_ ? dx[r] : (dx[r] - f0[r]) / h[pc];
          }
        }
      }
    }

    // J w == f(w) for linear systems (checks that the pattern is complete)
    template <typename Sys>
    bool consistent(Sys& sys, double t) const
    {
      State w(n_), fw(n_);
      double s = 0.37;
      for (int i = 0; i < n_; ++i) { s = std::fmod(s * 7.31 + 0.113, 1.0); w[i] = 0.25 + s; }
      sys.eval(w, fw, t);
      double err = 0.0, nrm = 0.0;
      for (int pr = 0; pr < n_; ++pr) {
        double acc = 0.0;
        const int lo = std::max(0, pr - kl_);
        const int hi = std::min(n_ - 1, pr + ku_);
        for (int pc = lo; pc <= hi; ++pc) acc += band(pr, pc) * w[perm_[pc]];
        const double f = fw[perm_[pr]];
        err = std::max(err, std::fabs(acc - f));
        nrm = std::max(nrm, std::fabs(f));
      }
      return err <= 1e-9 * (nrm + 1e-300) + 1e-300;
    }

    int n_;
    bool linear_, time_dependent_, dense_;
    bool have_pattern_ = false;
    bool evaluated_ = false;
    int kl_ = 0, ku_ = 0;
    std::vector<int> perm_, iperm_;
    std::vector<double> jb_;
  };


  // ---------------------------------------------------------------------------
  // Rosenbrock method, as boost::odeint::rosenbrock4 (coefficients identical)
  // ---------------------------------------------------------------------------
  struct rosenbrock4_coefficients
  {
    static constexpr double gamma = 0.25, d1 = 0.25, d2 = -0.1043, d3 = 0.1035,
      d4 = 0.3620000000000023e-01, c2 = 0.386, c3 = 0.21, c4 = 0.63,
      c21 = -0.5668800000000000e+01, a21 = 0.1544000000000000e+01,
      c31 = -0.2430093356833875e+01, c32 = -0.2063599157091915e+00,
      a31 = 0.9466785280815826e+00, a32 = 0.2557011698983284e+00,
      c41 = -0.1073529058151375e+00, c42 = -0.9594562251023355e+01,
      c43 = -0.2047028614809616e+02, a41 = 0.3314825187068521e+01,
      a42 = 0.2896124015972201e+01, a43 = 0.9986419139977817e+00,
      c51 = 0.7496443313967647e+01, c52 = -0.1024680431464352e+02,
      c53 = -0.3399990352819905e+02, c54 = 0.1170890893206160e+02,
      a51 = 0.1221224509226641e+01, a52 = 0.6019134481288629e+01,
      a53 = 0.1253708332932087e+02, a54 = -0.6878860361058950e+00,
      c61 = 0.8083246795921522e+01, c62 = -0.7981132988064893e+01,
      c63 = -0.3152159432874371e+02, c64 = 0.1631930543123136e+02,
      c65 = -0.6058818238834054e+01;
  };


  // Integrates sys from t0 to t1, replacing y. State: any random access vector
  // with size(), operator[] and a constructor from size (ublas::vector,
  // std::vector<double>).
  //
  // Generic driver. 'prepare(x, t, dt)' must evaluate the Jacobian J at (x, t)
  // and factorize (1/(gamma dt) I - J); 'solve(b)' solves with that matrix in
  // place.
  template <typename Sys, typename State, typename Prepare, typename Solve>
  void integrate_rosenbrock4_impl(Sys& sys, State& y, double t0, double t1,
                                  double atol, double rtol,
                                  Prepare&& prepare, Solve&& solve)
  {
    using C = rosenbrock4_coefficients;
    const int n = static_cast<int>(y.size());
    if (n == 0 || t1 == t0) return;
    State x(n), xnew(n), xerr(n), dxdt(n), dxdtnew(n), dfdt(n), xtmp(n),
      g1(n), g2(n), g3(n), g4(n), g5(n);
    for (int i = 0; i < n; ++i) x[i] = y[i];

    const double dir = (t1 >= t0) ? 1.0 : -1.0;
    double t = t0;
    double dt = 0.1 * (t1 - t0);
    bool first_step = true, last_rejected = false;
    double err_old = 0.0, dt_old = 0.0;
    const double safe = 0.9, fac1 = 5.0, fac2 = 1.0 / 6.0;
    int consecutive_failures = 0;

    while (dir * (t1 - t) > 0.0) {
      // clip step to the end of the interval (as integrate_adaptive)
      bool at_end = false;
      if (dir * (t + dt - t1) >= 0.0) { dt = t1 - t; at_end = true; }

      // ---- one Rosenbrock step from (x, t) with step dt -------------------
      sys(x, dxdt, t);
      if (Sys::stiff_time_dependent) {
        const double h = 1e-6 * std::max(1.0, std::fabs(t));
        sys.eval(x, xtmp, t + h);
        sys.eval(x, xnew, t - h);
        for (int i = 0; i < n; ++i) dfdt[i] = (xtmp[i] - xnew[i]) / (2.0 * h);
      } else {
        for (int i = 0; i < n; ++i) dfdt[i] = 0.0;
      }
      prepare(x, t, dt);

      for (int i = 0; i < n; ++i) g1[i] = dxdt[i] + dt * C::d1 * dfdt[i];
      solve(g1);

      for (int i = 0; i < n; ++i) xtmp[i] = x[i] + C::a21 * g1[i];
      sys(xtmp, dxdtnew, t + C::c2 * dt);
      for (int i = 0; i < n; ++i) g2[i] = dxdtnew[i] + dt * C::d2 * dfdt[i] + C::c21 * g1[i] / dt;
      solve(g2);

      for (int i = 0; i < n; ++i) xtmp[i] = x[i] + C::a31 * g1[i] + C::a32 * g2[i];
      sys(xtmp, dxdtnew, t + C::c3 * dt);
      for (int i = 0; i < n; ++i) g3[i] = dxdtnew[i] + dt * C::d3 * dfdt[i] + (C::c31 * g1[i] + C::c32 * g2[i]) / dt;
      solve(g3);

      for (int i = 0; i < n; ++i) xtmp[i] = x[i] + C::a41 * g1[i] + C::a42 * g2[i] + C::a43 * g3[i];
      sys(xtmp, dxdtnew, t + C::c4 * dt);
      for (int i = 0; i < n; ++i) g4[i] = dxdtnew[i] + dt * C::d4 * dfdt[i] + (C::c41 * g1[i] + C::c42 * g2[i] + C::c43 * g3[i]) / dt;
      solve(g4);

      for (int i = 0; i < n; ++i) xtmp[i] = x[i] + C::a51 * g1[i] + C::a52 * g2[i] + C::a53 * g3[i] + C::a54 * g4[i];
      sys(xtmp, dxdtnew, t + dt);
      for (int i = 0; i < n; ++i) g5[i] = dxdtnew[i] + (C::c51 * g1[i] + C::c52 * g2[i] + C::c53 * g3[i] + C::c54 * g4[i]) / dt;
      solve(g5);

      for (int i = 0; i < n; ++i) xtmp[i] += g5[i];
      sys(xtmp, dxdtnew, t + dt);
      for (int i = 0; i < n; ++i) xerr[i] = dxdtnew[i] + (C::c61 * g1[i] + C::c62 * g2[i] + C::c63 * g3[i] + C::c64 * g4[i] + C::c65 * g5[i]) / dt;
      solve(xerr);
      for (int i = 0; i < n; ++i) xnew[i] = xtmp[i] + xerr[i];

      // ---- step size control (rosenbrock4_controller) ----------------------
      double err = 0.0;
      for (int i = 0; i < n; ++i) {
        const double sk = atol + rtol * std::max(std::fabs(x[i]), std::fabs(xnew[i]));
        err += xerr[i] * xerr[i] / sk / sk;
      }
      err = std::sqrt(err / static_cast<double>(n));
      if (!std::isfinite(err)) err = 1e10;

      double fac = std::max(fac2, std::min(fac1, std::pow(err, 0.25) / safe));
      double dt_new = dt / fac;
      if (err <= 1.0) {
        if (first_step) {
          first_step = false;
        } else {
          double fac_pred = (dt_old / dt) * std::pow(err * err / err_old, 0.25) / safe;
          fac_pred = std::max(fac2, std::min(fac1, fac_pred));
          fac = std::max(fac, fac_pred);
          dt_new = dt / fac;
        }
        dt_old = dt;
        err_old = std::max(0.01, err);
        if (last_rejected) dt_new = (dt >= 0.0) ? std::min(dt_new, dt) : std::max(dt_new, dt);
        t = at_end ? t1 : t + dt;
        for (int i = 0; i < n; ++i) x[i] = xnew[i];
        dt = dt_new;
        last_rejected = false;
        consecutive_failures = 0;
      } else {
        dt = dt_new;
        last_rejected = true;
        if (++consecutive_failures > 100 || std::fabs(dt) < 1e-300) {
          throw std::runtime_error("daisie_stiff: step size underflow in rosenbrock4");
        }
      }
    }
    for (int i = 0; i < n; ++i) y[i] = x[i];
  }


  // Integrates sys from t0 to t1, replacing y, with the banded Jacobian
  // machinery above (pattern detection by probing).
  template <typename Sys, typename State>
  void integrate_rosenbrock4(Sys& sys, State& y, double t0, double t1,
                             double atol, double rtol, int dense_limit = 0,
                             int max_sparse_n = 6000)
  {
    const int n = static_cast<int>(y.size());
    if (n == 0 || t1 == t0) return;
    {
      // The pattern of the Jacobian is detected by probing with all n unit
      // vectors, which is quadratic in n; a dense Jacobian is cubic in n.
      const bool dense = !Sys::stiff_linear || n <= dense_limit;
      if ((dense && n > 500) || (!dense && n > max_sparse_n)) {
        throw std::runtime_error(
          "odeint::rosenbrock4 is not available for systems of this size "
          "(n = " + std::to_string(n) + "); use another method, e.g. "
          "odeint::runge_kutta_cash_karp54 or lsodes.");
      }
    }
    sparse_jacobian<State> jac(n, Sys::stiff_linear, Sys::stiff_time_dependent, dense_limit);
    band_lu lu;
    State work(n);
    auto prepare = [&](const State& x, double t, double dt) {
      jac.evaluate(sys, x, t);
      lu.init(n, jac.kl(), jac.ku());
      jac.build_matrix(lu, 1.0 / rosenbrock4_coefficients::gamma / dt);
      lu.factorize();
    };
    auto solve = [&](State& b) {
      const auto& perm = jac.perm();
      for (int p = 0; p < n; ++p) work[p] = b[perm[p]];
      lu.solve(work);
      for (int p = 0; p < n; ++p) b[perm[p]] = work[p];
    };
    integrate_rosenbrock4_impl(sys, y, t0, t1, atol, rtol, prepare, solve);
  }

}  // namespace daisie_stiff

#endif
