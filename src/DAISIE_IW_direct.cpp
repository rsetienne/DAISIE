//
//  DAISIE_IW_direct.cpp
//
//  The island-wide (IW) likelihood with the right-hand side of the ODE and the event step evaluated "on the fly"
//  (option rhs = "direct" of DAISIE_loglik_IW() and DAISIE_ML_IW()).
//
//  WHAT THE TABLES VERSION DOES (rhs = "tables", src/DAISIE_IW.cpp and R/DAISIE_loglik_IW.R).  For every ODE segment R builds
//  eight coefficient arrays c1..c8 of the size of the state space (lxm x lxe x sysdim doubles each), six helper arrays for the
//  event step, and the matrix ki of size sysdim x sysdim (also copied in C++), and the event step between two segments works
//  on further arrays of the size of the state space.
//
//  WHAT THIS FILE DOES.  Nothing of the size of the state space is stored except the state vector itself (and the vectors that
//  the odeint stepper needs).  A state is (i, j, a): i = n1 in 0..lxm-1 and j = n2 in 0..lxe-1 (the numbers of the two kinds of
//  species that are not tracked), a in 0..2^s-1 the bit pattern of the s open clades.  It is stored at
//        index = i + lxm * (j + lxe * a)          (0-based; 1 + index in R)
//  which is the layout of the R array with dim c(lxm, lxe, sysdim) (column-major).  All coefficients depend on (i, j) and on the
//  number l0 of zero bits of a only, so they are computed from the state with a few operations, with the same formulas, in the
//  same order of the operations, as DAISIE_IW_pars() and nndivdep() in R.  The only stored tables are the three diversity-
//  dependence factors as a function of i + j (length lxm + lxe - 1).  The matrix ki, whose row a has a one for every state
//  that is a with one set bit cleared, is replaced by a loop over the set bits of a.
//
//  The right-hand side is in "gather" form: every state computes its derivative from its (at most 7 + s) source states,
//        dx(i,j,a) = c1 x(i-1,j,a) + c2 x(i+1,j,a) + c3 x(i,j+1,a) + c4 x(i+1,j-1,a) + c5 x(i+1,j-2,a) + c6 x(i,j-1,a)
//                    - c7 x(i,j,a) + sum over the set bits b of a of ( laa x(i,j,a-2^b) + c8 x(i,j-1,a-2^b) ),
//  with x = 0 outside the box (the mass that leaves the box is dropped, as in the tables version, where x is padded with zeros).
//  This is exactly the stencil of cpp_daisie_iw<>::rhs in DAISIE_IW.cpp; the rows (a, j) are independent, which makes the
//  evaluation parallel without synchronisation.  (A "scatter" version, a loop over the source states, would need atomic
//  additions to be parallel; it gives the same result.)
//
//  All sizes and indices are 64 bit (std::int64_t).  R checks that the number of states is below 2^31 before it calls this
//  code (R/DAISIE_loglik_IW_direct.R).
//
//  COMPACT STORAGE (rhs = "compact").  With a finite K' every transition that increases n = i + j has the rate factor
//  max(0, 1 - (n + k) / K'), which is exactly 0 when n + k >= K'.  Starting from (0, 0) at k = 0, and because the events (which
//  multiply by the same factor) increase k by 1, the states with probability > 0 in the segment with k events so far all have
//        i + j <= ceiling(K') - k =: tri.
//  The states with i + j > tri are exactly 0 in the tables version (and in the version above), so they need not be stored.  The
//  states of a segment are then stored per configuration a, per j, with the i = 0..len[j]-1 that satisfy i + j <= tri and
//  i < lxm (len[j] = min(lxm, tri - j + 1); the rows with len[j] <= 0 are absent):
//        index = a * P + off[j] + i,     off[j] = len[0] + ... + len[j-1],     P = sum of len[j].
//  With tri >= lxm + lxe - 2 (in particular tri = Inf) this is the box layout above, so "direct" is the special case tri = Inf.
//  The arithmetic is that of the box version (the states that are left out would only add exact zeros), in the same order, so
//  the result is the same to the last bit.  The event step reads the layout of segment k (tri) and writes that of segment k + 1
//  (tri - 1; the same if tri = Inf).  The sparse Jacobian of odeint::rosenbrock4 assumes the box layout, so it can only be used
//  when the layout is the box.
//
#include "config.h"
#include "DAISIE_odeint.h"
#define EIGEN_USE_THREADS
#include <unsupported/Eigen/CXX11/ThreadPool>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// defined in DAISIE_IW.cpp (number of threads of the IW code, DAISIE_IW_num_threads())
extern unsigned daisie_odeint_iw_num_threads_;

// defined in DAISIE_IW_sparse.cpp
void daisie_iw_sparse_rosenbrock4(const std::function<void(const double*, double*)>& rhs,
  long d0, long d1, long S, std::vector<double>& y, double t0, double t1, double atol, double rtol);

namespace {

  using i64 = std::int64_t;

  inline int popcount64(std::uint64_t v)
  {
    int c = 0;
    while (v) { v &= v - 1; ++c; }
    return c;
  }

  // pmax(0, v) of R: NaN stays NaN
  inline double pmax0(double v)
  {
    return (v != v) ? v : (v > 0.0 ? v : 0.0);
  }

  // "no bound": a triangle bound that is larger than any box
  const i64 tri_inf = std::numeric_limits<i64>::max() / 4;

  // The states of one configuration a: the rows j = 0..nrow-1, each with i = 0..len[j]-1 (see COMPACT STORAGE above)
  struct iw_layout
  {
    i64 lxm = 0, lxe = 0, tri = 0, nrow = 0, P = 0;
    std::vector<i64> len, off;

    iw_layout() {}
    iw_layout(i64 lxm_, i64 lxe_, i64 tri_) : lxm(lxm_), lxe(lxe_), tri(tri_)
    {
      for (i64 j = 0; j < lxe; ++j) {
        const i64 l = std::min(lxm, tri - j + 1);
        if (l <= 0) break;                       // len[j] does not increase with j
        len.push_back(l);
      }
      nrow = static_cast<i64>(len.size());
      off.assign(static_cast<size_t>(nrow) + 1, 0);
      for (i64 j = 0; j < nrow; ++j) off[static_cast<size_t>(j) + 1] = off[static_cast<size_t>(j)] + len[static_cast<size_t>(j)];
      P = off[static_cast<size_t>(nrow)];
    }
    bool is_box() const { return P == lxm * lxe; }
    // 0-based index of (i, j) within a configuration, -1 if the state is not stored
    i64 pos(i64 i, i64 j) const
    {
      if (i < 0 || j < 0 || j >= nrow || i >= len[static_cast<size_t>(j)]) return -1;
      return off[static_cast<size_t>(j)] + i;
    }
  };

  inline i64 tri_from_double(double d)
  {
    if (!(d == d)) throw std::runtime_error("daisie_iw_direct: tri is NaN");
    if (d >= static_cast<double>(tri_inf)) return tri_inf;
    if (d < -1.0) return -1;
    return static_cast<i64>(std::floor(d));
  }

  // The parameters of one ODE segment (k events so far, s open clades) and the only tables:
  // the diversity-dependence factors as functions of nn = i + j.
  struct iw_par
  {
    double lac, mu, Kprime, gam, laa, M, k;
    bool recol_dd;             // TRUE: recolonisation is diversity dependent (DAISIE <= 5.0.3), see DAISIE_recol_dd()
    i64 lxm, lxe, S, P, n;     // box, number of configurations 2^s, states per configuration, number of states
    int s;
    iw_layout lay, lay_next;   // the layout of this segment, and of the next (after an event)
    std::vector<double> dd, ddm, ddp;       // divdepfac, divdepfacmin1, divdepfacplus1 at nn = i + j

    explicit iw_par(Rcpp::List p)
    {
      lac = Rcpp::as<double>(p["lac"]);
      mu = Rcpp::as<double>(p["mu"]);
      Kprime = Rcpp::as<double>(p["Kprime"]);
      gam = Rcpp::as<double>(p["gam"]);
      laa = Rcpp::as<double>(p["laa"]);
      M = Rcpp::as<double>(p["M"]);
      k = Rcpp::as<double>(p["k"]);
      recol_dd = p.containsElementNamed("recol_dd") ? Rcpp::as<bool>(p["recol_dd"]) : false;
      const double dlxm = Rcpp::as<double>(p["lxm"]), dlxe = Rcpp::as<double>(p["lxe"]), ds = Rcpp::as<double>(p["s"]);
      if (!(dlxm >= 1.0 && dlxe >= 1.0 && ds >= 0.0 && ds <= 40.0) || dlxm != std::floor(dlxm) || dlxe != std::floor(dlxe) || ds != std::floor(ds)) {
        throw std::runtime_error("daisie_iw_direct: invalid dimensions lxm, lxe or s");
      }
      lxm = static_cast<i64>(dlxm); lxe = static_cast<i64>(dlxe); s = static_cast<int>(ds);
      S = i64(1) << s;
      const i64 imax = std::numeric_limits<i64>::max();
      if (lxm > imax / lxe) throw std::runtime_error("daisie_iw_direct: lxm * lxe overflows 64 bit");
      // tri: the states with i + j > tri are not stored (default: all, the box)
      const i64 tri = p.containsElementNamed("tri") ? tri_from_double(Rcpp::as<double>(p["tri"])) : tri_inf;
      lay = iw_layout(lxm, lxe, tri);
      lay_next = iw_layout(lxm, lxe, (tri >= tri_inf) ? tri : tri - 1);
      if (lay.P < 1) throw std::runtime_error("daisie_iw_direct: no state is left (tri < 0)");
      P = lay.P;
      if (P > imax / S) throw std::runtime_error("daisie_iw_direct: the number of states overflows 64 bit");
      n = P * S;
      // as nndivdep(): divdepfac = pmax(0, 1 - (nn + k) / Kprime), with k - 1 and k + 1 for the two others
      const i64 nnn = lxm + lxe - 1;
      dd.resize(static_cast<size_t>(nnn)); ddm.resize(static_cast<size_t>(nnn)); ddp.resize(static_cast<size_t>(nnn));
      for (i64 q = 0; q < nnn; ++q) {
        const double nn = static_cast<double>(q);
        dd[static_cast<size_t>(q)] = pmax0(1.0 - (nn + k) / Kprime);
        ddm[static_cast<size_t>(q)] = pmax0(1.0 - (nn + k - 1.0) / Kprime);
        ddp[static_cast<size_t>(q)] = pmax0(1.0 - (nn + k + 1.0) / Kprime);
      }
    }

    // index of the state (i, j, a), 0-based, in 64 bit; -1 if it is not stored
    i64 index(i64 i, i64 j, i64 a) const
    {
      const i64 q = lay.pos(i, j);
      return (q < 0 || a < 0 || a >= S) ? -1 : a * P + q;
    }
  };


  // Splits the range 0..n over the threads of the pool; the calling thread does the first chunk.
  template <typename F>
  void parallel_chunks(Eigen::ThreadPool* pool, int nthreads, i64 n, F&& f)
  {
    if (!pool || nthreads <= 1 || n < 2 * static_cast<i64>(nthreads)) {
      f(i64(0), n);
      return;
    }
    Eigen::Barrier bar(static_cast<unsigned int>(nthreads - 1));
    for (int c = 1; c < nthreads; ++c) {
      const i64 lo = n * c / nthreads, hi = n * (c + 1) / nthreads;
      pool->Schedule([&f, lo, hi, &bar]() { f(lo, hi); bar.Notify(); });
    }
    f(i64(0), n / nthreads);
    bar.Wait();
  }


  // The right-hand side without tables.
  class cpp_daisie_iw_direct
  {
  public:
    explicit cpp_daisie_iw_direct(Rcpp::List pars) : p_(pars) {}

    // extents of the state tensor, as cpp_daisie_iw<>::extents() (only meaningful for the box layout)
    std::array<i64, 3> extents() const { return { p_.lxm, p_.lxe, p_.S }; }
    i64 size() const { return p_.n; }
    bool is_box() const { return p_.lay.is_box(); }

    void rhs(const double* x, double* dx, Eigen::ThreadPool* pool, int nthreads) const
    {
      parallel_chunks(pool, nthreads, p_.lay.nrow * p_.S,
                      [&](i64 r0, i64 r1) { rows(x, dx, r0, r1); });
    }

  private:
    // the rows r = a * nrow + j (all stored states (., j, a)) from r0 to r1
    void rows(const double* x, double* dx, i64 r0, i64 r1) const
    {
      const iw_layout& L = p_.lay;
      const i64 nrow = L.nrow, P = p_.P;
      const double lac = p_.lac, mu = p_.mu, gam = p_.gam, laa = p_.laa, M = p_.M, k = p_.k;
      const double* dd = p_.dd.data();
      const double* ddm = p_.ddm.data();
      const double* nbj[64];
      const double* nbjm1[64];
      for (i64 r = r0; r < r1; ++r) {
        const i64 a = r / nrow, j = r - a * nrow;
        const int l0 = p_.s - popcount64(static_cast<std::uint64_t>(a));       // number of clades with bit 0
        const i64 len = L.len[static_cast<size_t>(j)];                          // the states (i, j) with i < len are stored
        const i64 lm1 = (j >= 1) ? L.len[static_cast<size_t>(j - 1)] : 0;
        const i64 lm2 = (j >= 2) ? L.len[static_cast<size_t>(j - 2)] : 0;
        const i64 lp1 = (j + 1 < nrow) ? L.len[static_cast<size_t>(j + 1)] : 0;
        const double* base = x + a * P;
        const double* xr = base + L.off[static_cast<size_t>(j)];                                              // x(., j, a)
        const double* xrm1 = (j >= 1) ? base + L.off[static_cast<size_t>(j - 1)] : nullptr;                  // x(., j - 1, a)
        const double* xrm2 = (j >= 2) ? base + L.off[static_cast<size_t>(j - 2)] : nullptr;                  // x(., j - 2, a)
        const double* xrp1 = (j + 1 < nrow) ? base + L.off[static_cast<size_t>(j + 1)] : nullptr;            // x(., j + 1, a)
        double* out = dx + a * P + L.off[static_cast<size_t>(j)];
        // the states a with one set bit cleared, in order of increasing index (as the columns of the matrix ki)
        int nnb = 0;
        for (int b = p_.s - 1; b >= 0; --b) {
          if ((a >> b) & 1) {
            const double* q = x + (a ^ (i64(1) << b)) * P;
            nbj[nnb] = q + L.off[static_cast<size_t>(j)];
            nbjm1[nnb] = (j >= 1) ? q + L.off[static_cast<size_t>(j - 1)] : nullptr;
            ++nnb;
          }
        }
        const double dj = static_cast<double>(j);
        const double jprev = (j >= 1) ? dj - 1.0 : 0.0;                         // nn[nilm, nil2lxe - 1, allc]
        for (i64 i = 0; i < len; ++i) {
          const double di = static_cast<double>(i);
          const i64 nn = i + j;
          const double d0 = dd[nn], dm = ddm[nn];
          const double iprev = (i >= 1) ? di - 1.0 : 0.0;                       // nn[nil2lxm - 1, nile, allc]
          // the coefficients c1..c8 of DAISIE_IW_pars(), with the same order of the operations
          const double c1 = gam * dm * pmax0(M - l0 - iprev);
          const double c2 = mu * (di + 1.0);
          const double c3 = mu * (dj + 1.0);
          const double c4 = laa * (di + 1.0);
          const double c5 = lac * dm * (di + 1.0);
          const double c6 = lac * dm * (2.0 * (k - l0) + jprev);
          // colonisation out of the state: by the absent mainland species (diversity dependent) and recolonisation of the l0
          // observed non-endemics (not diversity dependent unless recol_dd), as immig_out in DAISIE_IW_pars()
          const double immig_out = p_.recol_dd ? gam * d0 * pmax0(M - di) : gam * d0 * pmax0(M - l0 - di) + gam * l0;
          const double c7 = immig_out + mu * (k + (di + dj)) + laa * (l0 + di) + lac * d0 * (k + (di + dj));
          const double c8 = 2.0 * lac * dm;
          // the sources of the terms; zero where they are not stored (outside the box, or outside the triangle where they are 0)
          const double x_im1 = (i >= 1) ? xr[i - 1] : 0.0;
          const double x_ip1 = (i + 1 < len) ? xr[i + 1] : 0.0;
          const double x_jp1 = (xrp1 && i < lp1) ? xrp1[i] : 0.0;
          const double x_ip1_jm1 = (xrm1 && i + 1 < lm1) ? xrm1[i + 1] : 0.0;
          const double x_ip1_jm2 = (xrm2 && i + 1 < lm2) ? xrm2[i + 1] : 0.0;
          const double x_jm1 = (xrm1 && i < lm1) ? xrm1[i] : 0.0;
          double v = c1 * x_im1 + c2 * x_ip1 + c3 * x_jp1 + c4 * x_ip1_jm1 + c5 * x_ip1_jm2 + c6 * x_jm1 - c7 * xr[i];
          if (nnb > 0) {
            double acc = 0.0;
            for (int t = 0; t < nnb; ++t) {
              acc += laa * nbj[t][i] + c8 * ((nbjm1[t] && i < lm1) ? nbjm1[t][i] : 0.0);
            }
            v += acc;
          }
          out[i] = v;
        }
      }
    }

    iw_par p_;
  };


  struct iw_direct_wrapper
  {
    using jacobian = daisie_odeint::jacobian_policy::const_from_linear_rhs<iw_direct_wrapper>;

    std::unique_ptr<Eigen::ThreadPool> pool;
    int nthreads;
    cpp_daisie_iw_direct sys;

    explicit iw_direct_wrapper(Rcpp::List pars) : nthreads(1), sys(pars)
    {
      nthreads = static_cast<int>(std::max(1u, daisie_odeint_iw_num_threads_));
      if (nthreads > 1) pool.reset(new Eigen::ThreadPool(nthreads));
    }

    static constexpr bool stiff_linear = true;
    static constexpr bool stiff_time_dependent = false;

    int size() const { return static_cast<int>(sys.size()); }

    // Rosenbrock integration with the sparse LU of DAISIE_IW_sparse.cpp, as for the tables version
    void integrate_rosenbrock4(state_type& y, double t0, double t1, double atol, double rtol)
    {
      if (!sys.is_box()) {
        throw std::runtime_error("odeint::rosenbrock4 needs the box layout of the states, not the compact storage of rhs = \"compact\"");
      }
      const auto e = sys.extents();
      std::vector<double> v(y.begin(), y.end());
      auto f = [this](const double* x, double* dx) { sys.rhs(x, dx, pool.get(), nthreads); };
      daisie_iw_sparse_rosenbrock4(f, static_cast<long>(e[0]), static_cast<long>(e[1]), static_cast<long>(e[2]), v, t0, t1, atol, rtol);
      for (std::size_t i = 0; i < v.size(); ++i) y[i] = v[i];
    }

    void eval(const state_type& x, state_type& dxdt, double t) { (*this)(x, dxdt, t); }

    // odeint interface
    void operator()(const state_type& x, state_type& dxdt, double)
    {
      sys.rhs(x.data().begin(), dxdt.data().begin(), pool.get(), nthreads);
    }
  };

} // anonymous namespace


// Driver for the boost::odeint solver for the IW model, right-hand side without tables
//
// As \code{daisie_odeint_iw()}, but \code{rpars} is a short list with the scalars lac, mu, Kprime, gam, laa, M, k (the
// number of colonisations so far), lxm, lxe and s (the number of open clades; the state space has lxm * lxe * 2^s states).

RcppExport SEXP daisie_odeint_iw_direct(SEXP ry, SEXP rtimes, SEXP rpars, SEXP Stepper, SEXP atolint, SEXP reltolint)
{
BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  auto pars = Rcpp::as<Rcpp::List>(rpars);
  iw_direct_wrapper iw(pars);
  const Rcpp::NumericVector ry_(ry);
  if (static_cast<i64>(ry_.size()) != iw.sys.size()) {
    throw std::runtime_error("daisie_odeint_iw_direct: the length of the state vector is not lxm * lxe * 2^s");
  }
  auto y = Rcpp::as<state_type>(ry);
  auto times = Rcpp::as<std::vector<double>>(rtimes);
  auto stepper = Rcpp::as<std::string>(Stepper);
  auto atol = Rcpp::as<double>(atolint);
  auto rtol = Rcpp::as<double>(reltolint);
  if ("odeint::rosenbrock4" == stepper) {
    iw.integrate_rosenbrock4(y, times[0], times[1], atol, rtol);
  } else {
    daisie_odeint::integrate(stepper, std::ref(iw), y, times[0], times[1], atol, rtol);
  }
  rcpp_result_gen = Rcpp::wrap(y);
  return rcpp_result_gen;
END_RCPP
}


// One evaluation of the right-hand side without tables (used to test it against the tables)

RcppExport SEXP daisie_iw_rhs_eval_direct(SEXP ry, SEXP rpars)
{
BEGIN_RCPP
  auto pars = Rcpp::as<Rcpp::List>(rpars);
  iw_direct_wrapper iw(pars);
  const Rcpp::NumericVector x(ry);
  if (static_cast<i64>(x.size()) != iw.sys.size()) {
    throw std::runtime_error("daisie_iw_rhs_eval_direct: wrong length of x");
  }
  Rcpp::NumericVector dx(x.size());
  iw.sys.rhs(x.begin(), dx.begin(), iw.pool.get(), iw.nthreads);
  return dx;
END_RCPP
}


// The index of states, 1-based, computed in 64 bit (used to test that R and C++ number the states in the same way)

RcppExport SEXP daisie_iw_state_index(SEXP ri, SEXP rj, SEXP ra, SEXP rlxm, SEXP rlxe)
{
BEGIN_RCPP
  const Rcpp::NumericVector vi(ri), vj(rj), va(ra);
  const i64 lxm = static_cast<i64>(Rcpp::as<double>(rlxm)), lxe = static_cast<i64>(Rcpp::as<double>(rlxe));
  Rcpp::NumericVector out(vi.size());
  for (R_xlen_t t = 0; t < vi.size(); ++t) {
    const i64 i = static_cast<i64>(vi[t]), j = static_cast<i64>(vj[t]), a = static_cast<i64>(va[t]);
    out[t] = static_cast<double>(1 + i + lxm * (j + lxe * a));
  }
  return out;
END_RCPP
}


// The layout of the states of one configuration for a bound tri (i + j <= tri; Inf for the box): the lengths len[j] of the
// rows, their offsets and the number P of states per configuration (used to test that R and C++ number the states in the same way)

RcppExport SEXP daisie_iw_layout(SEXP rlxm, SEXP rlxe, SEXP rtri)
{
BEGIN_RCPP
  const i64 lxm = static_cast<i64>(Rcpp::as<double>(rlxm)), lxe = static_cast<i64>(Rcpp::as<double>(rlxe));
  const iw_layout L(lxm, lxe, tri_from_double(Rcpp::as<double>(rtri)));
  Rcpp::NumericVector len(L.nrow), off(L.nrow);
  for (i64 j = 0; j < L.nrow; ++j) { len[j] = static_cast<double>(L.len[static_cast<size_t>(j)]); off[j] = static_cast<double>(L.off[static_cast<size_t>(j)]); }
  return Rcpp::List::create(Rcpp::Named("len") = len, Rcpp::Named("off") = off, Rcpp::Named("P") = static_cast<double>(L.P));
END_RCPP
}


// The 1-based index of the states (i, j, a) in the layout of tri, in 64 bit; NA if the state is not stored

RcppExport SEXP daisie_iw_state_index_tri(SEXP ri, SEXP rj, SEXP ra, SEXP rlxm, SEXP rlxe, SEXP rtri)
{
BEGIN_RCPP
  const Rcpp::NumericVector vi(ri), vj(rj), va(ra);
  const i64 lxm = static_cast<i64>(Rcpp::as<double>(rlxm)), lxe = static_cast<i64>(Rcpp::as<double>(rlxe));
  const iw_layout L(lxm, lxe, tri_from_double(Rcpp::as<double>(rtri)));
  Rcpp::NumericVector out(vi.size());
  for (R_xlen_t t = 0; t < vi.size(); ++t) {
    const i64 q = L.pos(static_cast<i64>(vi[t]), static_cast<i64>(vj[t]));
    out[t] = (q < 0) ? NA_REAL : static_cast<double>(1 + static_cast<i64>(va[t]) * L.P + q);
  }
  return out;
END_RCPP
}


// The event between two ODE segments of the IW likelihood, without arrays of the size of the state space in R
//
// \code{rtype} 1: colonisation; the number of configurations doubles (a new open clade is the new highest bit; the new
// states have bit 0).  2: speciation (multiplication by lac * divdepfac).  3: speciation that is the first of a clade; as 2,
// and then the bit \code{rbit} (0-based) is removed by adding the states in which it is 0 and 1.  The scalars k, s and tri of
// \code{rpars} are those of the segment that has just been solved: the state vector \code{ry} has the layout of tri, and the
// result has the layout of the next segment (tri - 1, or the same box if tri is infinite).

RcppExport SEXP daisie_iw_event(SEXP ry, SEXP rtype, SEXP rpars, SEXP rbit)
{
BEGIN_RCPP
  auto pars = Rcpp::as<Rcpp::List>(rpars);
  const iw_par p(pars);
  const int type = Rcpp::as<int>(rtype);
  const Rcpp::NumericVector y(ry);
  if (static_cast<i64>(y.size()) != p.n) throw std::runtime_error("daisie_iw_event: wrong length of the state vector");
  const iw_layout& Ls = p.lay;           // layout of the segment that has been solved
  const iw_layout& Ld = p.lay_next;      // layout of the next segment
  if (Ld.P < 1) throw std::runtime_error("daisie_iw_event: no state is left in the next segment (tri < 0)");
  const i64 Ps = Ls.P, Pd = Ld.P;
  const double lac = p.lac, gam = p.gam, M = p.M;
  const i64 imax = std::numeric_limits<i64>::max();
  if (Pd > imax / (2 * p.S)) throw std::runtime_error("daisie_iw_event: the number of states overflows 64 bit");
  if (type == 1) {
    Rcpp::NumericVector out(static_cast<R_xlen_t>(2 * p.S * Pd));  // zero: the new configurations have bit 1 (the first half has bit 0)
    for (i64 a = 0; a < p.S; ++a) {
      const int l0 = p.s - popcount64(static_cast<std::uint64_t>(a));
      const double denom = M - l0;
      for (i64 j = 0; j < Ld.nrow; ++j) {
        const i64 src = a * Ps + Ls.off[static_cast<size_t>(j)], dst = a * Pd + Ld.off[static_cast<size_t>(j)];
        const i64 lens = Ls.len[static_cast<size_t>(j)], lend = Ld.len[static_cast<size_t>(j)];
        for (i64 i = 0; i < lend; ++i) {
          const double di = static_cast<double>(i);
          const i64 nn = i + j;
          const double mfac = (di + 1.0) / denom;
          const double omf = (M - di - l0) / denom;
          const double p1 = y[static_cast<R_xlen_t>(src + i)];
          const double p2 = (i + 1 < lens) ? y[static_cast<R_xlen_t>(src + i + 1)] : 0.0;        // probs2: shifted by one in n1
          // the second term: recolonisation of an unobserved non-endemic (not diversity dependent unless recol_dd)
          const double recolfac = p.recol_dd ? p.ddp[static_cast<size_t>(nn)] : 1.0;
          out[static_cast<R_xlen_t>(dst + i)] = gam * p.dd[static_cast<size_t>(nn)] * omf * p1 + gam * recolfac * mfac * p2;
        }
      }
    }
    return out;
  }
  if (type == 2) {
    Rcpp::NumericVector out(static_cast<R_xlen_t>(p.S * Pd));
    for (i64 a = 0; a < p.S; ++a) {
      for (i64 j = 0; j < Ld.nrow; ++j) {
        const i64 src = a * Ps + Ls.off[static_cast<size_t>(j)], dst = a * Pd + Ld.off[static_cast<size_t>(j)];
        for (i64 i = 0; i < Ld.len[static_cast<size_t>(j)]; ++i) {
          out[static_cast<R_xlen_t>(dst + i)] = lac * p.dd[static_cast<size_t>(i + j)] * y[static_cast<R_xlen_t>(src + i)];
        }
      }
    }
    return out;
  }
  if (type == 3) {
    const int bit = Rcpp::as<int>(rbit);
    if (p.s < 1 || bit < 0 || bit >= p.s) throw std::runtime_error("daisie_iw_event: invalid bit to remove");
    const i64 half = p.S / 2, low = (i64(1) << bit) - 1;
    Rcpp::NumericVector out(static_cast<R_xlen_t>(half * Pd));
    for (i64 a2 = 0; a2 < half; ++a2) {
      const i64 a0 = ((a2 >> bit) << (bit + 1)) | (a2 & low);
      const i64 a1 = a0 | (i64(1) << bit);
      for (i64 j = 0; j < Ld.nrow; ++j) {
        const i64 b0 = a0 * Ps + Ls.off[static_cast<size_t>(j)], b1 = a1 * Ps + Ls.off[static_cast<size_t>(j)];
        const i64 bo = a2 * Pd + Ld.off[static_cast<size_t>(j)];
        for (i64 i = 0; i < Ld.len[static_cast<size_t>(j)]; ++i) {
          const double f = lac * p.dd[static_cast<size_t>(i + j)];
          out[static_cast<R_xlen_t>(bo + i)] = f * y[static_cast<R_xlen_t>(b0 + i)] + f * y[static_cast<R_xlen_t>(b1 + i)];
        }
      }
    }
    return out;
  }
  throw std::runtime_error("daisie_iw_event: unknown type of event");
END_RCPP
}


// The check and normalisation of the probabilities after a step (as checkprobs2() in R, in C++)
//
// Negative values are set to 0.  Returns a list with the (normalised) vector, the sum and a status: 0 if all is well, 1 if the
// vector contains NaN or NA, 2 if the sum is not positive (the vector is then returned without normalisation).

RcppExport SEXP daisie_iw_normalize(SEXP ry)
{
BEGIN_RCPP
  const Rcpp::NumericVector y(ry);
  const R_xlen_t n = y.size();
  long double sum = 0.0L;
  bool nan = false;
  for (R_xlen_t t = 0; t < n; ++t) {
    const double v = y[t];
    if (v != v) nan = true; else if (v > 0.0) sum += v;
  }
  int status = 0;
  double sp = static_cast<double>(sum);
  if (nan) status = 1; else if (!(sp > 0.0)) status = 2;
  Rcpp::NumericVector out(n);
  if (status == 0) {
    for (R_xlen_t t = 0; t < n; ++t) { const double v = y[t]; out[t] = (v > 0.0 ? v : 0.0) / sp; }
  } else {
    for (R_xlen_t t = 0; t < n; ++t) { const double v = y[t]; out[t] = (v != v) ? v : (v > 0.0 ? v : 0.0); }
  }
  return Rcpp::List::create(Rcpp::Named("probs") = out, Rcpp::Named("sum") = sp, Rcpp::Named("status") = status);
END_RCPP
}
