// Sparse Rosenbrock4 integration for the IW model (Eigen SparseLU).
// Kept apart from DAISIE_IW.cpp (Eigen Tensor + odeint) so that no single
// object file becomes too large for MinGW.
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>
#include "DAISIE_stiff.h"

namespace {

  using vec = std::vector<double>;

  struct iw_system
  {
    static constexpr bool stiff_linear = true;
    static constexpr bool stiff_time_dependent = false;
    const std::function<void(const double*, double*)>& f;
    long d0, d1, S;

    int size() const { return static_cast<int>(d0 * d1 * S); }
    void operator()(const vec& x, vec& dx, double) const { f(x.data(), dx.data()); }
    void eval(const vec& x, vec& dx, double t) const { (*this)(x, dx, t); }

    // Jacobian by structured probing. The right-hand side couples cell
    // (i, j, .) with cells (i + di, j + dj, .), di in {-1,0,1}, dj in
    // {-2,..,1}, and all slices of the same cell. Columns whose cell indices
    // differ by multiples of (3, 4) never write into the same row.
    Eigen::SparseMatrix<double> jacobian() const
    {
      const int n = size();
      auto idx = [&](long i, long j, long a) { return static_cast<int>(i + d0 * (j + d1 * a)); };
      constexpr long Pi = 3, Pj = 4;
      constexpr long shift_i[3] = {-1, 0, 1};
      constexpr long shift_j[4] = {-2, -1, 0, 1};
      std::vector<Eigen::Triplet<double>> trip;
      vec xv(n), dxv(n);
      for (long a = 0; a < S; ++a)
        for (long ri = 0; ri < Pi; ++ri)
          for (long rj = 0; rj < Pj; ++rj) {
            std::fill(xv.begin(), xv.end(), 0.0);
            for (long i = ri; i < d0; i += Pi)
              for (long j = rj; j < d1; j += Pj) xv[idx(i, j, a)] = 1.0;
            (*this)(xv, dxv, 0.0);
            for (long i = ri; i < d0; i += Pi)
              for (long j = rj; j < d1; j += Pj) {
                const int col = idx(i, j, a);
                for (long di : shift_i) {
                  const long i2 = i - di;
                  if (i2 < 0 || i2 >= d0) continue;
                  for (long dj : shift_j) {
                    const long j2 = j - dj;
                    if (j2 < 0 || j2 >= d1) continue;
                    for (long b = 0; b < S; ++b) {
                      const int row = idx(i2, j2, b);
                      if (dxv[row] != 0.0) trip.emplace_back(row, col, dxv[row]);
                    }
                  }
                }
              }
          }
      Eigen::SparseMatrix<double> J(n, n);
      J.setFromTriplets(trip.begin(), trip.end());

      // safeguard: J w must reproduce the right-hand side
      Eigen::VectorXd w(n);
      double r = 0.37;
      for (int i = 0; i < n; ++i) { r = std::fmod(r * 7.31 + 0.113, 1.0); w[i] = 0.25 + r; }
      vec wv(w.data(), w.data() + n), fw(n);
      (*this)(wv, fw, 0.0);
      const Eigen::VectorXd Jw = J * w;
      double err = 0.0, nrm = 0.0;
      for (int i = 0; i < n; ++i) {
        err = std::max(err, std::fabs(Jw[i] - fw[i]));
        nrm = std::max(nrm, std::fabs(fw[i]));
      }
      if (!(err <= 1e-9 * nrm + 1e-300))
        throw std::runtime_error("daisie_odeint_iw: Jacobian structure not recognized");
      return J;
    }
  };

} // anonymous namespace

void daisie_iw_sparse_rosenbrock4(const std::function<void(const double*, double*)>& rhs,
  long d0, long d1, long S, vec& y, double t0, double t1, double atol, double rtol)
{
  iw_system sys{rhs, d0, d1, S};
  const int n = sys.size();
  const Eigen::SparseMatrix<double> J = sys.jacobian();
  Eigen::SparseMatrix<double> Id(n, n);
  Id.setIdentity();
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> lu;
  bool analyzed = false;
  Eigen::VectorXd b(n);
  auto prepare = [&](const vec&, double, double dt) {
    Eigen::SparseMatrix<double> M = Id * (1.0 / daisie_stiff::rosenbrock4_coefficients::gamma / dt) - J;
    M.makeCompressed();
    if (!analyzed) { lu.analyzePattern(M); analyzed = true; }
    lu.factorize(M);
    if (lu.info() != Eigen::Success)
      throw std::runtime_error("daisie_odeint_iw: sparse LU factorization failed");
  };
  auto solve = [&](vec& rhsv) {
    for (int i = 0; i < n; ++i) b[i] = rhsv[i];
    const Eigen::VectorXd sol = lu.solve(b);
    for (int i = 0; i < n; ++i) rhsv[i] = sol[i];
  };
  daisie_stiff::integrate_rosenbrock4_impl(sys, y, t0, t1, atol, rtol, prepare, solve);
}
