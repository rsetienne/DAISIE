//
//  Copyright (c) 2023, Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

//' @export daisie_odeint_iw

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]

#include "config.h"
#include "DAISIE_odeint.h"
#define EIGEN_USE_THREADS
#include <RcppEigen.h>
#include <Eigen/SparseLU>
#include <unsupported/Eigen/CXX11/Tensor>
#include <utility>
#include <array>
#include <memory>
#include <functional>
#include <thread>


using namespace Eigen;


// num_threads
unsigned daisie_odeint_iw_num_threads_ = std::max(1u, std::thread::hardware_concurrency());
using namespace daisie_odeint::jacobian_policy;


namespace {

  using index_v = EIGEN_DEFAULT_DENSE_INDEX_TYPE;


  template <int I>
  using index_t = std::array<index_v, I>;


  template <index_v I, index_v J = I>
  using index_array_t = std::array<index_t<I>, J>;


  template <int I>
  using padding_t = std::array<std::pair<index_v, index_v>, I>;


  template <int Rank>
  index_t<Rank> dim_to_index(DoubleVector v)
  {
    IntegerVector iv = v.attr("dim");
    auto ret = index_t<Rank>{};
    for (size_t i = 0; i < Rank; ++i) {
      ret[i] = iv[i];
    }
    return ret;
  }


  template <int Rank>
  index_t<Rank> iofs(int i0, int i1);


  template <>
  index_t<2> iofs<2>(int i0, int i1)
  {
    return index_t<2>{i0, i1};
  }


  template <>
  index_t<3> iofs<3>(int i0, int i1)
  {
    return index_t<3>{i0, i1, 0};
  }


  template <int Rank>
  class cpp_daisie_iw
  {
    static constexpr int rank = Rank;
    using index = index_t<Rank>;
    using tensor = Tensor<double, Rank>;
    using tmap = TensorMap<tensor>;
    using ctensor = const Tensor<const double, Rank>;
    using ctmap = TensorMap<ctensor>;
    using matrix = Tensor<double, 2>;
    using cmatrix = Tensor<const double, 2>;
    using cmmap = TensorMap<cmatrix>;

  public:
    explicit cpp_daisie_iw(List pars);
    void rhs(const double* x, double* rdx, ThreadPoolDevice* dev);

    // extents of the state tensor (third extent is 1 for rank 2)
    std::array<index_v, 3> extents() const
    {
      const auto d = c_[0].dimensions();
      return { d[0], d[1], (Rank > 2) ? d[2] : index_v(1) };
    }

  private:
    static ctmap cmapt(List pars, const char* name)
    {
      DoubleVector v = pars[name];
      auto dim = dim_to_index<Rank>(v);
      return ctmap(v.begin(), dim);
    }

    double laa_;
    std::array<ctmap, 8> c_;
    matrix ki_;
  };


  template <int Rank>
  cpp_daisie_iw<Rank>::cpp_daisie_iw(List pars) :
    laa_(pars["laa"]),
    c_{
      cmapt(pars, "c1"),
      cmapt(pars, "c2"),
      cmapt(pars, "c3"),
      cmapt(pars, "c4"),
      cmapt(pars, "c5"),
      cmapt(pars, "c6"),
      cmapt(pars, "c7"),
      cmapt(pars, "c8")
    }
  {
    if (rank > 2) {
      DoubleVector ki = pars["ki"];
      auto dim = dim_to_index<2>(ki);
      ki_ = cmmap(ki.begin(), dim);
    }
  }


# define xx_slice(x, y) xx.slice(iofs<rank>(x,y), dim_c)


  template <>
  void cpp_daisie_iw<2>::rhs(const double* rx, double* rdx, ThreadPoolDevice* dev)
  {
    const auto dim_c = c_[0].dimensions();
    tmap dx(rdx, dim_c);
    ctmap x(rx, dim_c);
    auto xxpad = padding_t<rank>{{ {1,1}, {2,1} }};
    const auto xx = x.pad(xxpad);
    auto ddx =
      c_[0] * xx_slice(0,2) +
      c_[1] * xx_slice(2,2) +
      c_[2] * xx_slice(1,3) +
      c_[3] * xx_slice(2,1) +
      c_[4] * xx_slice(2,0) +
      c_[5] * xx_slice(1,1) -
      c_[6] * xx_slice(1,2);
    if (dev) {
      dx.device(*dev) = ddx;
    }
    else {
      dx = ddx;
    }
  }


  template <>
  void cpp_daisie_iw<3>::rhs(const double* rx, double* rdx, ThreadPoolDevice* dev)
  {
    const auto dim_c = c_[0].dimensions();
    tmap dx(rdx, dim_c);
    ctmap x(rx, dim_c);
    const auto product_dims = padding_t<1>{{ {2,1} }};
    auto xxpad = padding_t<rank>{{ {1,1}, {2,1}, {0,0} }};
    const auto xx = x.pad(xxpad);
    auto ddx =
      c_[0] * xx_slice(0,2) +
      c_[1] * xx_slice(2,2) +
      c_[2] * xx_slice(1,3) +
      c_[3] * xx_slice(2,1) +
      c_[4] * xx_slice(2,0) +
      c_[5] * xx_slice(1,1) -
      c_[6] * xx_slice(1,2) +
      (laa_ * xx_slice(1,2) + c_[7] * xx_slice(1,1)).contract(ki_, product_dims);
    if (dev) {
      dx.device(*dev) = ddx;
    }
    else {
      dx = ddx;
    }
  }


  struct daisie_iw_wrapper
  {
    using jacobian = const_from_linear_rhs<daisie_iw_wrapper>;

    std::unique_ptr<ThreadPool> pool;
    std::unique_ptr<ThreadPoolDevice> dev;

    std::unique_ptr<cpp_daisie_iw<2>> iw2;
    std::unique_ptr<cpp_daisie_iw<3>> iw3;

    daisie_iw_wrapper(List pars)
    {
      if (1 != daisie_odeint_iw_num_threads_) {
        pool.reset(new ThreadPool(daisie_odeint_iw_num_threads_));
        dev.reset(new ThreadPoolDevice(pool.get(), daisie_odeint_iw_num_threads_));
      }
      int sysdim = pars["sysdim"];
      if (1 == sysdim) {
        iw2 = std::make_unique<cpp_daisie_iw<2>>(pars);
      }
      else {
        iw3 = std::make_unique<cpp_daisie_iw<3>>(pars);
      }
    }

    static constexpr bool stiff_linear = true;
    static constexpr bool stiff_time_dependent = false;

    std::array<index_v, 3> extents() const
    {
      return iw2 ? iw2->extents() : iw3->extents();
    }

    int size() const
    {
      const auto e = extents();
      return static_cast<int>(e[0] * e[1] * e[2]);
    }

    // Jacobian of the (linear) system as a sparse matrix.
    //
    // The state is a tensor x(i, j, a) (column major). The right-hand side
    // couples the output cell (i, j, .) with the input cells
    // (i + di, j + dj, .), di in {-1, 0, 1}, dj in {-2, -1, 0, 1} (the shifts
    // of the padded slices in rhs()), and, through the contraction with ki_,
    // all slices a, b of the same cells. Two columns whose cell indices differ
    // by multiples of (3, 4) in (i, j) therefore never write into the same
    // row, and the 12 * (number of slices) probes with sums of such columns
    // recover the whole Jacobian.
    Eigen::SparseMatrix<double> sparse_jacobian()
    {
      const auto e = extents();
      const index_v d0 = e[0], d1 = e[1], S = e[2];
      const int n = size();
      auto idx = [&](index_v i, index_v j, index_v a) {
        return static_cast<int>(i + d0 * (j + d1 * a));
      };
      constexpr index_v Pi = 3, Pj = 4;
      constexpr index_v shift_i[3] = {-1, 0, 1};
      constexpr index_v shift_j[4] = {-2, -1, 0, 1};

      std::vector<Eigen::Triplet<double>> trip;
      state_type xv(n), dxv(n);
      for (index_v a = 0; a < S; ++a) {
        for (index_v ri = 0; ri < Pi; ++ri) {
          for (index_v rj = 0; rj < Pj; ++rj) {
            std::fill(xv.begin(), xv.end(), 0.0);
            for (index_v i = ri; i < d0; i += Pi)
              for (index_v j = rj; j < d1; j += Pj) xv[idx(i, j, a)] = 1.0;
            (*this)(xv, dxv, 0.0);
            for (index_v i = ri; i < d0; i += Pi) {
              for (index_v j = rj; j < d1; j += Pj) {
                const int col = idx(i, j, a);
                for (index_v di : shift_i) {
                  const index_v i2 = i - di;
                  if (i2 < 0 || i2 >= d0) continue;
                  for (index_v dj : shift_j) {
                    const index_v j2 = j - dj;
                    if (j2 < 0 || j2 >= d1) continue;
                    for (index_v b = 0; b < S; ++b) {
                      const int row = idx(i2, j2, b);
                      const double v = dxv[row];
                      if (v != 0.0) trip.emplace_back(row, col, v);
                    }
                  }
                }
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
      state_type wv(n), fw(n);
      for (int i = 0; i < n; ++i) wv[i] = w[i];
      (*this)(wv, fw, 0.0);
      Eigen::VectorXd Jw = J * w;
      double err = 0.0, nrm = 0.0;
      for (int i = 0; i < n; ++i) {
        err = std::max(err, std::fabs(Jw[i] - fw[i]));
        nrm = std::max(nrm, std::fabs(fw[i]));
      }
      if (!(err <= 1e-9 * nrm + 1e-300)) {
        throw std::runtime_error("daisie_odeint_iw: Jacobian structure not recognized");
      }
      return J;
    }

    // Rosenbrock integration with a sparse LU (Eigen, COLAMD ordering).
    void integrate_rosenbrock4(state_type& y, double t0, double t1, double atol, double rtol)
    {
      const int n = size();
      const Eigen::SparseMatrix<double> J = sparse_jacobian();
      Eigen::SparseMatrix<double> Id(n, n);
      Id.setIdentity();
      Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> lu;
      bool analyzed = false;
      Eigen::VectorXd rhs(n);
      auto prepare = [&](const state_type&, double, double dt) {
        Eigen::SparseMatrix<double> M = Id * (1.0 / daisie_stiff::rosenbrock4_coefficients::gamma / dt) - J;
        M.makeCompressed();
        if (!analyzed) { lu.analyzePattern(M); analyzed = true; }
        lu.factorize(M);
        if (lu.info() != Eigen::Success) {
          throw std::runtime_error("daisie_odeint_iw: sparse LU factorization failed");
        }
      };
      auto solve = [&](state_type& b) {
        for (int i = 0; i < n; ++i) rhs[i] = b[i];
        const Eigen::VectorXd sol = lu.solve(rhs);
        for (int i = 0; i < n; ++i) b[i] = sol[i];
      };
      daisie_stiff::integrate_rosenbrock4_impl(*this, y, t0, t1, atol, rtol, prepare, solve);
    }

    // right-hand side, used for Jacobian probing by the sparse Rosenbrock
    void eval(const state_type& x, state_type& dxdt, double t)
    {
      (*this)(x, dxdt, t);
    }

    // odeint interface
    void operator()(const state_type& x, state_type& dxdt, double)
    {
      (iw2) ? iw2->rhs(x.data().begin(), dxdt.data().begin(), dev.get())
            : iw3->rhs(x.data().begin(), dxdt.data().begin(), dev.get());
    }
  };

} // anonymous namespace


//' Driver for the boost::odeint solver for the IW model
//'
//' @name daisie_odeint_iw
//' @export daisie_odeint_iw
//' @return Object of type `state_type`, which itself is
//' `vector_t`<double>, with the result of the
//' integration depending on the runmod chosen.
RcppExport SEXP daisie_odeint_iw(SEXP ry, SEXP rtimes, SEXP rpars, SEXP Stepper, SEXP atolint, SEXP reltolint) {
BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  auto y = as<state_type>(ry);
  auto times = as<std::vector<double>>(rtimes);
  auto pars = as<List>(rpars);
  auto stepper = as<std::string>(Stepper);
  auto atol = as<double>(atolint);
  auto rtol = as<double>(reltolint);

  daisie_iw_wrapper iw(pars);
  if ("odeint::rosenbrock4" == stepper) {
    iw.integrate_rosenbrock4(y, times[0], times[1], atol, rtol);
  } else {
    daisie_odeint::integrate(stepper, std::ref(iw), y, times[0], times[1], atol, rtol);
  }

  rcpp_result_gen = y;
  return rcpp_result_gen;
END_RCPP
}


RcppExport SEXP daisie_odeint_iw_num_threads(SEXP rnum_threads) {
BEGIN_RCPP
  auto num_threads = as<int>(rnum_threads);
  if (0 <= num_threads) {
    daisie_odeint_iw_num_threads_ = (0 == num_threads)
      ? std::max(1u, std::thread::hardware_concurrency())
      : std::max(1u, std::min(static_cast<unsigned>(num_threads), std::thread::hardware_concurrency()));
  }
  return wrap(daisie_odeint_iw_num_threads_);
END_RCPP
}
