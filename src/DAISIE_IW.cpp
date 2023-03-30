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
  daisie_odeint::integrate(stepper, std::ref(iw), y, times[0], times[1], atol, rtol);

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
