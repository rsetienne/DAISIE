// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]


//' @export daisie_odeint_iw


#define EIGEN_USE_THREADS
#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <RcppEigen.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <boost/numeric/odeint.hpp>
#include <utility>
#include <memory>
#include <functional>
#include <thread>


using namespace Rcpp;
using namespace Eigen;
using namespace boost::numeric::odeint;


namespace {

  using index_v = EIGEN_DEFAULT_DENSE_INDEX_TYPE;

  template <int Rank>
  using index_t = DSizes<index_v, Rank>;


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
  index_t<Rank> iofs(index_v i0, index_v i1, index_v i2);


  template <>
  index_t<3> iofs<3>(index_v i0, index_v i1, index_v i2)
  {
    return index_t<3>{i0, i1, i2};
  }


  template <>
  index_t<4> iofs<4>(index_v i0, index_v i1, index_v i2)
  {
    return index_t<4>{i0, i1, i2, 0};
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
    static tmap mapt(List pars, const char* name)
    {
      DoubleVector v = pars[name];
      auto dim = dim_to_index<Rank>(v);
      return tmap(v.begin(), dim);
    }
    static ctmap cmapt(List pars, const char* name)
    {
      DoubleVector v = pars[name];
      auto dim = dim_to_index<Rank>(v);
      return ctmap(v.begin(), dim);
    }

    double laa_;
    std::array<ctmap, 12> c_;
    matrix ki_;
  };


  template <int Rank>
  cpp_daisie_iw<Rank>::cpp_daisie_iw(List pars)
    : laa_(pars["laa"]),
      c_{
        cmapt(pars, "c1"),
        cmapt(pars, "c2"),
        cmapt(pars, "c3"),
        cmapt(pars, "c4"),
        cmapt(pars, "c5"),
        cmapt(pars, "c6"),
        cmapt(pars, "c7"),
        cmapt(pars, "c89"),
        cmapt(pars, "c10"),
        cmapt(pars, "c11"),
        cmapt(pars, "c12"),
        cmapt(pars, "c13")
      }
  {
    if (rank == 4) {
      DoubleVector ki = pars["ki"];
      auto dim = dim_to_index<2>(ki);
      ki_ = cmmap(ki.begin(), dim);
    }
  }


  template <>
  void cpp_daisie_iw<3>::rhs(const double* rx, double* rdx, ThreadPoolDevice* dev)
  {
    const auto dim_c = c_[0].dimensions();
    tmap dx(rdx, dim_c);
    std::array<std::pair<int, int>, rank> xxpad = {
      std::make_pair(1,1),
      std::make_pair(1,1),
      std::make_pair(2,1)
    };
    ctmap x(rx, dim_c);
    auto xx = x.pad(xxpad);
    auto xce = xx.slice(iofs<rank>(1,1,1), dim_c);
    auto x12 = xx.slice(iofs<rank>(1,1,2), dim_c);
    auto ddx =
      c_[0] * xx.slice(iofs<rank>(0,1,2), dim_c) +
      c_[1] * xx.slice(iofs<rank>(1,0,2), dim_c) +
      c_[2] * xx.slice(iofs<rank>(2,1,2), dim_c) +
      c_[3] * xx.slice(iofs<rank>(1,2,2), dim_c) +
      c_[4] * xx.slice(iofs<rank>(1,1,3), dim_c) +
      c_[5] * xx.slice(iofs<rank>(2,1,0), dim_c) +
      c_[6] * xx.slice(iofs<rank>(1,2,0), dim_c) +
      c_[7] * xce +
      c_[8] * xx.slice(iofs<rank>(2,1,1), dim_c) +
      c_[9] * xx.slice(iofs<rank>(1,2,1), dim_c) +
      c_[10] * x12;
    dx.device(*dev) = ddx;
  }


  template <>
  void cpp_daisie_iw<4>::rhs(const double* rx, double* rdx, ThreadPoolDevice* dev)
  {
    const auto dim_c = c_[0].dimensions();
    tmap dx(rdx, dim_c);
    std::array<std::pair<int, int>, rank> xxpad = {
      std::make_pair(1,1),
      std::make_pair(1,1),
      std::make_pair(2,1),
      std::make_pair(0,0)
    };
    ctmap x(rx, dim_c);
    auto xx = x.pad(xxpad);
    auto xce = xx.slice(iofs<rank>(1,1,1), dim_c);
    auto x12 = xx.slice(iofs<rank>(1,1,2), dim_c);
    auto ddx =
      c_[0] * xx.slice(iofs<rank>(0,1,2), dim_c) +
      c_[1] * xx.slice(iofs<rank>(1,0,2), dim_c) +
      c_[2] * xx.slice(iofs<rank>(2,1,2), dim_c) +
      c_[3] * xx.slice(iofs<rank>(1,2,2), dim_c) +
      c_[4] * xx.slice(iofs<rank>(1,1,3), dim_c) +
      c_[5] * xx.slice(iofs<rank>(2,1,0), dim_c) +
      c_[6] * xx.slice(iofs<rank>(1,2,0), dim_c) +
      c_[7] * xce +
      c_[8] * xx.slice(iofs<rank>(2,1,1), dim_c) +
      c_[9] * xx.slice(iofs<rank>(1,2,1), dim_c) +
      c_[10] * x12;
    const array<std::pair<int, int>, 1> product_dims = { std::make_pair(3, 1) };
    dx.device(*dev) = ddx + (laa_ * x12 + c_[11] * xce).contract(ki_, product_dims);
  }


  struct daisie_iw_wrapper
  {
    std::unique_ptr<ThreadPool> pool;
    std::unique_ptr<ThreadPoolDevice> dev;

    std::unique_ptr<cpp_daisie_iw<3>> iw3;
    std::unique_ptr<cpp_daisie_iw<4>> iw4;

    daisie_iw_wrapper(List pars)
    {
      pool.reset(new ThreadPool(std::thread::hardware_concurrency()));
      dev.reset(new ThreadPoolDevice(pool.get(), std::thread::hardware_concurrency()));

      int sysdim = pars["sysdim"];
      if (1 == sysdim) {
        iw3 = std::make_unique<cpp_daisie_iw<3>>(pars);
      }
      else {
        iw4 = std::make_unique<cpp_daisie_iw<4>>(pars);
      }
    }

    // odeint interface
    void operator()(const std::vector<double>& x, std::vector<double>& dxdt, double)
    {
      (iw3) ? iw3->rhs(x.data(), dxdt.data(), dev.get())
            : iw4->rhs(x.data(), dxdt.data(), dev.get());
    }
  };


  template <typename Stepper, typename IWrap>
  void integrate(double atol, double rtol, IWrap iw, std::vector<double>& y, double t0, double t1)
  {
    integrate_adaptive(make_controlled<Stepper>(atol, rtol), iw, y, t0, t1, 0.1 * (t1 - t0));
  }

}


//' Driver for the boost::odeint solver
//'
//' @name daisie_odeint_iw
RcppExport SEXP daisie_odeint_iw(SEXP ry, SEXP rtimes, SEXP rpars, SEXP Stepper, SEXP atolint, SEXP reltolint) {
  BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    auto y = as<std::vector<double>>(ry);
    auto times = as<std::vector<double>>(rtimes);
    auto pars = as<List>(rpars);
    auto stepper = as<std::string>(Stepper);
    auto atol = as<double>(atolint);
    auto rtol = as<double>(reltolint);

    daisie_iw_wrapper iw(pars);
    if ("odeint::runge_kutta_cash_karp54" == stepper) {
      integrate<runge_kutta_cash_karp54<std::vector<double>>>(atol, rtol, std::ref(iw), y, times[0], times[1]);
    }
    else if ("odeint::runge_kutta_fehlberg78" == stepper) {
      integrate<runge_kutta_fehlberg78<std::vector<double>>>(atol, rtol, std::ref(iw), y, times[0], times[1]);
    }
    else if ("odeint::runge_kutta_dorpi5" == stepper) {
      integrate<runge_kutta_dopri5<std::vector<double>>>(atol, rtol, std::ref(iw), y, times[0], times[1]);
    }
    else if ("odeint::bulirsch_stoer" == stepper) {
      using stepper_t = bulirsch_stoer<std::vector<double>>;
      integrate_adaptive(stepper_t(atol, rtol), std::ref(iw), y, times[0], times[1], 0.1 * (times[1] - times[0]));
    }
    else {
      throw std::runtime_error("daisie_odeint_iw: unknown stepper");
    }
    rcpp_result_gen = y;
    return rcpp_result_gen;
  END_RCPP
}
