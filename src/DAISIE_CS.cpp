// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(BH)]]


//' @export daisie_odeint_cs


#include "DAISIE_odeint.h"


namespace {


  // maximal number of steps the solver is executing.
  // prevents odeint from getting stuckle
  // at-hoc - 'solution'.
  static constexpr int default_max_cs_steps = 1'000'000;
  static int max_cs_steps = default_max_cs_steps;

  // step-size factor for adams_bashforth_moulton integration
  static constexpr double default_abm_factor = 0.0001;
  static double abm_factor = default_abm_factor;

  //
  class padded_vector_view
  {
  public:
    padded_vector_view(int pad, const double* data, int n) :
    data_(data), pad_(pad), n_(n)
    {
    }

    // return 0.0 for indices 'i' outside [pad, pad + n)
    double operator[](int i) const
    {
      const auto ii = i - pad_;
      return (ii >= 0 && i < n_) ? *(data_ + ii) : 0.0;
    }

  private:
    const double* data_ = nullptr;  // this->operator[pad_] == *data_
    int pad_ = 0;
    int n_ = 0;
  };


  // common parameter
  struct param_t
  {
    param_t(int LX, int KK, state_type&& par) :
      lx(LX), kk(KK), P(std::move(par))
    {
    }

    int lx;
    int kk;
    state_type P;
    mutable int steps = 0;
  };


  class cpp_daisie_cs_runmod
  {
  public:
    cpp_daisie_cs_runmod(param_t&& par) :
      p_(std::move(par))
    {
    }

    // odeint interface
    void operator()(const state_type& x, state_type& dx, double) const
    {
      if (++p_.steps > max_cs_steps) throw std::runtime_error("cpp_daisie_cs_runmod: too many steps");

      // using padded views instead of vectors:
      const auto xx1 = padded_vector_view(2, x.data().begin(), p_.lx);
      const auto xx2 = padded_vector_view(2, x.data().begin() + p_.lx, p_.lx);
      const auto xx3 = x[2 * p_.lx];

      // using views instead of vectors:
      const auto chunk = p_.lx + 4 + 2 * p_.kk;
      const auto laavec = p_.P.data().begin();
      const auto lacvec = p_.P.data().begin() + chunk;
      const auto muvec = p_.P.data().begin() + 2 * chunk;
      const auto gamvec = p_.P.data().begin() + 3 * chunk;
      const auto nn = p_.P.data().begin() + 4 * chunk;

      // using offsets into our views instead of vectors:
      const int il1 = 2 + p_.kk - 1;
      const int il2 = 2 + p_.kk + 1;
      const int il3in3 = 2 + p_.kk;
      const int il4 = 2 + p_.kk - 2;
      const int in1 = 2 + 2 * p_.kk - 1;
      const int in2ix2 = 2 + 1;
      const int ix1 = 2 - 1;
      const int ix3 = 2;
      const int ix4 = 2 - 2;

      // using views into output vector:
      auto dx1 = dx.data().begin();
      auto dx2 = dx1 + p_.lx;
      for (int i = 0; i < p_.lx; ++i) {
        dx1[i] = laavec[il1 + i + 1] * xx2[ix1 + i]
               + lacvec[il4 + i + 1] * xx2[ix4 + i]
               + muvec[il2 + i + 1] * xx2[ix3 + i]
               + lacvec[il1 + i] * nn[in1 + i] * xx1[ix1 + i]
               + muvec[il2 + i] * nn[in2ix2 + i] * xx1[in2ix2 + i]
               - (muvec[il3in3 + i] + lacvec[il3in3 + i]) * nn[il3in3 + i] * xx1[ix3 + i]
               - gamvec[il3in3 + i] * xx1[ix3 + i];
        dx2[i] = gamvec[il3in3 + i] * xx1[ix3 + i]
               + lacvec[il1 + i + 1] * nn[in1 + i] * xx2[ix1 + i]
               + muvec[il2 + i + 1] * nn[in2ix2 + i] * xx2[in2ix2 + i]
               - (muvec[il3in3 + i + 1] + lacvec[il3in3 + i + 1]) * nn[il3in3 + i + 1] * xx2[ix3 + i]
               - laavec[il3in3 + i + 1] * xx2[ix3 + i];
      }

      auto dx3 = dx2 + p_.lx;
      dx3[0] = 0;
    }

  private:
    const param_t p_;
  };

  class cpp_daisie_cs_runmod_1
  {
  public:
    cpp_daisie_cs_runmod_1(param_t&& p) :
    p_(p)
    {
    }

    // odeint interface
    void operator()(const state_type& x, state_type& dx, double) const
    {
      if (++p_.steps > max_cs_steps) throw std::runtime_error("cpp_daisie_cs_runmod_1: too many steps");

      // using padded views instead of vectors:
      const auto xx1 = padded_vector_view(2, x.data().begin(), p_.lx);
      const auto xx2 = padded_vector_view(2, x.data().begin() + p_.lx, p_.lx);
      const auto xx3 = padded_vector_view(2, x.data().begin() + 2 * p_.lx, p_.lx);
      const auto xx4 = padded_vector_view(2, x.data().begin() + 3 * p_.lx, p_.lx);

      // using views instead of vectors:
      const auto chunk = p_.lx + 4 + 2 * p_.kk;
      const auto laavec = p_.P.data().begin();
      const auto lacvec = p_.P.data().begin() + chunk;
      const auto muvec = p_.P.data().begin() + 2 * chunk;
      const auto gamvec = p_.P.data().begin() + 3 * chunk;
      const auto nn = p_.P.data().begin() + 4 * chunk;

      // using offsets into our views instead of vectors:
      const int il1 = 2 + p_.kk - 1;
      const int il2 = 2 + p_.kk + 1;
      const int il3in3 = 2 + p_.kk;
      const int il4 = 2 + p_.kk - 2;
      const int in1 = 2 + 2 * p_.kk - 1;
      const int in2ix2 = 2 + 1;   // spilt in in2, ix2 at no cost?
      const int in4ix1 = 2 - 1;   // split in in4, ix1 at no cost?
      const int ix3 = 2;
      const int ix4 = 2 - 2;

      // using views into output vector:
      auto dx1 = dx.data().begin();
      auto dx2 = dx1 + p_.lx;
      auto dx3 = dx2 + p_.lx;
      auto dx4 = dx3 + p_.lx;
      for (int i = 0; i < p_.lx; ++i) {
        dx1[i] = lacvec[il1 + i] * nn[in1 + i] * xx1[in4ix1 + i]
        + laavec[il1 + i + 1] * xx2[in4ix1 + i]
        + lacvec[il4 + i + 1] * xx2[ix4 + i]
        + muvec[il2 + i] * nn[in2ix2 + i] * xx1[in2ix2 + i]
        + muvec[il3in3 + i + 1] * xx2[ix3 + i]
        - (muvec[il3in3 + i] + lacvec[il3in3 + i]) * nn[il3in3 + i] * xx1[ix3 + i]
        - gamvec[il3in3 + i] * xx1[ix3 + i];
        dx2[i] = gamvec[il3in3 + i] * xx1[ix3 + i]
        + gamvec[il3in3 + i] * xx3[ix3 + i]
        + gamvec[il3in3 + i + 1] * xx4[ix3 + i]
        + lacvec[il1 + i + 1] * nn[in1 + i] * xx2[in4ix1 + i]
        + muvec[il2 + i + 1] * nn[in2ix2 + i] * xx2[in2ix2 + i]
        - (muvec[il3in3 + i + 1] + lacvec[il3in3 + i + 1]) * nn[il3in3 + i + 1] * xx2[ix3 + i]
        - laavec[il3in3 + i + 1] * xx2[ix3 + i];
        dx3[i] = lacvec[il1 + i] * nn[in1 + i] * xx3[in4ix1 + i]
        + laavec[il1 + i + 1] * xx4[in4ix1 + i]
        + lacvec[il4 + i + 1] * xx4[ix4 + i]
        + muvec[il2 + i] * nn[in2ix2 + i] * xx3[in2ix2 + i]
        + muvec[il3in3 + i + 1] * xx4[ix3 + i]
        - (lacvec[il3in3 + i] + muvec[il3in3 + i]) * nn[il3in3 + i] * xx3[ix3 + i]
        - gamvec[il3in3 + i] * xx3[ix3 + i];
        dx4[i] = lacvec[il1 + i + 1] * nn[in1 + i] * xx4[in4ix1 + i]
        + muvec[il2 + i + 1] * nn[in2ix2 + i] * xx4[in2ix2 + i]
        - (lacvec[il3in3 + i + 1] + muvec[il3in3 + i + 1]) * nn[il3in3 + i + 1] * xx4[ix3 + i]
        - gamvec[il3in3 + i + 1] * xx4[ix3 + i];
      }
    }

  private:
    const param_t p_;
  };


  class cpp_daisie_cs_runmod_2
  {
  public:
    cpp_daisie_cs_runmod_2(param_t&& p) :
      p_(p)
    {
    }

    // odeint interface
    void operator()(const state_type& x, state_type& dx, double) const
    {
      if (++p_.steps > max_cs_steps) throw std::runtime_error("cpp_daisie_cs_runmod_2: too many steps");

      // using padded views instead of vectors:
      const auto xx1 = padded_vector_view(2, x.data().begin(), p_.lx);
      const auto xx2 = padded_vector_view(2, x.data().begin() + p_.lx, p_.lx);
      const auto xx3 = padded_vector_view(2, x.data().begin() + 2 * p_.lx, p_.lx);

      // using views instead of vectors:
      const auto chunk = p_.lx + 4 + 2 * p_.kk;
      const auto laavec = p_.P.data().begin();
      const auto lacvec = p_.P.data().begin() + chunk;
      const auto muvec = p_.P.data().begin() + 2 * chunk;
      const auto gamvec = p_.P.data().begin() + 3 * chunk;
      const auto nn = p_.P.data().begin() + 4 * chunk;

      // using offsets into our views instead of vectors:
      const int il1 = 2 + p_.kk - 1;
      const int il2 = 2 + p_.kk + 1;
      const int il3in3 = 2 + p_.kk;
      const int il4 = 2 + p_.kk - 2;
      const int in1 = 2 + 2 * p_.kk - 1;
      const int in2ix2 = 2 + 1;   // spilt in in2, ix2 at no cost?
      const int in4ix1 = 2 - 1;   // split in in4, ix1 at no cost?
      const int ix3 = 2;
      const int ix4 = 2 - 2;

      // using views into output vector:
      auto dx1 = dx.data().begin();
      auto dx2 = dx1 + p_.lx;
      auto dx3 = dx2 + p_.lx;
      for (int i = 0; i < p_.lx; ++i) {
        dx1[i] = laavec[il1 + i + 1] * xx2[in4ix1 + i]
               + lacvec[il4 + i + 1] * xx2[ix4 + i]
               + muvec[il2 + i + 1] * xx2[ix3 + i]
               + lacvec[il1 + i] * nn[in1 + i] * xx1[in4ix1 + i]
               + muvec[il2 + i] * nn[in2ix2 + i] * xx1[in2ix2 + i]
               - (muvec[il3in3 + i] + lacvec[il3in3 + i]) * nn[il3in3 + i] * xx1[ix3 + i]
               - gamvec[il3in3 + i] * xx1[ix3 + i];
        if (1 == p_.kk) {
          dx1[i] += laavec[il3in3 + i] * xx3[ix3 + i] + 2.0 * lacvec[il1 + i] * xx3[in4ix1 + i];
        }
        dx2[i] = gamvec[il3in3 + i] * xx1[ix3 + i]
               + lacvec[il1 + i + 1] * nn[in1 + i] * xx2[in4ix1 + i]
               + muvec[il2 + i + 1] * nn[in2ix2 + i] * xx2[in2ix2 + i]
               - (muvec[il3in3 + i + 1] + lacvec[il3in3 + i + 1]) * nn[il3in3 + i + 1] * xx2[ix3 + i]
               - laavec[il3in3 + i + 1] * xx2[ix3 + i];
        dx3[i] = lacvec[il1 + i] * nn[in4ix1 + i] * xx3[in4ix1 + i]
               + muvec[il2 + i] * nn[in2ix2 + i] * xx3[in2ix2 + i]
               - (lacvec[il3in3 + i] + muvec[il3in3 + i]) * nn[il3in3 + i] * xx3[ix3 + i]
               - (laavec[il3in3 + i] + gamvec[il3in3 + i]) * xx3[ix3 + i];
      }
    }

  private:
    const param_t p_;
  };

} // anonymous namespace


//' Driver for the boost::odeint solver
//'
//' @name daisie_odeint_cs
RcppExport SEXP daisie_odeint_cs(SEXP rrunmod, SEXP ry, SEXP rtimes, SEXP rlx, SEXP rkk, SEXP rpar, SEXP Stepper, SEXP atolint, SEXP reltolint) {
BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  auto runmod = as<std::string>(rrunmod);
  auto y = as<state_type>(ry);
  auto times = as<std::vector<double>>(rtimes);
  auto lx = as<int>(rlx);
  auto kk = as<int>(rkk);
  auto stepper = as<std::string>(Stepper);
  auto atol = as<double>(atolint);
  auto rtol = as<double>(reltolint);

  auto p = param_t(lx, kk, as<state_type>(rpar));
  if (runmod == "daisie_runmod") {
    cpp_daisie_cs_runmod rhs(std::move(p));
    daisie_odeint::integrate(stepper, std::ref(rhs), y, times[0], times[1], atol, rtol);
  }
  else if (runmod == "daisie_runmod1") {
    cpp_daisie_cs_runmod_1 rhs(std::move(p));
    daisie_odeint::integrate(stepper, std::ref(rhs), y, times[0], times[1], atol, rtol);
  }
  else if (runmod == "daisie_runmod2") {
    cpp_daisie_cs_runmod_2 rhs(std::move(p));
    daisie_odeint::integrate(stepper, std::ref(rhs), y, times[0], times[1], atol, rtol);
  }
  else {
    throw std::runtime_error("daisie_odeint_cs: unknown runmod");
  }

  rcpp_result_gen = y;
  return rcpp_result_gen;
END_RCPP
}


RcppExport SEXP daisie_odeint_cs_max_steps(SEXP rmax_steps) {
  BEGIN_RCPP
  max_cs_steps = (0 < as<int>(rmax_steps)) ? as<int>(rmax_steps) : default_max_cs_steps;
  return wrap(max_cs_steps);
  END_RCPP
}


namespace daisie_odeint {

  // step-size factor for adams_bashforth_moulton integration
  constexpr double default_abm_factor = 0.0001;
  double abm_factor = default_abm_factor;

}


// misplaced
RcppExport SEXP daisie_odeint_abm_factor(SEXP rfactor) {
  BEGIN_RCPP
  daisie_odeint::abm_factor = (0 < as<double>(rfactor)) ? as<double>(rfactor) : daisie_odeint::default_abm_factor;
  return wrap(daisie_odeint::abm_factor);
  END_RCPP
}

