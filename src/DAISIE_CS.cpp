//
//  Copyright (c) 2023, Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

// [[Rcpp::depends(BH)]]



#include "config.h"
#include "DAISIE_odeint.h"

using namespace daisie_odeint::jacobian_policy;


namespace {


  // maximal number of steps the solver is executing.
  // prevents odeint from getting stuckle
  // at-hoc - 'solution'.
  static constexpr int default_max_cs_steps = 1000000;
  static int max_cs_steps = default_max_cs_steps;


  // step-size factor for adams_bashforth_moulton integration
  static constexpr double default_abm_factor = 0.0001;
  double abm_factor = default_abm_factor;


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
    using jacobian = const_from_linear_rhs<cpp_daisie_cs_runmod>;

    cpp_daisie_cs_runmod(param_t&& par) :
      p_(std::move(par))
    {
    }

    // odeint interface
    void operator()(const state_type& x, state_type& dx, double /*t*/) const
    {
      if (++p_.steps > max_cs_steps) {
        throw std::runtime_error("cpp_daisie_cs_runmod: too many steps");
      }

      const auto xx1 = padded_vector_view<2>(x.data().begin(), p_.lx);
      const auto xx2 = padded_vector_view<2>(x.data().begin() + p_.lx, p_.lx);
      const auto chunk = p_.lx + 4 + 2 * p_.kk;
      const auto laavec = p_.P.data().begin();
      const auto lacvec = p_.P.data().begin() + chunk;
      const auto muvec = p_.P.data().begin() + 2 * chunk;
      const auto gamvec = p_.P.data().begin() + 3 * chunk;
      const auto nn = p_.P.data().begin() + 4 * chunk;

      // using offsets into our views instead of vectors:
      const int nil2lx = 2;
      const int il1 = nil2lx + p_.kk - 1;
      const int il2 = nil2lx + p_.kk + 1;
      const int il3 = nil2lx + p_.kk;
      const int il4 = nil2lx + p_.kk - 2;

      const int in1 = nil2lx + 2 * p_.kk - 1;
      const int in2 = nil2lx + 1;
      const int in3 = nil2lx + p_.kk;

      const int ix1 = nil2lx - 1;
      const int ix2 = nil2lx + 1;
      const int ix3 = nil2lx;
      const int ix4 = nil2lx - 2;

      auto dx1 = dx.data().begin();
      auto dx2 = dx1 + p_.lx;
      auto dx3 = dx2 + p_.lx;

      for (int i = 0; i < p_.lx; ++i) {
        dx1[i] = laavec[il1 + i + 1] * xx2[ix1 + i]
               + lacvec[il4 + i + 1] * xx2[ix4 + i]
               + muvec[il2 + i + 1] * xx2[ix3 + i]
               + lacvec[il1 + i] * nn[in1 + i] * xx1[ix1 + i]
               + muvec[il2 + i] * nn[in2 + i] * xx1[ix2 + i]
               - (muvec[il3 + i] + lacvec[il3 + i]) * nn[in3 + i] * xx1[ix3 + i]
               - gamvec[il3 + i] * xx1[ix3 + i];
        dx2[i] = gamvec[il3 + i] * xx1[ix3 + i]
               + lacvec[il1 + i + 1] * nn[in1 + i] * xx2[ix1 + i]
               + muvec[il2 + i + 1] * nn[in2 + i] * xx2[ix2 + i]
               - (muvec[il3 + 1 + i] + lacvec[il3 + 1 + i]) * nn[in3 + i + 1] * xx2[ix3 + i]
               - laavec[il3 + i] * xx2[ix3 + i];
      }
      dx3[0] = 0.0;
    }

  private:
    const param_t p_;
  };

  class cpp_daisie_cs_runmod_1
  {
  public:
    using jacobian = const_from_linear_rhs<cpp_daisie_cs_runmod_1>;

    cpp_daisie_cs_runmod_1(param_t&& par) :
      p_(std::move(par))
    {
    }

    // odeint interface
    void operator()(const state_type& x, state_type& dx, double /*t*/) const
    {
      if (++p_.steps > max_cs_steps) {
        throw std::runtime_error("cpp_daisie_cs_runmod_1: too many steps");
      }

      const auto xx1 = padded_vector_view<2>(x.data().begin(), p_.lx);
      const auto xx2 = padded_vector_view<2>(x.data().begin() + p_.lx, p_.lx);
      const auto xx3 = padded_vector_view<2>(x.data().begin() + 2 * p_.lx, p_.lx);
      const auto xx4 = padded_vector_view<2>(x.data().begin() + 3 * p_.lx, p_.lx);

      const auto chunk = p_.lx + 4 + 2 * p_.kk;
      const auto laavec = p_.P.data().begin();
      const auto lacvec = p_.P.data().begin() + chunk;
      const auto muvec = p_.P.data().begin() + 2 * chunk;
      const auto gamvec = p_.P.data().begin() + 3 * chunk;
      const auto nn = p_.P.data().begin() + 4 * chunk;

      // using offsets into our views instead of vectors:
      const int nil2lx = 2;
      const int il1 = nil2lx + p_.kk - 1;
      const int il2 = nil2lx + p_.kk + 1;
      const int il3 = nil2lx + p_.kk;
      const int il4 = nil2lx + p_.kk - 2;

      const int in1 = nil2lx + 2 * p_.kk - 1;
      const int in2 = nil2lx + 1;
      const int in3 = nil2lx + p_.kk;

      const int ix1 = nil2lx - 1;
      const int ix2 = nil2lx + 1;
      const int ix3 = nil2lx;
      const int ix4 = nil2lx - 2;

      // using views into output vector:
      auto dx1 = dx.data().begin();
      auto dx2 = dx1 + p_.lx;
      auto dx3 = dx2 + p_.lx;
      auto dx4 = dx3 + p_.lx;

      for (int i = 0; i < p_.lx; ++i) {
        dx1[i] = lacvec[il1 + i] * nn[in1 + i] * xx1[ix1 + i]
               + laavec[il1 + i + 1] * xx2[ix1 + i]
               + lacvec[il4 + i + 1] * xx2[ix4 + i]
               + muvec[il2 + i] * nn[in2 + i] * xx1[ix2 + i]
               + muvec[il3 + i + 1] * xx2[ix3 + i]
               - (muvec[il3 + i] + lacvec[il3 + i]) * nn[in3 + i] * xx1[ix3 + i]
               - gamvec[il3 + i] * xx1[ix3 + i];
        dx2[i] = gamvec[il3 + i] * xx1[ix3 + i]
               + gamvec[il3 + i] * xx3[ix3 + i]
               + gamvec[il3 + i + 1] * xx4[ix3 + i]
               + lacvec[il1 + i + 1] * nn[in1 + i] * xx2[ix1 + i]
               + muvec[il2 + i + 1] * nn[in2 + i] * xx2[ix2 + i]
               - (muvec[il3 + 1 + i] + lacvec[il3 + 1 + i]) * nn[in3 + i + 1] * xx2[ix3 + i]
               - laavec[il3 + i] * xx2[ix3 + i];
        dx3[i] = lacvec[il1 + i] * nn[in1 + i] * xx3[ix1 + i]
               + laavec[il1 + i + 1] * xx4[ix1 + i]
               + lacvec[il4 + i + 1] * xx4[ix4 + i]
               + muvec[il2 + i] * nn[in2 + i] * xx3[ix2 + i]
               + muvec[il3 + i + 1] * xx4[ix3 + i]
               - (lacvec[il3 + i] + muvec[il3 + i]) * nn[in3 + i] * xx3[ix3 + i]
               - gamvec[il3 + i] * xx3[ix3 + i];
        dx4[i] = lacvec[il1 + i + 1] * nn[in1 + i] * xx4[ix1 + i]
               + muvec[il2 + i + 1] * nn[in2 + i] * xx4[ix2 + i]
               - (lacvec[il3 + i + 1] + muvec[il3 + i + 1]) * nn[in3 + i + 1] * xx4[ix3 + i]
               - laavec[il3 + i + 1] * xx4[ix3 + i]
               - gamvec[il3 + i + 1] * xx4[ix3 + i];
      }
    }

  private:
    const param_t p_;
  };


  class cpp_daisie_cs_runmod_2
  {
  public:
    using jacobian = const_from_linear_rhs<cpp_daisie_cs_runmod_2>;

    cpp_daisie_cs_runmod_2(param_t&& par) :
      p_(std::move(par))
    {
    }

    // odeint interface
    void operator()(const state_type& x, state_type& dx, double /*t*/) const
    {
      if (++p_.steps > max_cs_steps) {
        throw std::runtime_error("cpp_daisie_cs_runmod_2: too many steps");
      }

      const auto xx1 = padded_vector_view<2>(x.data().begin(), p_.lx);
      const auto xx2 = padded_vector_view<2>(x.data().begin() + p_.lx, p_.lx);
      const auto xx3 = padded_vector_view<2>(x.data().begin() + 2 * p_.lx, p_.lx);

      const auto chunk = p_.lx + 4 + 2 * p_.kk;
      const auto laavec = p_.P.data().begin();
      const auto lacvec = p_.P.data().begin() + chunk;
      const auto muvec = p_.P.data().begin() + 2 * chunk;
      const auto gamvec = p_.P.data().begin() + 3 * chunk;
      const auto nn = p_.P.data().begin() + 4 * chunk;

      // using offsets into our views instead of vectors:
      const int nil2lx = 2;
      const int il1 = nil2lx + p_.kk - 1;
      const int il2 = nil2lx + p_.kk + 1;
      const int il3 = nil2lx + p_.kk;
      const int il4 = nil2lx + p_.kk - 2;

      const int in1 = nil2lx + 2 * p_.kk - 1;
      const int in2 = nil2lx + 1;
      const int in3 = nil2lx + p_.kk;
      const int in4 = nil2lx-1;

      const int ix1 = nil2lx - 1;
      const int ix2 = nil2lx + 1;
      const int ix3 = nil2lx;
      const int ix4 = nil2lx - 2;

      // using views into output vector:
      auto dx1 = dx.data().begin();
      auto dx2 = dx1 + p_.lx;
      auto dx3 = dx2 + p_.lx;

      const auto kk = (1 == p_.kk) ? 1.0 : 0.0;         // make the loop body branch-free
      for (int i = 0; i < p_.lx; ++i) {
        dx1[i] = laavec[il1 + i + 1] * xx2[ix1 + i]
               + lacvec[il4 + i + 1] * xx2[ix4 + i]
               + muvec[il2 + i + 1] * xx2[ix3 + i]
               + lacvec[il1 + i] * nn[in1 + i] * xx1[ix1 + i]
               + muvec[il2 + i] * nn[in2 + i] * xx1[ix2 + i]
               - (muvec[il3 + i] + lacvec[il3 + i]) * nn[in3 + i] * xx1[ix3 + i]
               - gamvec[il3 + i] * xx1[ix3 + i]
               + kk * (laavec[il3 + i] * xx3[ix3 + i] + 2.0 * lacvec[il1 + i] * xx3[ix1 + i]);
        dx2[i] = gamvec[il3 + i] * xx1[ix3 + i]
               + lacvec[il1 + i + 1] * nn[in1 + i] * xx2[ix1 + i]
               + muvec[il2 + i + 1] * nn[in2 + i] * xx2[ix2 + i]
               - (muvec[il3 + 1 + i] + lacvec[il3 + 1 + i]) * nn[in3 + i + 1] * xx2[ix3 + i]
               - laavec[il3 + i] * xx2[ix3 + i];
        dx3[i] = lacvec[il1 + i] * nn[in4 + i] * xx3[ix1 + i]
               + muvec[il2 + i] * nn[in2 + i] * xx3[ix2 + i]
               - (lacvec[il3 + i] + muvec[il3 + i]) * nn[in3 + i] * xx3[ix3 + i]
               - (laavec[il3 + i] + gamvec[il3 + i]) * xx3[ix3 + i];
      }
    }

  private:
    const param_t p_;
  };



  class cpp_daisie_cs_runmod_3
  {
  public:
    using jacobian = const_from_linear_rhs<cpp_daisie_cs_runmod_3>;

    cpp_daisie_cs_runmod_3(param_t&& par) :
      p_(std::move(par))
    {
    }

    // odeint interface
    void operator()(const state_type& x, state_type& dx, double /*t*/) const
    {
      if (++p_.steps > max_cs_steps) {
        throw std::runtime_error("cpp_daisie_cs_runmod_3: too many steps");
      }
      using pvv = padded_vector_view<2>;
      using pmv = padded_mat_view<2>;
      auto lx1 = p_.lx;
      auto lx2 = p_.kk;       // lx1 == lx2 might be violated in the future
      auto p = p_.P.data().begin();
      auto a1 = p; p += lx1;
      auto a2 = p; p += lx1;
      auto a3 = p; p += lx1;
      auto b1 = const_mat_view(p, lx2); p += lx1 * lx2;
      auto b2 = const_mat_view(p, lx2); p += lx1 * lx2;
      auto b3 = const_mat_view(p, lx2); p += lx1 * lx2;
      auto b4 = p; p += lx1;
      auto b5 = p; p += lx1;
      auto b6 = p; p += lx1;
      auto b7 = const_mat_view(p, lx2); p += lx1 * lx2;
      auto b8 = const_mat_view(p, lx2); p += lx1 * lx2;
      auto b9 = const_mat_view(p, lx2); p += lx1 * lx2;
      auto b10 = const_mat_view(p, lx2); p += lx1 * lx2;
      auto b11 = const_mat_view(p, lx2); p += lx1 * lx2;
      auto c1 = const_mat_view(p, lx2); p += lx1 * lx2;
      auto c2 = p; p += lx1;
      auto c3 = const_mat_view(p, lx2); p += lx1 * lx2;
      auto c4 = const_mat_view(p, lx2); p += lx1 * lx2;
      auto c5 = const_mat_view(p, lx2); p += lx1 * lx2;
      auto c6 = const_mat_view(p, lx2); p += lx1 * lx2;
      auto c7 = const_mat_view(p, lx2); 
    
      p = x.data().begin();
      auto x1 = pvv(p, lx1); p += lx1;
      auto x2 = pmv(p, lx1, lx2); p += lx1 * lx2;
      auto x3 = pmv(p, lx1, lx2);
    
      auto dx1 = dx.data().begin();
      auto dx2 = mat_view(dx.data().begin() + lx1, lx2);
      auto dx3 = mat_view(dx.data().begin() + lx1 + lx1 * lx2, lx2);
      for (int i = 0; i < lx1; ++i) {
        dx1[i] = a1[i] * x1[i + 1] +
                 a2[i] * x1[i + 3] -
                 a3[i] * x1[i + 2];  
      }
      for (int i = 0; i < lx1; ++i) {
        for (int j = 0; j < lx2; ++j) {
          dx2(i, j) = b1(i, j) * x3(i + 2, j + 1) +
                      b2(i, j) * x3(i + 2, j + 0) +
                      b3(i, j) * x3(i + 2, j + 2) +
                      b7(i, j) * x2(i + 1, j + 2) +
                      b8(i, j) * x2(i + 2, j + 1) +
                      b9(i, j) * x2(i + 3, j + 2) +
                      b10(i, j) * x2(i + 2, j + 3) -
                      b11(i, j) * x2(i + 2, j + 2);
          dx3(i, j) = c1(i, j) * x2(i + 2, j + 2) +
                      c3(i, j) * x3(i + 1, j + 2) +
                      c4(i, j) * x3(i + 2, j + 1) +
                      c5(i, j) * x3(i + 3, j + 2) +
                      c6(i, j) * x3(i + 2, j + 3) -
                      c7(i, j) * x3(i + 2, j + 2);
              
        }
        dx2(i, 0) += b4[i] * x1[i + 2] +
                     b5[i] * x1[i + 1] +
                     b6[i] * x1[i];
        dx3(i, 0) += c2[i] * x1[i + 2];
      }
    }

  private:
    param_t p_;
  };

} // anonymous namespace


//' Driver for the boost::odeint solver for the CS model
//'
//' @name daisie_odeint_cs
//' @export daisie_odeint_cs
//' @return Object of type `state_type`, which itself is
//' `vector_t`<double>, with the result of the
//' integration depending on the runmod chosen.
RcppExport SEXP daisie_odeint_cs(SEXP rrunmod, SEXP ry, SEXP rtimes, SEXP rlx, SEXP rkk, SEXP rpar, SEXP Stepper, SEXP ratol, SEXP rrtol) {
BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  auto runmod = as<std::string>(rrunmod);
  auto y = as<state_type>(ry);
  auto times = as<std::vector<double>>(rtimes);
  auto lx = as<int>(rlx);
  auto kk = as<int>(rkk);
  auto stepper = as<std::string>(Stepper);
  auto atol = as<double>(ratol);
  auto rtol = as<double>(rrtol);

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
  else if (runmod == "daisie_runmod3") {
    cpp_daisie_cs_runmod_3 rhs(std::move(p));
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

