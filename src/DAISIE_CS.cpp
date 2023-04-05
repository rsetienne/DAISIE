//
//  Copyright (c) 2023, Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

#include "config.h"
#include "DAISIE_odeint.h"


//' Driver for the boost::odeint solver
//'
//' @export
// [[Rcpp::export]]
state_type DAISIE_odeint_cs(std::string runmod, 
                            state_type y, 
                            std::vector<double> times, 
                            int lx, 
                            int kk, 
                            state_type par,
                            std::string Stepper, 
                            double atol, 
                            double rtol);

//' CS iteration control
//'
//' Sets or retrieves the max. number of iterations used by the odeint solver.
//'
//' @param max_steps \code{max_steps}: sets max. iterations to \code{max_steps}. \cr
//' @return current max. iterations
//'
//' @export
// [[Rcpp::export]]
int DAISIE_CS_max_steps(int rmax_steps);


//` adams_bashforth and adams_bashforth_moulton integration control
//'
//' Sets or retrieves the factor to calculate the step-size used by the odeint::adams_bashforth[_moulton] solvers.
//'
//' @param factor sets step-size to \code{factor * (t1 - t0)}. \cr
//' @return current factor
//'
//' @export
// [[Rcpp::export]]
double DAISIE_abm_factor(double rfactor);


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
      constexpr int nil2lx = 2;
      const int il1 = nil2lx + p_.kk - 1;
      const int il2 = nil2lx + p_.kk + 1;
      const int il3 = nil2lx + p_.kk;
      const int il4 = nil2lx + p_.kk - 2;

      const int in1 = nil2lx + 2 * p_.kk - 1;
      const int in2 = nil2lx + 1;
      const int in3 = nil2lx + p_.kk;

      constexpr int ix1 = nil2lx - 1;
      constexpr int ix2 = nil2lx + 1;
      constexpr int ix3 = nil2lx;
      constexpr int ix4 = nil2lx - 2;

      auto dx1 = dx.data().begin();
      auto dx2 = dx1 + p_.lx;
      auto dx3 = dx2 + p_.lx;

      for (int i = 0; i < p_.lx; ++i) {
        // common subexpression elimination left to the compiler
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
      constexpr int nil2lx = 2;
      const int il1 = nil2lx + p_.kk - 1;
      const int il2 = nil2lx + p_.kk + 1;
      const int il3 = nil2lx + p_.kk;
      const int il4 = nil2lx + p_.kk - 2;

      const int in1 = nil2lx + 2 * p_.kk - 1;
      const int in2 = nil2lx + 1;
      const int in3 = nil2lx + p_.kk;

      constexpr int ix1 = nil2lx - 1;
      constexpr int ix2 = nil2lx + 1;
      constexpr int ix3 = nil2lx;
      constexpr int ix4 = nil2lx - 2;

      // using views into output vector:
      auto dx1 = dx.data().begin();
      auto dx2 = dx1 + p_.lx;
      auto dx3 = dx2 + p_.lx;
      auto dx4 = dx3 + p_.lx;

      for (int i = 0; i < p_.lx; ++i) {
        // common subexpression elimination left to the compiler
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
      constexpr int nil2lx = 2;
      const int il1 = nil2lx + p_.kk - 1;
      const int il2 = nil2lx + p_.kk + 1;
      const int il3 = nil2lx + p_.kk;
      const int il4 = nil2lx + p_.kk - 2;

      const int in1 = nil2lx + 2 * p_.kk - 1;
      const int in2 = nil2lx + 1;
      const int in3 = nil2lx + p_.kk;
      const int in4 = nil2lx-1;

      constexpr int ix1 = nil2lx - 1;
      constexpr int ix2 = nil2lx + 1;
      constexpr int ix3 = nil2lx;
      constexpr int ix4 = nil2lx - 2;

      // using views into output vector:
      auto dx1 = dx.data().begin();
      auto dx2 = dx1 + p_.lx;
      auto dx3 = dx2 + p_.lx;

      const auto kk = (1 == p_.kk) ? 1.0 : 0.0;         // make the loop body branch-free
      for (int i = 0; i < p_.lx; ++i) {
        // common subexpression elimination left to the compiler
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

} // anonymous namespace


state_type DAISIE_odeint_cs(std::string runmod, 
                            state_type y, 
                            std::vector<double> times, 
                            int lx, 
                            int kk, 
                            state_type par,
                            std::string Stepper, 
                            double atol, 
                            double rtol) 
{
  auto p = param_t(lx, kk, std::move(par));
  if (runmod == "daisie_runmod") {
    cpp_daisie_cs_runmod rhs(std::move(p));
    daisie_odeint::integrate(Stepper, std::ref(rhs), y, times[0], times[1], atol, rtol);
  }
  else if (runmod == "daisie_runmod1") {
    cpp_daisie_cs_runmod_1 rhs(std::move(p));
    daisie_odeint::integrate(Stepper, std::ref(rhs), y, times[0], times[1], atol, rtol);
  }
  else if (runmod == "daisie_runmod2") {
    cpp_daisie_cs_runmod_2 rhs(std::move(p));
    daisie_odeint::integrate(Stepper, std::ref(rhs), y, times[0], times[1], atol, rtol);
  }
  else {
    throw std::runtime_error("DAISIE_odeint_cs: unknown runmod");
  }
  return y;
}


int DAISIE_CS_max_steps(int rmax_steps) 
{
  max_cs_steps = (0 < rmax_steps) ? rmax_steps : default_max_cs_steps;
  return max_cs_steps;
}


namespace daisie_odeint {

  // step-size factor for adams_bashforth_moulton integration
  constexpr double default_abm_factor = 0.0001;
  double abm_factor = default_abm_factor;

}


double DAISIE_abm_factor(double rfactor) 
{
  daisie_odeint::abm_factor = (0 < rfactor) ? rfactor : daisie_odeint::default_abm_factor;
  return daisie_odeint::abm_factor;
}
