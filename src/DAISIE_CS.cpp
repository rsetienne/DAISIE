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
#include "DAISIE_time_rates.h"
#include "DAISIE_stiff.h"

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

    static constexpr bool stiff_linear = true;
    static constexpr bool stiff_time_dependent = false;

    // odeint interface
    void operator()(const state_type& x, state_type& dx, double t) const
    {
      if (++p_.steps > max_cs_steps) {
        throw std::runtime_error("cpp_daisie_cs_runmod: too many steps");
      }
      eval(x, dx, t);
    }

    // right-hand side without step counting (used for Jacobian probing)
    void eval(const state_type& x, state_type& dx, double /*t*/) const
    {

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

    static constexpr bool stiff_linear = true;
    static constexpr bool stiff_time_dependent = false;

    // odeint interface
    void operator()(const state_type& x, state_type& dx, double t) const
    {
      if (++p_.steps > max_cs_steps) {
        throw std::runtime_error("cpp_daisie_cs_runmod_1: too many steps");
      }
      eval(x, dx, t);
    }

    // right-hand side without step counting (used for Jacobian probing)
    void eval(const state_type& x, state_type& dx, double /*t*/) const
    {

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

    static constexpr bool stiff_linear = true;
    static constexpr bool stiff_time_dependent = false;

    // odeint interface
    void operator()(const state_type& x, state_type& dx, double t) const
    {
      if (++p_.steps > max_cs_steps) {
        throw std::runtime_error("cpp_daisie_cs_runmod_2: too many steps");
      }
      eval(x, dx, t);
    }

    // right-hand side without step counting (used for Jacobian probing)
    void eval(const state_type& x, state_type& dx, double /*t*/) const
    {

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



  // Right-hand side of the rhs3 system (clade with a second, later colonist),
  // shared by the time-independent and the time-dependent class. P holds the
  // coefficient arrays a1 ... c7 in the layout of DAISIE_loglik_rhs_precomp2().
  inline void rhs3_body(int lx1, int lx2, const double* P, const state_type& x, state_type& dx)
  {
      using pvv = padded_vector_view<2>;
      using pmv = padded_mat_view<2>;
      auto p = P;
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

  class cpp_daisie_cs_runmod_3
  {
  public:
    using jacobian = const_from_linear_rhs<cpp_daisie_cs_runmod_3>;

    cpp_daisie_cs_runmod_3(param_t&& par) :
      p_(std::move(par))
    {
    }

    static constexpr bool stiff_linear = true;
    static constexpr bool stiff_time_dependent = false;

    // odeint interface
    void operator()(const state_type& x, state_type& dx, double t) const
    {
      if (++p_.steps > max_cs_steps) {
        throw std::runtime_error("cpp_daisie_cs_runmod_3: too many steps");
      }
      eval(x, dx, t);
    }

    // right-hand side without step counting (used for Jacobian probing)
    void eval(const state_type& x, state_type& dx, double /*t*/) const
    {
      rhs3_body(p_.lx, p_.kk, p_.P.data().begin(), x, dx);   // lx1 == lx2 might be violated in the future
    }

  private:
    param_t p_;
  };


  // Jacobian policy for the time-dependent linear system dx/dt = A(t) x:
  // J(t) is evaluated column-wise from unit vectors at the requested time and
  // dfdt by a central finite difference in t.
  template <typename RHS>
  struct time_jacobian
  {
    explicit time_jacobian(RHS& rhs) : rhs_(rhs) {}

    void operator()(const vector_t<double>& x, matrix_t<double>& J, double t, vector_t<double>& dfdt)
    {
      const size_t n = x.size();
      auto single = vector_t<double>(n, 0.0);
      auto dxdt = vector_t<double>(n, 0.0);
      rhs_.set_time(t);   // rate vectors are updated once per Jacobian
      for (size_t i = 0; i < n; ++i) {
        single[i] = 1.0;
        rhs_.apply(single, dxdt);
        for (size_t r = 0; r < n; ++r) J(r, i) = dxdt[r];
        single[i] = 0.0;
      }
      const double h = 1e-6 * std::max(1.0, std::fabs(t));
      auto fp = vector_t<double>(n, 0.0);
      auto fm = vector_t<double>(n, 0.0);
      rhs_.eval(x, fp, t + h);
      rhs_.eval(x, fm, t - h);
      rhs_.set_time(t);
      for (size_t r = 0; r < n; ++r) dfdt[r] = (fp[r] - fm[r]) / (2.0 * h);
    }

    RHS& rhs_;
  };


  // ---------------------------------------------------------------------------
  // Time-dependent (island ontogeny / sea level) CS model.
  //
  // C++ analogs of DAISIE_loglik_rhs_time(), _time1() and _time2() in
  // R/DAISIE_loglik_CS_time.R. The state-vector layout and the equations are
  // those of the time-independent classes above; the difference is that the
  // per-capita rate vectors (laavec, lacvec, muvec, gamvec) depend on the
  // island area A(t) and are recomputed in every call of the right-hand side.
  // ---------------------------------------------------------------------------

  struct time_param_t
  {
    time_param_t(int LX, int KK, const state_type& par) :
      lx(LX), kk(KK), pars(par.data().begin()), rates(LX, KK)
    {
      if (par.size() != static_cast<size_t>(daisie_time::pars_t::size)) {
        throw std::runtime_error("daisie_odeint_cs_time: expected 18 model parameters");
      }
    }

    int lx;
    int kk;
    daisie_time::pars_t pars;
    daisie_time::rate_vectors rates;
    mutable int steps = 0;
  };


  // Shared plumbing: step counter, rate update, offsets into the padded views.
  // The numerical Jacobian of the (linear in x, but time-dependent) system is
  // provided by 'time_jacobian' below, which needs the non-counting eval().
  template <typename Derived>
  class cpp_daisie_cs_time_base
  {
  public:
    using jacobian = time_jacobian<Derived>;
    static constexpr bool stiff_linear = true;
    static constexpr bool stiff_time_dependent = true;

    explicit cpp_daisie_cs_time_base(time_param_t&& par) : p_(std::move(par)) {}

    void operator()(const state_type& x, state_type& dx, double t) const
    {
      if (++p_.steps > max_cs_steps) {
        throw std::runtime_error("cpp_daisie_cs_time: too many steps");
      }
      eval(x, dx, t);
    }

    // The following three members are used by 'time_jacobian' and do not
    // count steps: set_time() updates the area-dependent rate vectors,
    // apply() evaluates the right-hand side with the current rates, and
    // eval() does both.
    void set_time(double t) const { p_.rates.update(t, p_.pars); }

    void apply(const state_type& x, state_type& dx) const
    {
      static_cast<const Derived*>(this)->rhs(x, dx);
    }

    void eval(const state_type& x, state_type& dx, double t) const
    {
      set_time(t);
      apply(x, dx);
    }

  protected:
    mutable time_param_t p_;
  };


  class cpp_daisie_cs_time_runmod : public cpp_daisie_cs_time_base<cpp_daisie_cs_time_runmod>
  {
  public:
    using cpp_daisie_cs_time_base<cpp_daisie_cs_time_runmod>::cpp_daisie_cs_time_base;

    void rhs(const state_type& x, state_type& dx) const
    {
      const int lx = p_.lx;
      const int kk_ = p_.kk;
      const auto& R = p_.rates;
      const double* laavec = R.laavec.data();
      const double* lacvec = R.lacvec.data();
      const double* muvec = R.muvec.data();
      const double* gamvec = R.gamvec.data();
      const double* nn = R.nn.data();
      const int nil2lx = 2;
      const int il1 = nil2lx + kk_ - 1;
      const int il2 = nil2lx + kk_ + 1;
      const int il3 = nil2lx + kk_;
      const int il4 = nil2lx + kk_ - 2;
      const int in1 = nil2lx + 2 * kk_ - 1;
      const int in2 = nil2lx + 1;
      const int in3 = nil2lx + kk_;
      const int ix1 = nil2lx - 1;
      const int ix2 = nil2lx + 1;
      const int ix3 = nil2lx;
      const int ix4 = nil2lx - 2;
      const auto xx1 = padded_vector_view<2>(x.data().begin(), lx);
      const auto xx2 = padded_vector_view<2>(x.data().begin() + lx, lx);
      auto dx1 = dx.data().begin();
      auto dx2 = dx1 + lx;
      auto dx3 = dx2 + lx;
      for (int i = 0; i < lx; ++i) {
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
               - laavec[il3 + i + 1] * xx2[ix3 + i];
      }
      dx3[0] = 0.0;
    }
  };

  class cpp_daisie_cs_time_runmod_1 : public cpp_daisie_cs_time_base<cpp_daisie_cs_time_runmod_1>
  {
  public:
    using cpp_daisie_cs_time_base<cpp_daisie_cs_time_runmod_1>::cpp_daisie_cs_time_base;

    void rhs(const state_type& x, state_type& dx) const
    {
      const int lx = p_.lx;
      const int kk_ = p_.kk;
      const auto& R = p_.rates;
      const double* laavec = R.laavec.data();
      const double* lacvec = R.lacvec.data();
      const double* muvec = R.muvec.data();
      const double* gamvec = R.gamvec.data();
      const double* nn = R.nn.data();
      const int nil2lx = 2;
      const int il1 = nil2lx + kk_ - 1;
      const int il2 = nil2lx + kk_ + 1;
      const int il3 = nil2lx + kk_;
      const int il4 = nil2lx + kk_ - 2;
      const int in1 = nil2lx + 2 * kk_ - 1;
      const int in2 = nil2lx + 1;
      const int in3 = nil2lx + kk_;
      const int ix1 = nil2lx - 1;
      const int ix2 = nil2lx + 1;
      const int ix3 = nil2lx;
      const int ix4 = nil2lx - 2;
      const auto xx1 = padded_vector_view<2>(x.data().begin(), lx);
      const auto xx2 = padded_vector_view<2>(x.data().begin() + lx, lx);
      const auto xx3 = padded_vector_view<2>(x.data().begin() + 2 * lx, lx);
      const auto xx4 = padded_vector_view<2>(x.data().begin() + 3 * lx, lx);
      auto dx1 = dx.data().begin();
      auto dx2 = dx1 + lx;
      auto dx3 = dx2 + lx;
      auto dx4 = dx3 + lx;
      for (int i = 0; i < lx; ++i) {
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
               - laavec[il3 + i + 1] * xx2[ix3 + i];
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
  };

  class cpp_daisie_cs_time_runmod_2 : public cpp_daisie_cs_time_base<cpp_daisie_cs_time_runmod_2>
  {
  public:
    using cpp_daisie_cs_time_base<cpp_daisie_cs_time_runmod_2>::cpp_daisie_cs_time_base;

    void rhs(const state_type& x, state_type& dx) const
    {
      const int lx = p_.lx;
      const int kk_ = p_.kk;
      const auto& R = p_.rates;
      const double* laavec = R.laavec.data();
      const double* lacvec = R.lacvec.data();
      const double* muvec = R.muvec.data();
      const double* gamvec = R.gamvec.data();
      const double* nn = R.nn.data();
      const int nil2lx = 2;
      const int il1 = nil2lx + kk_ - 1;
      const int il2 = nil2lx + kk_ + 1;
      const int il3 = nil2lx + kk_;
      const int il4 = nil2lx + kk_ - 2;
      const int in1 = nil2lx + 2 * kk_ - 1;
      const int in2 = nil2lx + 1;
      const int in3 = nil2lx + kk_;
      const int ix1 = nil2lx - 1;
      const int ix2 = nil2lx + 1;
      const int ix3 = nil2lx;
      const int ix4 = nil2lx - 2;
      const int in4 = nil2lx - 1;
      const auto xx1 = padded_vector_view<2>(x.data().begin(), lx);
      const auto xx2 = padded_vector_view<2>(x.data().begin() + lx, lx);
      const auto xx3 = padded_vector_view<2>(x.data().begin() + 2 * lx, lx);
      auto dx1 = dx.data().begin();
      auto dx2 = dx1 + lx;
      auto dx3 = dx2 + lx;
      const auto k1 = (1 == kk_) ? 1.0 : 0.0;   // branch-free
      for (int i = 0; i < lx; ++i) {
        dx1[i] = laavec[il1 + i + 1] * xx2[ix1 + i]
               + lacvec[il4 + i + 1] * xx2[ix4 + i]
               + muvec[il2 + i + 1] * xx2[ix3 + i]
               + lacvec[il1 + i] * nn[in1 + i] * xx1[ix1 + i]
               + muvec[il2 + i] * nn[in2 + i] * xx1[ix2 + i]
               - (muvec[il3 + i] + lacvec[il3 + i]) * nn[in3 + i] * xx1[ix3 + i]
               - gamvec[il3 + i] * xx1[ix3 + i]
               + k1 * (laavec[il3 + i] * xx3[ix3 + i] + 2.0 * lacvec[il1 + i] * xx3[ix1 + i]);
        dx2[i] = gamvec[il3 + i] * xx1[ix3 + i]
               + lacvec[il1 + i + 1] * nn[in1 + i] * xx2[ix1 + i]
               + muvec[il2 + i + 1] * nn[in2 + i] * xx2[ix2 + i]
               - (muvec[il3 + 1 + i] + lacvec[il3 + 1 + i]) * nn[in3 + i + 1] * xx2[ix3 + i]
               - laavec[il3 + i + 1] * xx2[ix3 + i];
        dx3[i] = lacvec[il1 + i] * nn[in4 + i] * xx3[ix1 + i]
               + muvec[il2 + i] * nn[in2 + i] * xx3[ix2 + i]
               - (lacvec[il3 + i] + muvec[il3 + i]) * nn[in3 + i] * xx3[ix3 + i]
               - (laavec[il3 + i] + gamvec[il3 + i]) * xx3[ix3 + i];
      }
    }
  };


  // Time-dependent analog of DAISIE_loglik_rhs3 (DAISIE_loglik_rhs_time3() in
  // R/DAISIE_loglik_CS_time.R). The coefficients a1 ... c7 of
  // DAISIE_loglik_rhs_precomp2() are rebuilt at every time t from the
  // per-capita rates at the area A(t); the number of species in a cell
  // (n1, n2) is N = n1 + n2 (k = 0). Here lx1 == lx2 == lx and kk == lx
  // (so that the rate vectors cover N + 3).
  class cpp_daisie_cs_time_runmod_3 : public cpp_daisie_cs_time_base<cpp_daisie_cs_time_runmod_3>
  {
  public:
    explicit cpp_daisie_cs_time_runmod_3(time_param_t&& par) :
      cpp_daisie_cs_time_base<cpp_daisie_cs_time_runmod_3>(std::move(par))
    {
      const size_t lx = static_cast<size_t>(p_.lx);
      P_.assign(7 * lx + 14 * lx * lx, 0.0);
    }

    void rhs(const state_type& x, state_type& dx) const
    {
      const int lx = p_.lx;
      const auto& R = p_.rates;
      const double* lacv = R.lacvec.data();   // lacv[m + 2]: per-capita rate with m species
      const double* gamv = R.gamvec.data();
      const double mu = R.muvec[0];
      const double laa = R.laavec[0];
      double* p = P_.data();
      double* a1 = p; p += lx;
      double* a2 = p; p += lx;
      double* a3 = p; p += lx;
      double* b1 = p; p += lx * lx;
      double* b2 = p; p += lx * lx;
      double* b3 = p; p += lx * lx;
      double* b4 = p; p += lx;
      double* b5 = p; p += lx;
      double* b6 = p; p += lx;
      double* b7 = p; p += lx * lx;
      double* b8 = p; p += lx * lx;
      double* b9 = p; p += lx * lx;
      double* b10 = p; p += lx * lx;
      double* b11 = p; p += lx * lx;
      double* c1 = p; p += lx * lx;
      double* c2 = p; p += lx;
      double* c3 = p; p += lx * lx;
      double* c4 = p; p += lx * lx;
      double* c5 = p; p += lx * lx;
      double* c6 = p; p += lx * lx;
      double* c7 = p;
      auto lac = [&](int m) { return lacv[std::max(m, 0) + 2]; };
      auto gam = [&](int m) { return gamv[std::max(m, 0) + 2]; };
      for (int i = 0; i < lx; ++i) {
        const double im1 = std::max(0, i - 1);
        a1[i] = lac(i) * im1;
        a2[i] = mu * (i + 1);
        a3[i] = (lac(i + 1) + mu) * i + gam(i + 1) + lac(i + 1) + laa + mu;
        b4[i] = mu;
        b5[i] = laa;
        b6[i] = lac(i - 1);
        c2[i] = gam(i + 1);
      }
      // matrices are column major: element (i, j) at i + lx * j
      for (int j = 0; j < lx; ++j) {
        const double jm1 = std::max(0, j - 1);
        for (int i = 0; i < lx; ++i) {
          const int k = i + lx * j;
          const int N = i + j;
          const double im1 = std::max(0, i - 1);
          const double L0 = lac(N), Lm = lac(N - 1), Lp = lac(N + 1);
          const double G0 = gam(N);
          b1[k] = laa;
          b2[k] = Lm;
          b3[k] = mu;
          b7[k] = Lm * im1;
          b8[k] = Lm * jm1;
          b9[k] = mu * (i + 1);
          b10[k] = mu * (j + 1);
          b11[k] = (L0 + mu) * N + G0;
          c1[k] = G0;
          c3[k] = L0 * im1;
          c4[k] = L0 * jm1;
          c5[k] = mu * (i + 1);
          c6[k] = mu * (j + 1);
          c7[k] = (Lp + mu) * N + laa + Lp + mu;
        }
      }
      rhs3_body(lx, lx, P_.data(), x, dx);
    }

  private:
    mutable std::vector<double> P_;
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


//' Driver for the boost::odeint solver for the time-dependent CS model
//'
//' C++ analog of the time-dependent integration done with `deSolve::ode()`
//' and `DAISIE_loglik_rhs_time*()` in `DAISIE_integrate_time()`.
//' Arguments as in `daisie_odeint_cs()`, except that `rpar` holds the 18
//' time-dependent model parameters (see `DAISIE_loglik_rhs_time()`).
//'
//' @name daisie_odeint_cs_time
//' @export daisie_odeint_cs_time
//' @return Object of type `state_type`, with the result of the integration
//' depending on the runmod chosen.
RcppExport SEXP daisie_odeint_cs_time(SEXP rrunmod, SEXP ry, SEXP rtimes, SEXP rlx, SEXP rkk, SEXP rpar, SEXP Stepper, SEXP ratol, SEXP rrtol) {
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

  auto p = time_param_t(lx, kk, as<state_type>(rpar));
  if (runmod == "daisie_runmod") {
    cpp_daisie_cs_time_runmod rhs(std::move(p));
    daisie_odeint::integrate(stepper, std::ref(rhs), y, times[0], times[1], atol, rtol);
  }
  else if (runmod == "daisie_runmod1") {
    cpp_daisie_cs_time_runmod_1 rhs(std::move(p));
    daisie_odeint::integrate(stepper, std::ref(rhs), y, times[0], times[1], atol, rtol);
  }
  else if (runmod == "daisie_runmod2") {
    cpp_daisie_cs_time_runmod_2 rhs(std::move(p));
    daisie_odeint::integrate(stepper, std::ref(rhs), y, times[0], times[1], atol, rtol);
  }
  else if (runmod == "daisie_runmod3") {
    cpp_daisie_cs_time_runmod_3 rhs(std::move(p));
    daisie_odeint::integrate(stepper, std::ref(rhs), y, times[0], times[1], atol, rtol);
  }
  else {
    throw std::runtime_error("daisie_odeint_cs_time: unknown runmod");
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

