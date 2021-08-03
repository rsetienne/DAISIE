#ifndef DAISIE_ODEINT_H_INCLUDED
#define DAISIE_ODEINT_H_INCLUDED


// boiler-plate code calling into boost::odeint


#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <boost/numeric/odeint.hpp>
#include <vector>
#include <stdexcept>


using namespace Rcpp;
using namespace boost::numeric::odeint;


// type of the ode state
using state_type = std::vector<double>;


namespace daisie_odeint {


  template <typename Stepper, typename Rhs>
  inline void do_integrate(double atol, double rtol, Rhs rhs, state_type& y, double t0, double t1)
  {
    integrate_adaptive(make_controlled<Stepper>(atol, rtol), rhs, y, t0, t1, 0.1 * (t1 - t0));
  }


  // wrapper around odeint::integrate
  // maps runtime stepper name -> compiletime odeint::stepper type
  template <typename Rhs>
  inline void integrate(
      const std::string& stepper,
      Rhs rhs,
      state_type& y,
      double t0,
      double t1,
      double atol,
      double rtol)
  {
    if ("odeint::runge_kutta_cash_karp54" == stepper) {
      do_integrate<runge_kutta_cash_karp54<state_type>>(atol, rtol, rhs, y, t0, t1);
    }
    else if ("odeint::runge_kutta_fehlberg78" == stepper) {
      do_integrate<runge_kutta_fehlberg78<state_type>>(atol, rtol, rhs, y, t0, t1);
    }
    else if ("odeint::runge_kutta_dopri5" == stepper) {
      do_integrate<runge_kutta_dopri5<state_type>>(atol, rtol, rhs, y, t0, t1);
    }
    else if ("odeint::bulirsch_stoer" == stepper) {
      // outlier in calling convention
      using stepper_t = bulirsch_stoer<state_type>;
      integrate_adaptive(stepper_t(atol, rtol), rhs, y, t0, t1, 0.1 * (t1 - t0));
    }
    else {
      throw std::runtime_error("DAISIE_odeint_helper::integrate: unknown stepper");
    }
  }

}

#endif
