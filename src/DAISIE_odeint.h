#ifndef DAISIE_ODEINT_H_INCLUDED
#define DAISIE_ODEINT_H_INCLUDED


// boiler-plate code calling into boost::odeint


#define STRICT_R_HEADERS
#include <RcppCommon.h>

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(BH)]]

#include <boost/numeric/ublas/vector.hpp>

// Provide Forward Declarations
namespace Rcpp {

  namespace traits{

    // Setup non-intrusive extension via template specialization for
    // 'ublas' class boost::numeric::ublas

    // Support for wrap
    template <typename T> SEXP wrap(const boost::numeric::ublas::vector<T> & obj);

    // Support for as<T>
    template <typename T> class Exporter< boost::numeric::ublas::vector<T> >;

  }
}


#include <Rcpp.h>
#include <boost/numeric/odeint.hpp>
#include <vector>
#include <algorithm>
#include <stdexcept>


// boost::numeric::ublas wrapping from:
// https://gallery.rcpp.org/articles/custom-templated-wrap-and-as-for-seamingless-interfaces/
namespace Rcpp {

  namespace traits{

    // Defined wrap case
    template <typename T> SEXP wrap(const boost::numeric::ublas::vector<T> & obj){
      const int RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype ;

      return Rcpp::Vector< RTYPE >(obj.begin(), obj.end());
    };


    // Defined as< > case
    template<typename T> class Exporter< boost::numeric::ublas::vector<T> > {
      typedef typename boost::numeric::ublas::vector<T> OUT ;

      // Convert the type to a valid rtype.
      const static int RTYPE = Rcpp::traits::r_sexptype_traits< T >::rtype ;
      Rcpp::Vector<RTYPE> vec;

    public:
      Exporter(SEXP x) : vec(x) {
        if (TYPEOF(x) != RTYPE)
          throw std::invalid_argument("Wrong R type for mapped 1D array");
      }
      OUT get() {
        // Need to figure out a way to perhaps do a pointer pass?
        OUT x(vec.size());
        std::copy(vec.begin(), vec.end(), x.begin()); // have to copy data
        return x;
      }
    };

  }

}


using namespace Rcpp;
using namespace boost::numeric::odeint;


// type of the ode state
using state_type = boost::numeric::ublas::vector<double>;



namespace daisie_odeint {


extern double abm_factor;


template <typename Stepper, typename Rhs>
inline void do_integrate(double atol, double rtol, Rhs rhs, state_type& y, double t0, double t1)
{
  integrate_adaptive(make_controlled<Stepper>(atol, rtol), rhs, y, t0, t1, 0.1 * (t1 - t0));
}


template <size_t Steps, typename Rhs>
inline void abm(Rhs rhs, state_type& y, double t0, double t1)
{
  auto abm = adams_bashforth_moulton<Steps, state_type>{};
  abm.initialize(rhs, y, t0, abm_factor * (t1 - t0));
  integrate_const(abm, rhs, y, t0, t1, abm_factor * (t1 - t0));
}


template <size_t Steps, typename Rhs>
inline void ab(Rhs rhs, state_type& y, double t0, double t1)
{
  auto ab = adams_bashforth<Steps, state_type>{};
  ab.initialize(rhs, y, t0, abm_factor * (t1 - t0));
  integrate_const(ab, rhs, y, t0, t1, abm_factor * (t1 - t0));
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
    else if (0 == stepper.compare(0, stepper.size() - 2, "odeint::adams_bashforth")) {
      const char steps = stepper.back();
      switch (steps) {
      case '1': ab<1>(rhs, y, t0, t1); break;
      case '2': ab<2>(rhs, y, t0, t1); break;
      case '3': ab<3>(rhs, y, t0, t1); break;
      case '4': ab<4>(rhs, y, t0, t1); break;
      case '5': ab<5>(rhs, y, t0, t1); break;
      case '6': ab<6>(rhs, y, t0, t1); break;
      case '7': ab<7>(rhs, y, t0, t1); break;
      case '8': ab<8>(rhs, y, t0, t1); break;
      default: throw std::runtime_error("DAISIE_odeint_helper::integrate: unsupported steps for admam_bashforth_moulton");
      }
    }
    else if (0 == stepper.compare(0, stepper.size() - 2, "odeint::adams_bashforth_moulton")) {
      const char steps = stepper.back();
      switch (steps) {
      case '1': abm<1>(rhs, y, t0, t1); break;
      case '2': abm<2>(rhs, y, t0, t1); break;
      case '3': abm<3>(rhs, y, t0, t1); break;
      case '4': abm<4>(rhs, y, t0, t1); break;
      case '5': abm<5>(rhs, y, t0, t1); break;
      case '6': abm<6>(rhs, y, t0, t1); break;
      case '7': abm<7>(rhs, y, t0, t1); break;
      case '8': abm<8>(rhs, y, t0, t1); break;
      default: throw std::runtime_error("DAISIE_odeint_helper::integrate: unsupported steps for admam_bashforth_moulton");
      }
    }
    else {
      throw std::runtime_error("DAISIE_odeint_helper::integrate: unknown stepper");
    }
  }

}

#endif
