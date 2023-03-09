#ifndef ODEINT_CONFIG_H_INCLUDED
#define ODEINT_CONFIG_H_INCLUDED

// [[Rcpp::plugins(cpp14)]]

// Special case to make use of some steppers that would include
// boost/functional.hpp
// moved to Makevars[.win]
/*
#if __cplusplus >= 201703L
# ifdef _HAS_AUTO_PTR_ETC
#  undef _HAS_AUTO_PTR_ETC
# endif
# define _HAS_AUTO_PTR_ETC 0
#endif
 */

// Special case to make use of some steppers that would include
// boost/get_pointer.hpp
#ifndef BOOST_NO_AUTO_PTR
# define BOOST_NO_AUTO_PTR
#endif

// uncomment if unitialized member variable bulirsch_stoer::m_dt_last
// is fixed in boost (BH)
#define USE_BULRISCH_STOER_PATCH

#endif
