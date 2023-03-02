#ifndef ODEINT_CONFIG_H_INCLUDED
#define ODEINT_CONFIG_H_INCLUDED

// [[Rcpp::plugins(cpp14)]]

// Special case to make use of some steppers that would include
// boost/functional.hpp
#if __cplusplus >= 201703L
#ifdef _HAS_AUTO_PTR_ETC
#undef _HAS_AUTO_PTR_ETC
#endif
#define _HAS_AUTO_PTR_ETC 0
#endif

// Special case to make use of some steppers that would include
// boost/get_pointer.hpp
#ifdef BOOST_NO_AUTO_PTR
#undef BOOST_NO_AUTO_PTR
#endif
#define BOOST_NO_AUTO_PTR

#endif
