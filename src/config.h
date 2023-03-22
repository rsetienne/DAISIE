//
//  Copyright (c) 2023, Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef ODEINT_CONFIG_H_INCLUDED
#define ODEINT_CONFIG_H_INCLUDED

// [[Rcpp::plugins(cpp14)]]

// Special case to make use of some steppers that would include
// boost/get_pointer.hpp
#ifndef BOOST_NO_AUTO_PTR
# define BOOST_NO_AUTO_PTR
#endif

// Addresses unitialized member variable bulirsch_stoer<>::m_dt_last.
//
// The issue is *not* fixed in BOOST_VERSION 1.81.1.
// We need to check for fixes in upcomming boost (BH) releases.
//
// Uncomment if unitialized member variable bulirsch_stoer::m_dt_last
// is fixed in boost (BH):
#define USE_BULRISCH_STOER_PATCH

#endif
