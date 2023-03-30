//
//  Copyright (c) 2023, Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef DAISIE_TYPES_H_INCLUDED
#define DAISIE_TYPES_H_INCLUDED

#include "config.h"
#include <RcppCommon.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>


namespace ublas = boost::numeric::ublas;


template <typename T> using vector_t = ublas::vector<T>;
template <typename T> using matrix_t = ublas::matrix<T>;


// forward declarations Rcpp <-> boost::numeric::ublas
namespace Rcpp {

  namespace traits {

    template <typename T> SEXP wrap(const vector_t<T>&);
    template <typename T> SEXP wrap(const matrix_t<T>&);

    template <typename T> vector_t<T> as(SEXP);
    template <typename T> matrix_t<T> as(SEXP);

    template <typename T> class Exporter<vector_t<T>>;
    template <typename T> class Exporter<matrix_t<T>>;

  }

}


#include <Rcpp.h>


namespace Rcpp {

  namespace traits {

    template <typename T> inline SEXP wrap(const vector_t<T>& obj) {
      const int RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype;
      return Rcpp::Vector<RTYPE>(obj.begin(), obj.end());
    }


    template <typename T> inline SEXP wrap(const matrix_t<T>& obj) {
      const size_t nr = static_cast<size_t>(obj.size1());
      const size_t nc = static_cast<size_t>(obj.size2());
      const int RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype;
      Rcpp::Matrix<RTYPE> rmat(nr, nc);
      for (size_t i = 0; i < nr; ++i) {
        for (size_t j = 0; j < nc; ++j) {
          rmat(i, j) = obj(i, j);
        }
      }
      return rmat;
    }


    template<typename T>
    class Exporter<vector_t<T>>
    {
    private:
      static constexpr int RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype;
      Rcpp::Vector<RTYPE> rvec;

    public:
      Exporter(SEXP x) : rvec(x) {
        if (TYPEOF(x) != RTYPE) {
          throw std::invalid_argument("Wrong R type for mapped 1D array");
        }
      }

      vector_t<T> get() {
        vector_t<T> x(rvec.size());
        std::copy(rvec.begin(), rvec.end(), x.begin());
        return x;
      }
    };


    template<typename T>
    class Exporter<matrix_t<T>>
    {
    private:
      static constexpr int RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype;
      Rcpp::Matrix<RTYPE> rmat;

    public:
      Exporter(SEXP x) : rmat(x) {
        if (TYPEOF(x) != RTYPE) {
          throw std::invalid_argument("Wrong R type for mapped 2D array");
        }
      }

      matrix_t<T> get() {
        const size_t nr = static_cast<size_t>(rmat.rows());
        const size_t nc = static_cast<size_t>(rmat.cols());
        matrix_t<T> x(nr, nc);
        for (size_t i = 0; i < nr; ++i) {
          for (size_t j = 0; j < nc; ++j) {
            x(i, j) = rmat(i, j);
          }
        }
        return x;
      }
    };

  }

}

#endif
