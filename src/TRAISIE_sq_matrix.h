//
//  Copyright (c) 2025, Thijs Janzen and Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <Rcpp.h>
#include <vector>

template <typename T>
class vector_view_t {
public:
  vector_view_t(T* data, size_t n) : first_(data), n_(n) {}

  size_t size() const noexcept { return n_; }
  T* begin() const noexcept { return first_; }
  T* end() const noexcept { return first_ + n_; }
  T& operator[](size_t i) const { return *(first_ + i); }
  void advance(size_t s) noexcept { first_ += s; }

private:
  T* first_ = nullptr;
  size_t n_ = 0;
};

struct sq_matrix {
  std::vector<double> data_;
  const size_t n_;

  explicit sq_matrix(const Rcpp::NumericMatrix& mat) : n_(mat.nrow()) {
    data_ = std::vector<double>(n_ * n_);
    for (size_t i = 0; i < n_; ++i) {
      for (size_t j = 0; j < n_; ++j) {
        data_[i * n_ + j] = mat(i, j);
      }
    }
  }

  sq_matrix(const Rcpp::NumericVector& vec,
            size_t n) : n_(n) {
    data_ = std::vector<double>(n_ * n_);
    for (size_t i = 0; i < n_; ++i) {
      for (size_t j = 0; j < n_; ++j) {
        data_[i * n_ + j] = vec(j);
      }
    }
  }

  sq_matrix(const sq_matrix& other) : n_(other.n_) {
    data_ = other.data_;
  }

  explicit sq_matrix(size_t n) : n_(n) {
    data_ = std::vector<double>(n_ * n_);
  }

  constexpr size_t index( size_t i, size_t j) const noexcept {
    return i * n_ + j;
  }

  constexpr double value(size_t i, size_t j) const noexcept {
    return data_[index(i,j )];
  }

  size_t size() const{
    return n_;
  }

  std::vector<double> operator*(const vector_view_t<const double>& v) const {
    std::vector<double> out(n_);

    for (size_t i = 0; i < n_; ++i) {
      double s = 0.0;
      for (size_t j = 0; j < n_; ++j) {
        s += v[j] * value(i, j);
      }
      out[i] = s;
    }

    return out;
  }


  std::vector<double> row_sums() const {
    std::vector<double> out(n_, 0.0);
    for (size_t i = 0; i < n_; ++i) {
      for (size_t j = 0; j < n_; ++j) {
        out[i] += value(i, j);
      }
    }
    return out;
  }

  std::vector<double> non_self() const {
    std::vector<double> out = row_sums();
    for (size_t i = 0; i < n_; ++i) {
      out[i] -= value(i, i);
    }
    return out;
  }

  double sum() const {
    return std::accumulate(data_.begin(), data_.end(), 0.0);
  }
};


inline sq_matrix element_mult(const sq_matrix& A, const sq_matrix& B) {
  auto n = A.size();

  sq_matrix out(n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      //out.set_value(i, j, temp);
      out.data_[ out.index(i, j) ] = A.value(i, j) * B.value(i, j);
    }
  }
  return out;
}

inline void element_mult(const sq_matrix& A,
                         const sq_matrix& B,
                         sq_matrix* out) {
  auto n = A.size();

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      out->data_[ out->index(i, j) ] = A.value(i, j) * B.value(i, j);
    }
  }
  return;
}

inline sq_matrix element_mult_one_minus(const sq_matrix& A,
                                        const sq_matrix& B) {

  auto n = A.size();

  sq_matrix out(n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      out.data_[ out.index(i, j)] = (1.0 - A.value(i, j)) * B.value(i, j);
    }
  }
  return out;
}

inline void element_mult_one_minus(const sq_matrix& A,
                                   const sq_matrix& B,
                                   sq_matrix* out) {
  auto n = A.size();

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      out->data_[ out->index(i, j)] = (1.0 - A.value(i, j)) * B.value(i, j);
    }
  }
  return;
}



inline void matrix_mult(const sq_matrix& A,
                        const sq_matrix& B,
                        sq_matrix* out) {
  auto n_ = A.n_;

  for (size_t i = 0; i < n_; ++i) {
    for (size_t j = 0; j < n_; ++j) {
      double s = 0.0;
      for (size_t k = 0; k < n_; ++k) {
        s += A.value(i, k) * B.value(k, j);
      }
      out->data_[ out->index(i, j)] = s;
    }
  }

  return;
}

inline void vector_mult(const sq_matrix& A,
                        const std::vector<double>& v,
                        std::vector<double>* out) {
  auto n_ = A.n_;
  for (size_t i = 0; i < n_; ++i) {
    double s = 0.0;
    for (size_t j = 0; j < n_; ++j) {
      s += v[j] * A.value(i, j);
    }
    (*out)[i] = s;
  }

  return;
}

inline void vector_view_mult(const sq_matrix& A,
                             const vector_view_t<const double>& v,
                             std::vector<double>* out)  {
  const std::size_t n = A.size();

  for (size_t i = 0; i < n; ++i) {
    double s = 0.0;
    for (size_t j = 0; j < n; ++j) {
      s += v[j] * A.value(i, j);
    }
    (*out)[i] = s;
  }

  return;
}



