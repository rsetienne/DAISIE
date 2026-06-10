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

class sq_matrix {
 private:
  std::vector<double> data_;
  const size_t n_;

 public:
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

  sq_matrix(const sq_matrix& other) : n_(other.get_n()) {
    data_ = other.get_data_();
  }

  explicit sq_matrix(size_t n) : n_(n) {}

  double value(size_t i, size_t j) const {
    return data_[i * n_ + j];
  }

  std::vector<double> get_data_() const {
    return data_;
  }

  size_t get_n() const {
    return n_;
  }

  void set_value(size_t i, size_t j, double val) {
    data_[i * n_ + j] = val;
  }


  sq_matrix operator*(const sq_matrix& other) const {
    assert(n_ == other.get_n());

    sq_matrix out(n_);

    for (size_t i = 0; i < n_; ++i) {
      for (size_t j = 0; j < n_; ++j) {
        double s = 0.0;
        for (size_t k = 0; k < n_; ++k) {
          s += value(i, k) * other.value(k, j);
        }
        out.set_value(i, j, s);
      }
    }

    return out;
  }

  std::vector<double> operator*(const std::vector<double>& v) const {
    std::vector<double> out(n_);

    assert(v.size() == n_);

    for (size_t i = 0; i < n_; ++i) {
      double s = 0.0;
      for (size_t j = 0; j < n_; ++j) {
        auto a = v[j];
        auto b = value(i, j);

        s += a * b;
      }
      out[i] = s;
    }

    return out;
  }

  std::vector<double> operator*(const vector_view_t<const double>& v) const {
    std::vector<double> out(n_);

    assert(v.size() == n_);

    for (size_t i = 0; i < n_; ++i) {
      double s = 0.0;
      for (size_t j = 0; j < n_; ++j) {
        auto a = v[j];
        auto b = value(i, j);

        s += a * b;
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
