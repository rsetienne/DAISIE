#include <Rcpp.h>
#include <boost/numeric/odeint.hpp>
#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>
#include <cmath>

namespace {

  using State = std::vector<double>;

  struct SeriesRHS {
    int K;
    int n;
    int interval;
    unsigned long evaluations = 0;

    State lc, mu, la, g, loss;
    double G;

    Rcpp::NumericMatrix q;
    Rcpp::NumericMatrix p;

    std::vector<std::vector<int>> left;
    std::vector<std::vector<int>> right;

    SeriesRHS(
      const Rcpp::List& parameter,
      const Rcpp::NumericVector& gamma,
      const Rcpp::List& pairs,
      int interval_number)
    : K(pairs.size()),
    n(Rcpp::as<Rcpp::NumericVector>(parameter[0]).size()),
    interval(interval_number),
    lc(Rcpp::as<State>(parameter[0])),
    mu(Rcpp::as<State>(parameter[1])),
    la(Rcpp::as<State>(parameter[3])),
    g(gamma.begin(), gamma.end()),
    loss(n),
    G(std::accumulate(g.begin(), g.end(), 0.0)),
    q(Rcpp::as<Rcpp::NumericMatrix>(parameter[4])),
    p(Rcpp::as<Rcpp::NumericMatrix>(parameter[5])),
    left(K),
    right(K) {

      if (K < 1 || n < 1 ||
          mu.size() != lc.size() ||
          la.size() != lc.size() ||
          g.size() != lc.size() ||
          q.nrow() != n || q.ncol() != n ||
          p.nrow() != n || p.ncol() != n) {
        Rcpp::stop("Invalid series parameter dimensions.");
      }

      if (interval < 1 || interval > 4) {
        Rcpp::stop("Invalid series interval.");
      }

      for (int s = 0; s < n; ++s) {
        loss[s] = lc[s] + mu[s] + la[s];

        for (int j = 0; j < n; ++j) {
          loss[s] += q(s, j);
        }
      }

      for (int k = 0; k < K; ++k) {
        Rcpp::List pair = pairs[k];

        Rcpp::IntegerVector a = pair["left"];
        Rcpp::IntegerVector b = pair["right"];

        if (a.size() != b.size()) {
          Rcpp::stop("Invalid coefficient multiplication pairs.");
        }

        for (int i = 0; i < a.size(); ++i) {
          if (a[i] == NA_INTEGER || b[i] == NA_INTEGER ||
              a[i] < 1 || a[i] > K ||
              b[i] < 1 || b[i] > K) {
            Rcpp::stop("Invalid coefficient row index.");
          }

          // Convert R's row numbers to C++ indices.
        left[k].push_back(a[i] - 1);
        right[k].push_back(b[i] - 1);
      }
    }
  }

  // R stores the coefficient matrix column by column.
  int position(int column, int degree) const {
    return column * K + degree;
  }

  double value(
      const State& x,
      int column,
      int degree) const {

    return x[position(column, degree)];
  }

  double product(
      const State& x,
      int column_a,
      int column_b,
      int degree) const {

    double result = 0.0;

    for (std::size_t i = 0; i < left[degree].size(); ++i) {
      result +=
        value(x, column_a, left[degree][i]) *
        value(x, column_b, right[degree][i]);
    }

    return result;
  }

  // Multiply a precomputed f_s series by a DA series.
  double product_f(
      const State& f,
      int state,
      const State& x,
      int da_column,
      int degree) const {

    double result = 0.0;

    for (std::size_t i = 0; i < left[degree].size(); ++i) {
      result +=
        f[position(state, left[degree][i])] *
        value(x, da_column, right[degree][i]);
    }

    return result;
  }

  void operator()(
      const State& x,
      State& dx,
      double /* time */) {

    ++evaluations;

    if (evaluations % 1024 == 0) {
      Rcpp::checkUserInterrupt();
    }

    if (evaluations > 2000000) {
      Rcpp::stop("Too many series RHS evaluations.");
    }

    dx.assign(x.size(), 0.0);

    const int DA2 = 5 * n;
    const int DA3 = 5 * n + 1;

    State f(K * n, 0.0);
    State h(K * n, 0.0);

    // Calculate E and DE derivatives and the f/h source series.
    for (int s = 0; s < n; ++s) {
      const int DE = s;
      const int E = 4 * n + s;
      const double qsum = loss[s] - lc[s] - mu[s] - la[s];

      for (int k = 0; k < K; ++k) {
        const double e = value(x, E, k);
        const double de = value(x, DE, k);

        const double ee = product(x, E, E, k);
        const double dee = product(x, DE, E, k);

        double q_e = 0.0;
        double q_de = 0.0;
        double pq_e = 0.0;
        double pq_de = 0.0;

        for (int j = 0; j < n; ++j) {
          const double qsj = q(s, j);
          const double pqsj = p(s, j) * qsj;

          q_e += qsj * value(x, 4 * n + j, k);
          q_de += qsj * value(x, j, k);

          pq_e += pqsj * value(x, 4 * n + j, k);
          pq_de += pqsj * value(x, j, k);
        }

        const double constant_mu = k == 0 ? mu[s] : 0.0;

        f[position(s, k)] =
          constant_mu + la[s] * e + lc[s] * ee + pq_e;

        h[position(s, k)] =
          la[s] * de + 2.0 * lc[s] * dee + pq_de;

        dx[position(E, k)] =
          constant_mu -
          (lc[s] + mu[s] + qsum) * e +
          lc[s] * ee +
          q_e;

        if (interval != 4) {
          dx[position(DE, k)] =
            -(lc[s] + mu[s] + qsum) * de +
            2.0 * lc[s] * dee +
            q_de;
        }
      }
    }

    // Calculate DM and DA derivatives.
    for (int k = 0; k < K; ++k) {
      double gM1 = 0.0;
      double gM2 = 0.0;
      double gM3 = 0.0;

      for (int j = 0; j < n; ++j) {
        gM1 += g[j] * value(x, n + j, k);
        gM2 += g[j] * value(x, 2 * n + j, k);
        gM3 += g[j] * value(x, 3 * n + j, k);
      }

      for (int s = 0; s < n; ++s) {
        const int M1 = n + s;
        const int M2 = 2 * n + s;
        const int M3 = 3 * n + s;

        double opqM1 = 0.0;
        double opqM2 = 0.0;
        double opqM3 = 0.0;

        for (int j = 0; j < n; ++j) {
          const double rate = (1.0 - p(s, j)) * q(s, j);

          opqM1 += rate * value(x, n + j, k);
          opqM2 += rate * value(x, 2 * n + j, k);
          opqM3 += rate * value(x, 3 * n + j, k);
        }

        if (interval != 4) {
          dx[position(M3, k)] =
            -(loss[s] + G) * value(x, M3, k) +
            product_f(f, s, x, DA3, k) +
            opqM3 +
            gM3;
        }

        if (interval == 2) {
          dx[position(M2, k)] =
            -(loss[s] + G) * value(x, M2, k) +
            product_f(h, s, x, DA3, k) +
            opqM2;
        }

        if (interval == 3) {
          dx[position(M1, k)] =
            -(loss[s] + G) * value(x, M1, k) +
            product_f(f, s, x, DA2, k) +
            opqM1 +
            gM2;

          dx[position(M2, k)] =
            -(loss[s] + G) * value(x, M2, k) +
            product_f(f, s, x, DA2, k) +
            product_f(h, s, x, DA3, k) +
            opqM2 +
            gM2;
        }

        if (interval == 4) {
          // The DA2 column stores DA1 during interval4.
          dx[position(M1, k)] =
            -(loss[s] + G) * value(x, M1, k) +
            product_f(f, s, x, DA2, k) +
            opqM1 +
            gM1;
        }
      }

      if (interval != 4) {
        dx[position(DA3, k)] =
          -G * value(x, DA3, k) + gM3;
      }

      if (interval == 3) {
        dx[position(DA2, k)] =
          -G * value(x, DA2, k) + gM2;
      }

      if (interval == 4) {
        dx[position(DA2, k)] =
          -G * value(x, DA2, k) + gM1;
      }
    }

    for (double derivative : dx) {
      if (!std::isfinite(derivative)) {
        Rcpp::stop("Non-finite series derivative.");
      }
    }
  }
};

} // anonymous namespace

RcppExport SEXP TRAISIE_series_cpp(
    SEXP initialSEXP,
    SEXP parameterSEXP,
    SEXP gammaSEXP,
    SEXP pairsSEXP,
    SEXP intervalSEXP,
    SEXP durationSEXP,
    SEXP atolSEXP,
    SEXP rtolSEXP) {

  BEGIN_RCPP

  State x = Rcpp::as<State>(initialSEXP);

  Rcpp::List parameter(parameterSEXP);
  Rcpp::NumericVector gamma(gammaSEXP);
  Rcpp::List pairs(pairsSEXP);

  const int interval = Rcpp::as<int>(intervalSEXP);
  const double duration = Rcpp::as<double>(durationSEXP);
  const double atol = Rcpp::as<double>(atolSEXP);
  const double rtol = Rcpp::as<double>(rtolSEXP);

  if (parameter.size() < 6) {
    Rcpp::stop("Six parameter components are required.");
  }

  if (!std::isfinite(duration) || duration < 0.0 ||
      !std::isfinite(atol) || atol <= 0.0 ||
      !std::isfinite(rtol) || rtol <= 0.0) {
    Rcpp::stop("Invalid integration duration or tolerances.");
  }

  SeriesRHS system(parameter, gamma, pairs, interval);

  const std::size_t expected =
    static_cast<std::size_t>(system.K) * (5 * system.n + 2);

  if (x.size() != expected) {
    Rcpp::stop("Invalid initial coefficient vector length.");
  }

  for (double coefficient : x) {
    if (!std::isfinite(coefficient)) {
      Rcpp::stop("Non-finite initial coefficient.");
    }
  }

  if (duration > 0.0) {
    using namespace boost::numeric::odeint;
    using Stepper = runge_kutta_fehlberg78<State>;

    auto controlled = make_controlled(atol, rtol, Stepper());

    integrate_adaptive(
      controlled,
      std::ref(system),
      x,
      0.0,
      duration,
      std::min(0.01, duration)
    );
  }

  for (double coefficient : x) {
    if (!std::isfinite(coefficient)) {
      Rcpp::stop("Non-finite final coefficient.");
    }
  }

  return Rcpp::wrap(x);

  END_RCPP
}
