//
//  C++ port of the island-area and per-capita rate functions used by the
//  time-dependent (ontogeny / sea-level) CS likelihood, i.e. of
//  island_area_vector(), island_area(), calc_Abeta(),
//  get_clado_rate_per_capita(), get_ext_rate_per_capita() and
//  get_immig_rate_per_capita() in R/DAISIE_rates.R and
//  R/DAISIE_loglik_CS_time.R.
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef DAISIE_TIME_RATES_H_INCLUDED
#define DAISIE_TIME_RATES_H_INCLUDED

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>


namespace daisie_time {

  // Parameters of the time-dependent CS model. Layout of the R vector
  // (first 18 elements of 'parsvec' of DAISIE_loglik_rhs_time*):
  //  1 lac0, 2 mu0, 3 K0, 4 gam0, 5 laa0, 6 d, 7 x,
  //  8-14 area_pars (max area, current area, proportional peak time,
  //       total island age, sea level amplitude, sea level frequency,
  //       island gradient angle),
  //  15 island_ontogeny, 16 sea_level, 17 total_time, 18 peak
  struct pars_t
  {
    static constexpr int size = 18;

    double lac0, mu0, K0, gam0, laa0, d, x;
    double area_pars[7];
    double island_ontogeny, sea_level, total_time, peak;

    template <typename It>
    explicit pars_t(It p)
      : lac0(p[0]), mu0(p[1]), K0(p[2]), gam0(p[3]), laa0(p[4]), d(p[5]), x(p[6]),
        island_ontogeny(p[14]), sea_level(p[15]), total_time(p[16]), peak(p[17])
    {
      for (int i = 0; i < 7; ++i) area_pars[i] = p[7 + i];
    }
  };


  // calc_Abeta()
  inline double calc_Abeta(double proptime, double proptime_max, double peak, double Amax)
  {
    const double f = proptime_max / (1.0 - proptime_max);
    const double a = f * peak / (1.0 + f);
    const double b = peak / (1.0 + f);
    return Amax * std::pow(proptime, a) * std::pow(1.0 - proptime, b) /
      (std::pow(a / (a + b), a) * std::pow(b / (a + b), b));
  }


  // island_area()
  inline double island_area(double timeval, const pars_t& p)
  {
    const double Amax = p.area_pars[0];
    const double Acurr = p.area_pars[1];
    const double proptime_max = p.area_pars[2];
    const double Tmax = p.area_pars[3];
    const double ampl = p.area_pars[4];
    const double freq = p.area_pars[5];
    const double theta = p.area_pars[6] * (M_PI / 180.0);
    const double proptime = timeval / Tmax;
    const double proptime_curr = p.total_time / Tmax;
    const bool ontogeny = (p.island_ontogeny == 1.0);
    const bool sea = (p.sea_level == 1.0);

    if (!ontogeny && !sea) {
      return 1.0;
    }
    if (ontogeny && !sea) {
      return calc_Abeta(proptime, proptime_max, p.peak, Amax);
    }
    if (!ontogeny && sea) {
      const double angular_freq = 2.0 * M_PI * freq;
      const double delta_sl = ampl * std::cos((proptime_curr - proptime) * angular_freq);
      const double r_curr = std::sqrt(Acurr / M_PI);
      const double h_curr = std::tan(theta) * r_curr;
      const double h_delta = std::max(0.0, h_curr - ampl + delta_sl);
      return M_PI * std::pow(h_delta / std::tan(theta), 2);
    }
    if (ontogeny && sea) {
      const double A_beta = calc_Abeta(proptime, proptime_max, p.peak, Amax);
      const double angular_freq = 2.0 * M_PI * freq;
      const double delta_sl = ampl * std::cos((proptime_curr - proptime) * angular_freq);
      const double r_curr = std::sqrt(A_beta / M_PI);
      const double h_curr = std::tan(theta) * r_curr;
      const double h_delta = std::max(0.0, h_curr - ampl + delta_sl);
      return M_PI * std::pow(h_delta / std::tan(theta), 2);
    }
    throw std::runtime_error("daisie_time::island_area: invalid island_ontogeny / sea_level");
  }


  // island_area_vector(): constant ontogeny (0 or NA) always gives area 1
  inline double island_area_vector(double timeval, const pars_t& p)
  {
    if (p.island_ontogeny == 0.0 || std::isnan(p.island_ontogeny)) {
      return 1.0;
    }
    return island_area(timeval, p);
  }


  // R's pmax(0, v, na.rm = TRUE): NaN -> 0
  inline double pmax0(double v)
  {
    return std::isnan(v) ? 0.0 : std::max(0.0, v);
  }


  // Per-capita rate vectors on the grid of species numbers
  // nn[j] = max(0, j - 2), j = 0 .. lnn - 1, lnn = lx + 4 + 2 * kk,
  // exactly as constructed in DAISIE_loglik_rhs_time*().
  // Storage is allocated once and reused for every rhs evaluation.
  class rate_vectors
  {
  public:
    rate_vectors(int lx, int kk)
      : lnn_(lx + 4 + 2 * kk),
        laavec(lnn_), lacvec(lnn_), muvec(lnn_), gamvec(lnn_), nn(lnn_)
    {
      for (int j = 0; j < lnn_; ++j) {
        nn[j] = std::max(0, j - 2);
      }
    }

    void update(double t, const pars_t& p)
    {
      const double A = island_area_vector(std::fabs(t), p);
      // get_ext_rate_per_capita(extcutoff = 1e6)
      const double mu = std::min(pmax0(p.mu0 * std::pow(A, -p.x)), 1000000.0);
      for (int j = 0; j < lnn_; ++j) {
        lacvec[j] = pmax0(p.lac0 * std::pow(A, p.d) * (1.0 - nn[j] / (p.K0 * A)));
        gamvec[j] = pmax0(p.gam0 * (1.0 - nn[j] / (p.K0 * A)));
        muvec[j] = mu;
        laavec[j] = p.laa0;
      }
    }

    int size() const { return lnn_; }

    int lnn_;
    std::vector<double> laavec, lacvec, muvec, gamvec, nn;
  };

}

#endif
