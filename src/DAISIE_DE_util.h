//
//  Copyright (c) 2026 Thijs Janzen
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#include <cstdlib>    // std::getenv, std::atoi
#include <vector>
#include <array>
#include <chrono>
#include <string>
#include <utility>
#include <algorithm>
#include <memory>
#include "config.h"    // NOLINT [build/include_subdir]
#include <Rcpp.h>
#include "DAISIE_DE_odeint.h"    // NOLINT [build/include_subdir]
#include "secsse_loglik.h"    // NOLINT [build/include_subdir]


template <typename ODE>
std::vector<double> solve_branch(std::unique_ptr<ODE> od,
                                 const std::vector<double>& states,
                                 const std::array<double, 2>& forTime,
                                 const std::string& method,
                                 double atol,
                                 double rtol) {
  auto t0 = std::min(forTime[0], forTime[1]);
  auto t1 = std::max(forTime[0], forTime[1]);

  auto states_out = std::vector<double>(states.begin(), states.end());

  auto workhorse = Integrator<ODE, odeintcpp::no_normalization>(std::move(od), method, atol, rtol);

  workhorse(states_out, t0, t1);

  return states_out;
}

template <typename ODE>
std::vector<std::vector<double>> solve_branch_times(std::unique_ptr<ODE> od,
                                                    const std::vector<double>& states,
                                                    const std::vector<double>& forTime,
                                                    const std::string& method,
                                                    double atol,
                                                    double rtol) {
  std::vector< std::vector<double > > states_out;
  std::vector<double> times(forTime.begin(), forTime.end());

  auto workhorse = Integrator<ODE, odeintcpp::no_normalization>(std::move(od), method, atol, rtol);

  std::vector<double> states_in(states.begin(), states.end());

  workhorse(states_in, times, &states_out);

  return states_out;
}

struct solver {
public:
  solver(const Rcpp::NumericVector& brts,
         size_t missnumspec,
         double lambda_c,
         double lambda_a,
         double mu,
         double gamma,
         std::string method,
         double atol,
         double rtol) :
  lambda_c_(lambda_c),
  lambda_a_(lambda_a),
  mu_(mu),
  gamma_(gamma),
  method_(method),
  atol_(atol),
  rtol_(rtol),
  rho((brts.size() - 1) / (brts.size() - 1 + missnumspec))
  {
    t0 = brts[0];
    t1 = brts[1];
    if (brts.size() > 2) t2 = brts[2];
    tp = 0;
    ti = std::vector<double>(brts.begin(), brts.end());
    std::sort(ti.begin(), ti.end());
    ti.pop_back(); ti.pop_back(); //remove last two entries
  }

  double calculate_likelihood(const std::string& type, size_t stac) {
    if (type == "pEC")      return pEC(stac);
    else if (type == "pES") return pES(stac);
    else if (type == "pNE") return pNE(stac);
    else throw std::invalid_argument("type should be one of 'pEC', 'pES', or 'pNE'");
  }

private:
  double t0, t1, t2, tp;
  std::vector<double> ti;
  const double lambda_c_;
  const double lambda_a_;
  const double mu_;
  const double gamma_;
  const std::string method_;
  const double atol_;
  const double rtol_;
  const double rho; // sampling fraction
  double pEC(size_t stac) const {
    std::vector<double> initial_conditions1 = {rho, 0.0, 1 - rho, 1.0};
    if (stac == 3) initial_conditions1      = {rho, 1.0, 1 - rho, 0.0};
    auto ti2 = ti;
    ti2.insert(ti2.begin(), 0);
    auto solution0 = solve_branch_times(std::make_unique<loglik::interval2_EC>(lambda_c_, lambda_a_, mu_, gamma_), initial_conditions1, ti2, method_, atol_, rtol_);

    for (size_t i = 1; i < ti2.size(); ++i) {
      std::array<double, 2> time = {ti2[i - 1], ti2[i]};
      auto solution1 = solve_branch(std::make_unique<loglik::interval2_EC>(lambda_c_, lambda_a_, mu_, gamma_), initial_conditions1, time, method_, atol_, rtol_);
      enum state {DE, DM3, E, DA3};
      initial_conditions1 = {lambda_c_ * solution0[i][DE] * solution1[DE], // DE = pars1[1] * solution0[, "DE"][idx + 1] * solution1[, "DE"][2],
                             0,                                            // DM3
                             solution0[i][E],                              // E = solution0[, "E"][idx + 1],
                                         1                                             // DA3
                            };
    }

    std::vector<double> initial_conditions2;
    if (stac == 6) {
      enum state {DE, DM3, E, DA3};

      initial_conditions2 = { initial_conditions1[DE],                          // DE = initial_conditions1["DE"][[1]],
                              0,                                                // DM1 = 0,
                              initial_conditions1[DE] * solution0.back()[DA3],  // DM2 = initial_conditions1["DE"][[1]] * solution0[, "DA3"][length(ti) + 1],
                              solution0.back()[DA3],                            // DM3 = solution0[, "DM3"][length(ti) + 1],
                              initial_conditions1[E],                           // E = initial_conditions1["E"][[1]],
                              0,                                                // DA2 = 0,
                              solution0.back()[DA3]                             // DA3 = solution0[, "DA3"][length(ti) + 1])
                            };
    } else {
      enum state {DE, DM3, E, DA3};

      initial_conditions2 = { initial_conditions1[DE],                          // DE = initial_conditions1["DE"][[1]],
                              initial_conditions1[DE] * solution0.back()[DA3],  // DM2 = initial_conditions1["DE"][[1]] * solution0[, "DA3"][length(ti) + 1],
                              solution0.back()[DM3],                            // DM3 = solution0[, "DM3"][length(ti) + 1],
                              initial_conditions1[E],                           // E = initial_conditions1["E"][[1]],
                              solution0.back()[DA3]                             // DA3 = solution0[, "DA3"][length(ti) + 1])
      };
    }

    std::array<double, 2> time2 = {t2, t1};

    std::vector<double> solution2 = stac == 6 ?
    solve_branch(std::make_unique<loglik::interval3_ES>(lambda_c_, lambda_a_, mu_, gamma_), initial_conditions2, time2, method_, atol_, rtol_)
      :
      solve_branch(std::make_unique<loglik::interval2_ES>(lambda_c_, lambda_a_, mu_, gamma_), initial_conditions2, time2, method_, atol_, rtol_);

    std::vector<double> initial_conditions3;
    if (stac == 6) {
      // solution2 returns (interval3_ES): DE, DM1, DM2, DM3, E, DA2, DA3
      enum ls {DE, DM1, DM2, DM3, E, DA2, DA3}; // local state
      initial_conditions3 = {solution2[ls::DA2],  // DA1 = solution2[, "DA2"][[2]],
                             solution2[ls::DM1],  // DM1 = solution2[, "DM1"][[2]],
                             solution2[ls::E]     // E   = solution2[, "E"][[2]]
      };
    } else {
      // solution2 returns (interval2_ES): DE, DM2, DM3, E, DA3
      enum ls {DE, DM2, DM3, E, DA3}; // local state
      initial_conditions3 = {gamma_ * solution2[ls::DM2],  // DA1 = pars1[4] * solution2[, "DM2"][[2]],
                             gamma_ * solution2[ls::DM2], // DM1 = pars1[4] * solution2[, "DM2"][[2]],
                             solution2[ls::E]  // E = solution2[, "E"][[2]])
      };
    }

    std::array<double, 2> time3 = {t1, t0};

    auto solution3 = solve_branch(std::make_unique<loglik::interval4>(lambda_c_, lambda_a_, mu_, gamma_), initial_conditions3, time3, method_, atol_, rtol_);
    // interval4 returns: DA1, DM1, E
    enum class ls {DA1, DM1, E}; // local state
    auto prob = solution3[static_cast<size_t>(ls::DA1)];

    return std::log(prob);
  }

  double pES(size_t stac) const {

    std::array<double, 2> time1 = {tp, t1};

    std::vector<double> initial_conditions1 = {rho, 0.0, 0.0, 1 - rho, 1.0};
    if (stac == 3)      initial_conditions1 = {rho, 0.0, 1.0, 1 - rho, 0.0};
    if (stac == 5)      initial_conditions1 = {rho, 0.0, 0.0, 0.0, 1 - rho, 0.0, 1.0};
    if (stac == 9)      {
      initial_conditions1 = {1.0, 0.0, 0.0, 0.0, 1.0};
      time1 = {tp, t2};
    }

    auto solution1 = stac == 5 ?
      solve_branch(std::make_unique<loglik::interval3_ES>(lambda_c_, lambda_a_, mu_, gamma_), initial_conditions1, time1, method_, atol_, rtol_)
      :
      solve_branch(std::make_unique<loglik::interval2_ES>(lambda_c_, lambda_a_, mu_, gamma_), initial_conditions1, time1, method_, atol_, rtol_);

    std::vector<double> initial_conditions3;
    if (stac == 9) {
      enum state {DE, DM2, DM3, E, DA3};
      std::vector<double> initial_conditions2 = { solution1[DE],  // DE  = solution1[, "DE"][[2]],
                                                  0,              // DM1 = 0
                                                  solution1[DM2], //      DM2 = solution1[, "DM2"][[2]]
                                                  solution1[DM3], // DM3 = solution1[, "DM3"][[2]],
                                                  solution1[E],   // E   =  solution1[, "E"][[2]],
                                                  0,              // DA2 = 0
                                                  solution1[DA3]  // DA3 = solution1[, "DA3"][[2]]
      };

      std::array<double, 2> time2 = {t2, t1};

      auto solution2 = solve_branch(std::make_unique<loglik::interval3_ES>(lambda_c_, lambda_a_, mu_, gamma_), initial_conditions2, time2, method_, atol_, rtol_);
      enum class ls {DE, DM1, DM2, DM3, E, DA2, DA3}; // local state
      initial_conditions3 = {solution2[static_cast<size_t>(ls::DA2)],  // DA1 = solution2[, "DA2"][[2]],
                             solution2[static_cast<size_t>(ls::DM1)],  // DM1 = solution2[, "DM1"][[2]],
                             solution2[static_cast<size_t>(ls::E)]     // E   = solution2[, "E"][[2]])
      };
    } else if (stac == 5) {
      enum state {DE, DM1, DM2, DM3, E, DA2, DA3};
      initial_conditions3 = {solution1[DA2], // DA1 = solution2[, "DA2"][[2]],
                             solution1[DM1], // DM1 = solution2[, "DM1"][[2]],
                             solution1[E]    // E   = solution2[, "E"][[2]])
      };
    } else {
      enum state {DE, DM2, DM3, E, DA3};
      initial_conditions3 = {gamma_ * solution1[DM2],  // DA1 = pars1[4] * solution1[, "DM2"][[2]],
                             gamma_ * solution1[DM2],  // DM1 = pars1[4] * solution1[, "DM2"][[2]],
                             solution1[E]              // E = solution1[, "E"][[2]])
      };
    }

    std::array<double, 2> time3 = {t1, t0};

    auto solution3 = solve_branch(std::make_unique<loglik::interval4>(lambda_c_, lambda_a_, mu_, gamma_), initial_conditions3, time3, method_, atol_, rtol_);

    // interval4 returns: DA1, DM1, E
    enum class ls {DA1, DM1, E}; // local state
    auto prob = solution3[static_cast<size_t>(ls::DA1)];

    return std::log(prob);
  }

  double pNE(size_t stac) const {

    std::array<double, 2> time1 = {tp, t1};

    std::vector<double> initial_conditions1 = {1.0, 0.0};
    if (stac == 1)      initial_conditions1 = {0.0, 1.0, 0.0, 0.0};
    if (stac == 8) time1 = {tp, t2};

    auto solution1 = stac == 1 ?
      solve_branch(std::make_unique<loglik::interval3_NE>(lambda_c_, lambda_a_, mu_, gamma_), initial_conditions1, time1, method_, atol_, rtol_)
      :
      solve_branch(std::make_unique<loglik::interval2_NE>(lambda_c_, lambda_a_, mu_, gamma_), initial_conditions1, time1, method_, atol_, rtol_);



    std::vector<double> initial_conditions3;
    if (stac == 8) {
      enum class ls1 {DM2, E}; // local state
      std::vector<double>     initial_conditions2 = { 0,                                        // DM1 = 0
                                                      solution1[static_cast<size_t>(ls1::DM2)], //      DM2 = solution1[, "DM2"][[2]]
                                                      solution1[static_cast<size_t>(ls1::E)],   // E   =  solution1[, "E"][[2]],
                                                      0,                                        // DA2 = 0
      };

      std::array<double, 2> time2 = {t2, t1};

      auto solution2 = solve_branch(std::make_unique<loglik::interval3_NE>(lambda_c_, lambda_a_, mu_, gamma_), initial_conditions2, time2, method_, atol_, rtol_);
      enum class ls {DM1, DM2, E, DA2}; // local state
      initial_conditions3 = {solution2[static_cast<size_t>(ls::DA2)], // DA1 = solution2[, "DA2"][[2]],
                             solution2[static_cast<size_t>(ls::DM1)], // DM1 = solution2[, "DM1"][[2]],
                             solution2[static_cast<size_t>(ls::E)]    // E   = solution2[, "E"][[2]])
      };
    } else if (stac == 1) {
      enum class ls {DM1, DM2, E, DA2}; // local state
      initial_conditions3 = {solution1[static_cast<size_t>(ls::DA2)],  // DA1 = solution2[, "DA2"][[2]],
                             solution1[static_cast<size_t>(ls::DM1)],  // DM1 = solution2[, "DM1"][[2]],
                             solution1[static_cast<size_t>(ls::E)]     // E   = solution2[, "E"][[2]])
      };
    } else {
      enum class ls {DM2, E}; // local state
      initial_conditions3 = {gamma_ * solution1[static_cast<size_t>(ls::DM2)],  // DA1 = pars1[4] * solution1[, "DM2"][[2]],
                             gamma_ * solution1[static_cast<size_t>(ls::DM2)],  // DM1 = pars1[4] * solution1[, "DM2"][[2]],
                             solution1[static_cast<size_t>(ls::E)]              // E = solution1[, "E"][[2]])
      };
    }

    std::array<double, 2> time3 = {t1, t0};

    auto solution3 = solve_branch(std::make_unique<loglik::interval4>(lambda_c_, lambda_a_, mu_, gamma_), initial_conditions3, time3, method_, atol_, rtol_);

    // interval4 returns: DA1, DM1, E
    enum class ls {DA1, DM1, E}; // local state
    auto prob = solution3[static_cast<size_t>(ls::DA1)];

    return std::log(prob);
  }
};
