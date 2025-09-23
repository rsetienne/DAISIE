#include <memory>
#include <Rcpp.h>

#include "config.h"

#include "util.h"
#include "DAISIE_sim.h"
#include "island.h"

std::unique_ptr<DAISIE_sim> create_daisie_sim(double total_time,
                                              std::vector<double>& pars,
                                              Rcpp::List hyper_pars,
                                              Rcpp::List area_pars,
                                              int seed,
                                              int mainland_n,
                                              int island_ontogeny,
                                              int sea_level) {
  if (island_ontogeny == 0 && sea_level == 0) {
    auto is = island_static();
    return std::unique_ptr<DAISIE_sim>(new DAISIE_sim_impl<island_static>(pars,
                                                                          hyper_pars["d"],
                                                                          hyper_pars["x"],
                                                                          mainland_n,
                                                                          total_time,
                                                                          seed,
                                                                          is));
  } else if (island_ontogeny == 1 && sea_level == 0) {
    auto is_b = island_beta(area_pars["max_area"],
                            area_pars["total_island_age"],
                            area_pars["proportional_peak_t"]);

    return std::unique_ptr<DAISIE_sim>(new DAISIE_sim_impl<island_beta>(pars,
                                                                        hyper_pars["d"],
                                                                        hyper_pars["x"],
                                                                        mainland_n,
                                                                        total_time,
                                                                        seed,
                                                                        is_b));
  } else if (island_ontogeny == 0 && sea_level == 1) {
    auto is_a = island_angular(area_pars["sea_level_frequency"],
                               area_pars["total_island_age"],
                               total_time,
                               area_pars["amplitude"],
                               area_pars["current_area"],
                               area_pars["island_gradient_angle"]);

    return std::unique_ptr<DAISIE_sim>(new DAISIE_sim_impl<island_angular>(pars,
                                                                           hyper_pars["d"],
                                                                           hyper_pars["x"],
                                                                           mainland_n,
                                                                           total_time,
                                                                           seed,
                                                                           is_a));
  } else if (island_ontogeny == 1 && sea_level == 1) {
    auto is_b_a = island_beta_angular(area_pars["sea_level_frequency"],
                                      area_pars["total_island_age"],
                                      total_time,
                                      area_pars["amplitude"],
                                      area_pars["current_area"],
                                      area_pars["island_gradient_angle"],
                                      area_pars["max_area"],
                                      area_pars["proportional_peak_t"]);

    return std::unique_ptr<DAISIE_sim>(new DAISIE_sim_impl<island_beta_angular>(pars,
                                                                                hyper_pars["d"],
                                                                                hyper_pars["x"],
                                                                                mainland_n,
                                                                                total_time,
                                                                                seed,
                                                                                is_b_a));
  } else {
    throw "not viable island ontogeny selected";
  }
}

Rcpp::List daisie_sim_rcpp_wrap(double total_time,
                                std::vector<double>& pars,
                                Rcpp::List hyper_pars,
                                Rcpp::List  area_pars,
                                int seed,
                                int mainland_n,
                                int island_ontogeny,
                                int sea_level) {

  std::unique_ptr<DAISIE_sim> sim = create_daisie_sim(total_time,
                                                      pars,
                                                      hyper_pars,
                                                      area_pars,
                                                      seed,
                                                      mainland_n,
                                                      island_ontogeny,
                                                      sea_level);

  sim->run();

  auto island_spec_for_R = make_island_spec_for_R(sim->island_spec_);
  auto stt_table_for_R = make_stt_table_for_R(sim->stt_table_);
  // place holder code to get it to compile for now.
  Rcpp::List output =
    Rcpp::List::create(Rcpp::Named("island_spec") = island_spec_for_R,
                       Rcpp::Named("stt_table") = stt_table_for_R);
  return output;
}


//' Driver for the daisie sim time loop
//'
//' @description This is the rcpp function to run DAISIE simulations
//' @name daisie_sim_rcpp
//' @export daisie_sim_rcpp
//' @return List with the stt_table and the island_spec matrix
RcppExport SEXP daisie_sim_rcpp(SEXP total_timeSEXP, SEXP parsSEXP, SEXP hyper_parsSEXP, SEXP area_parsSEXP, SEXP seedSEXP, SEXP mainland_nSEXP, SEXP island_ontogenySEXP, SEXP sea_levelSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< double >::type total_time(total_timeSEXP);
  Rcpp::traits::input_parameter< std::vector<double>& >::type pars(parsSEXP);
  Rcpp::traits::input_parameter< Rcpp::List >::type hyper_pars(hyper_parsSEXP);
  Rcpp::traits::input_parameter< Rcpp::List >::type area_pars(area_parsSEXP);
  Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
  Rcpp::traits::input_parameter< int >::type mainland_n(mainland_nSEXP);
  Rcpp::traits::input_parameter< int >::type island_ontogeny(island_ontogenySEXP);
  Rcpp::traits::input_parameter< int >::type sea_level(sea_levelSEXP);
  rcpp_result_gen = Rcpp::wrap(daisie_sim_rcpp_wrap(total_time, pars, hyper_pars, area_pars, seed, mainland_n, island_ontogeny, sea_level));
  return rcpp_result_gen;
  END_RCPP
}







