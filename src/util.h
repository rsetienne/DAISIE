#ifndef TRAISIE_UTIL_H
#define TRAISIE_UTIL_H

#include <vector>
#include <string>
#include <sstream>
#include "island_spec.h"

Rcpp::NumericMatrix make_stt_table_for_R(const std::vector< std::array< double, 4 >>& stt_table) {
  int num_rows = stt_table.size();
  int num_cols = 4;
  Rcpp::NumericMatrix out(num_rows, num_cols);
  for (size_t i = 0; i < num_rows; ++i) {
    for (size_t j = 0; j < num_cols; ++j) {
      out(i, j) = stt_table[i][j];
    }
  }
  return out;
}

std::string get_string(const extinction_type& an) {
  if (an == extinction_type::clado_extinct) {
    return "Clado_extinct";
  }
  if (an == extinction_type::immig_parent) {
    return "Immig_parent";
  }
  if (an == extinction_type::none) {
    return "NA";
  }
  return "NA";
}

std::string get_string(const species_type& st) {
  if (st == species_type::I) return "I";
  if (st == species_type::A) return "A";
  if (st == species_type::C) return "C";
  return "NA";
}

std::string get_anc_string(const std::string& anc) {
  if (anc.empty()) return "NA";
  if (anc.size() < 1) return "NA";
  return anc;
}

std::string d_to_string(const double& x, double precision = 10) {
  // alternative to std::to_string that allows for precision.
  std::ostringstream ss;
  ss << std::fixed << std::setprecision(precision) << x;
  return ss.str();
}

std::string get_ext_time(double e_t) {
  if (e_t == -1) return "NA";
  return d_to_string(e_t);
}

Rcpp::StringMatrix make_island_spec_for_R(const island_spec& is) {
  int num_rows = is.size();
  int num_cols = 7;
  Rcpp::StringMatrix out(num_rows, num_cols);
  for (size_t i = 0; i < is.size(); ++i) {
    out(i, 0) = std::to_string(1 + static_cast<int>(is[i].id));
    out(i, 1) = std::to_string(1 + static_cast<int>(is[i].parent));
    out(i, 2) = d_to_string(is[i].colonisation_time);
    out(i, 3) = get_string(is[i].type_species);
    out(i, 4) = get_anc_string(is[i].anc_type);
    out(i, 5) = get_ext_time(is[i].get_extinction_time());
    out(i, 6) = get_string(is[i].ext_type);
  }
  return out;
}

bool match_motif(const std::string& anc,
                 const std::string& motif,
                 int n) {
  
  if (motif.empty()) return false;
  if (anc.empty()) return false;
  
 
  if (anc.compare(0, n, motif) == 0) return true;
  
  return false;
}


#endif
