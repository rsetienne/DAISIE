// C++ engine for the DAISIE simulation cores.
//
// This engine coexists with the R implementations of the simulation cores
// (DAISIE_sim_core_cr, _cr_shift, _time_dep, _trait_dep and DAISIE_SR_sim_core)
// and is selected by use_rcpp = TRUE. It implements the Monte Carlo loop of
// each core (the part that dominates the run time) and returns the island_spec
// matrix and the stt_table; the finalisation (DAISIE_create_island etc.) stays
// in R, so both implementations return output of identical format.
//
// The engine uses R's random number generator and reproduces the way R draws
// from it (stats::rexp, stats::runif, sample() with and without prob), in the
// same order as the R code. Hence, for the same seed (set.seed()), the C++ and
// R implementations return identical simulations.

#include <Rcpp.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

namespace {

const double kInf = std::numeric_limits<double>::infinity();

// ---------------------------------------------------------------- RNG ------

// stats::rexp(1, rate)
inline double rexp1(double rate) {
  return (1.0 / rate) * exp_rand();
}

// stats::runif(1)
inline double runif1() {
  return unif_rand();
}

// sample.int(n, 1) (zero-based result)
inline int sample_index(int n) {
  return static_cast<int>(R_unif_index(static_cast<double>(n)));
}

// DDD::sample2(x, 1) for a vector x: a vector of length one is doubled, hence
// it always consumes a random number.
inline int sample2_index(int n) {
  if (n < 1) Rcpp::stop("cannot sample from an empty set");
  return sample_index(n < 2 ? 2 : n) % n;
}

// sample(1:n, 1, prob = p), identical to R's ProbSampleReplace (n < 200)
int sample_prob(const std::vector<double>& probs) {
  const int n = static_cast<int>(probs.size());
  std::vector<double> p(probs);
  std::vector<int> perm(n);
  double sum = 0.0;
  int npos = 0;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(p[i]) || p[i] < 0.0) {
      Rcpp::stop("NA in probability vector");
    }
    if (p[i] > 0.0) {
      ++npos;
      sum += p[i];
    }
  }
  if (npos == 0) {
    Rcpp::stop("too few positive probabilities");
  }
  for (int i = 0; i < n; ++i) {
    p[i] /= sum;
    perm[i] = i + 1;
  }
  revsort(p.data(), perm.data(), n);
  for (int i = 1; i < n; ++i) {
    p[i] += p[i - 1];
  }
  const double rU = unif_rand();
  int j;
  for (j = 0; j < n - 1; ++j) {
    if (rU <= p[j]) break;
  }
  return perm[j];
}

// ------------------------------------------------------- island state ------

// as.character(<double>) keeps 15 significant digits. The R code stores times
// in a character matrix; we round in the same way to obtain identical results.
inline double round15(double x) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.15g", x);
  return std::strtod(buf, nullptr);
}

enum Status { kNA = 0, kImmigParent = 1, kCladoExtinct = 2 };

struct Sp {
  int id;
  int anc;
  double col;       // colonisation time
  char type;        // 'I', 'A' or 'C'
  std::string br;   // branch code
  bool has_br;
  double brt;       // branching time
  bool has_brt;
  int status;
  int trait;        // 0 (no trait), 1 or 2
};

struct Island {
  std::vector<Sp> sp;
  int maxspecID;
  bool trait_mode;
  std::vector<std::vector<double>> stt;  // rows of the stt_table

  int n() const { return static_cast<int>(sp.size()); }
  int count_type(char t, int trait = 0) const {
    int k = 0;
    for (const auto& s : sp) {
      if (s.type == t && (trait == 0 || s.trait == trait)) ++k;
    }
    return k;
  }
  int count_trait(int trait) const {
    int k = 0;
    for (const auto& s : sp) {
      if (s.trait == trait) ++k;
    }
    return k;
  }
  std::vector<int> which_trait(int trait) const {
    std::vector<int> v;
    for (int i = 0; i < n(); ++i) {
      if (sp[i].trait == trait) v.push_back(i);
    }
    return v;
  }

  void push_stt(double time_before_present) {
    if (trait_mode) {
      stt.push_back({time_before_present,
                     static_cast<double>(count_type('I', 1)),
                     static_cast<double>(count_type('A', 1)),
                     static_cast<double>(count_type('C', 1)),
                     static_cast<double>(count_type('I', 2)),
                     static_cast<double>(count_type('A', 2)),
                     static_cast<double>(count_type('C', 2))});
    } else {
      stt.push_back({time_before_present,
                     static_cast<double>(count_type('I')),
                     static_cast<double>(count_type('A')),
                     static_cast<double>(count_type('C'))});
    }
  }

  // ---- events (mirror DAISIE_sim_update_state_*) ----

  void immigrate(int colonist, double timeval, int trait) {
    int there = -1;
    for (int i = 0; i < n(); ++i) {
      if (sp[i].id == colonist) {
        there = i;
        break;
      }
    }
    Sp s;
    s.id = colonist;
    s.anc = colonist;
    s.col = round15(timeval);
    s.type = 'I';
    s.has_br = false;
    s.has_brt = false;
    s.brt = 0.0;
    s.status = kNA;
    s.trait = trait;
    if (there < 0) {
      sp.push_back(s);
    } else {
      sp[there] = s;
    }
  }

  void extinct(int ext) {
    if (sp[ext].type == 'I' || sp[ext].type == 'A') {
      sp.erase(sp.begin() + ext);
      return;
    }
    // type C: find sisters
    std::vector<int> sisters;
    for (int i = 0; i < n(); ++i) {
      if (sp[i].anc == sp[ext].anc && sp[i].col == sp[ext].col) {
        sisters.push_back(i);
      }
    }
    std::vector<int> survivors;
    for (int i : sisters) {
      if (i != ext) survivors.push_back(i);
    }
    if (sisters.size() == 2) {
      Sp& s = sp[survivors[0]];
      s.type = 'A';
      s.has_br = false;
      s.has_brt = false;
      s.status = kCladoExtinct;
      sp.erase(sp.begin() + ext);
    } else if (sisters.size() >= 3) {
      const std::string& b = sp[ext].br;
      const size_t nsplits = b.size();
      const char mostrecent = b[nsplits - 1];
      const char sister_recent = (mostrecent == 'B') ? 'A' : 'B';
      std::string motif = b.substr(0, nsplits - 1) + sister_recent;
      std::vector<int> possible;
      for (int i : survivors) {
        if (sp[i].br.compare(0, nsplits, motif) == 0 &&
            sp[i].br.size() >= nsplits) {
          possible.push_back(i);
        }
      }
      if (mostrecent == 'A' && !possible.empty()) {
        double mn = kInf;
        for (int i : possible) mn = std::min(mn, sp[i].brt);
        for (int i : possible) {
          if (sp[i].brt == mn) {
            sp[i].brt = sp[ext].brt;
            sp[i].has_brt = sp[ext].has_brt;
          }
        }
      }
      for (int i : possible) {
        // remove the character at position nsplits (1-based)
        sp[i].br = sp[i].br.substr(0, nsplits - 1) + sp[i].br.substr(nsplits);
      }
      sp.erase(sp.begin() + ext);
    }
  }

  void anagenesis(int idx) {
    ++maxspecID;
    sp[idx].type = 'A';
    sp[idx].id = maxspecID;
    sp[idx].status = kImmigParent;
  }

  void cladogenesis(int idx, double timeval) {
    const double t = round15(timeval);
    Sp daughter;
    if (sp[idx].type == 'C') {
      const std::string old = sp[idx].br;
      sp[idx].id = maxspecID + 1;
      sp[idx].br = old + "A";
      sp[idx].status = kNA;
      daughter.br = old + "B";
    } else {
      sp[idx].type = 'C';
      sp[idx].id = maxspecID + 1;
      sp[idx].br = "A";
      sp[idx].has_br = true;
      sp[idx].brt = sp[idx].col;
      sp[idx].has_brt = true;
      sp[idx].status = kNA;
      daughter.br = "B";
    }
    daughter.id = maxspecID + 2;
    daughter.anc = sp[idx].anc;
    daughter.col = sp[idx].col;
    daughter.type = 'C';
    daughter.has_br = true;
    daughter.brt = t;
    daughter.has_brt = true;
    daughter.status = kNA;
    daughter.trait = sp[idx].trait;
    sp.push_back(daughter);
    maxspecID += 2;
  }
};

Island island_from_R(const Rcpp::CharacterMatrix& m, int maxspecID, bool trait) {
  Island is;
  is.maxspecID = maxspecID;
  is.trait_mode = trait;
  for (int i = 0; i < m.nrow(); ++i) {
    Sp s;
    s.id = std::atoi(Rcpp::as<std::string>(m(i, 0)).c_str());
    s.anc = std::atoi(Rcpp::as<std::string>(m(i, 1)).c_str());
    s.col = std::atof(Rcpp::as<std::string>(m(i, 2)).c_str());
    s.type = Rcpp::as<std::string>(m(i, 3))[0];
    s.has_br = !Rcpp::CharacterVector::is_na(m(i, 4));
    s.br = s.has_br ? Rcpp::as<std::string>(m(i, 4)) : std::string();
    s.has_brt = !Rcpp::CharacterVector::is_na(m(i, 5));
    s.brt = s.has_brt ? std::atof(Rcpp::as<std::string>(m(i, 5)).c_str()) : 0.0;
    s.status = kNA;
    if (!Rcpp::CharacterVector::is_na(m(i, 6))) {
      s.status = (Rcpp::as<std::string>(m(i, 6)) == "Immig_parent")
                     ? kImmigParent
                     : kCladoExtinct;
    }
    s.trait = 0;
    if (trait && m.ncol() >= 8) {
      s.trait = std::atoi(Rcpp::as<std::string>(m(i, 7)).c_str());
    }
    is.sp.push_back(s);
  }
  return is;
}

std::string fmt15(double x) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.15g", x);
  return std::string(buf);
}

Rcpp::CharacterMatrix island_to_R(const Island& is) {
  const int nc = is.trait_mode ? 8 : 7;
  Rcpp::CharacterMatrix m(is.n(), nc);
  for (int i = 0; i < is.n(); ++i) {
    const Sp& s = is.sp[i];
    m(i, 0) = std::to_string(s.id);
    m(i, 1) = std::to_string(s.anc);
    m(i, 2) = fmt15(s.col);
    m(i, 3) = std::string(1, s.type);
    m(i, 4) = s.has_br ? Rcpp::String(s.br) : Rcpp::String(NA_STRING);
    m(i, 5) = s.has_brt ? Rcpp::String(fmt15(s.brt)) : Rcpp::String(NA_STRING);
    if (s.status == kNA) {
      m(i, 6) = NA_STRING;
    } else {
      m(i, 6) = (s.status == kImmigParent) ? "Immig_parent" : "Clado_extinct";
    }
    if (is.trait_mode) m(i, 7) = std::to_string(s.trait);
  }
  return m;
}

Rcpp::NumericMatrix stt_to_R(const Island& is) {
  const int nc = is.trait_mode ? 7 : 4;
  Rcpp::NumericMatrix m(static_cast<int>(is.stt.size()), nc);
  for (int i = 0; i < m.nrow(); ++i) {
    for (int j = 0; j < nc; ++j) m(i, j) = is.stt[i][j];
  }
  return m;
}

// ------------------------------------------------------------- area --------

struct AreaPars {
  double max_area, current_area, proptime_max, total_island_age;
  double amplitude, frequency, theta_deg;
};

AreaPars area_from_R(const Rcpp::List& l) {
  AreaPars a;
  a.max_area = Rcpp::as<double>(l["max_area"]);
  a.current_area = Rcpp::as<double>(l["current_area"]);
  a.proptime_max = Rcpp::as<double>(l["proportional_peak_t"]);
  a.total_island_age = Rcpp::as<double>(l["total_island_age"]);
  a.amplitude = Rcpp::as<double>(l["sea_level_amplitude"]);
  a.frequency = Rcpp::as<double>(l["sea_level_frequency"]);
  a.theta_deg = Rcpp::as<double>(l["island_gradient_angle"]);
  return a;
}

double calc_Abeta(double proptime, double proptime_max, double peak,
                  double Amax) {
  const double f = proptime_max / (1.0 - proptime_max);
  const double a = f * peak / (1.0 + f);
  const double b = peak / (1.0 + f);
  return Amax * std::pow(proptime, a) * std::pow(1.0 - proptime, b) /
         (std::pow(a / (a + b), a) * std::pow(b / (a + b), b));
}

double island_area(double timeval, double total_time, const AreaPars& ap,
                   double peak, int ontogeny, int sea_level) {
  const double Tmax = ap.total_island_age;
  const double proptime = timeval / Tmax;
  const double proptime_curr = total_time / Tmax;
  const double theta = ap.theta_deg * (M_PI / 180.0);
  if (ontogeny == 0 && sea_level == 0) {
    return 1.0;
  }
  if (ontogeny == 1 && sea_level == 0) {
    return calc_Abeta(proptime, ap.proptime_max, peak, ap.max_area);
  }
  const double angular_freq = 2.0 * M_PI * ap.frequency;
  const double delta_sl =
      ap.amplitude * std::cos((proptime_curr - proptime) * angular_freq);
  double base_area;
  if (ontogeny == 0) {
    base_area = ap.current_area;
  } else {
    base_area = calc_Abeta(proptime, ap.proptime_max, peak, ap.max_area);
  }
  const double r_curr = std::sqrt(base_area / M_PI);
  const double h_curr = std::tan(theta) * r_curr;
  const double h_delta = std::max(0.0, h_curr - ap.amplitude + delta_sl);
  return M_PI * std::pow(h_delta / std::tan(theta), 2.0);
}

// ------------------------------------------------------------ rates --------

// pmax(0, v, na.rm = TRUE)
inline double pmax0(double v) {
  return (std::isnan(v) || v < 0.0) ? 0.0 : v;
}

struct Rates4 {
  double immig, ext, ana, clado;
};

double ext_rate(double mu, double x, double extcutoff, int num_spec, double A) {
  double pc = pmax0(mu * std::pow(A, -x));
  pc = std::min(pc, extcutoff);
  double r = num_spec * pc;
  return std::min(r, extcutoff);
}

double clado_rate(double lac, double d, int num_spec, double K, double A) {
  double pc = pmax0(lac * std::pow(A, d) * (1.0 - num_spec / (K * A)));
  return num_spec * pc;
}

double immig_rate(double gam, double A, int num_spec, double K,
                  double mainland_n) {
  return mainland_n * pmax0(gam * (1.0 - (num_spec / (K * A))));
}

Rates4 rates4(double gam, double laa, double lac, double mu, double d, double x,
              double extcutoff, double K, int num_spec, int num_imm,
              double mainland_n, double A) {
  Rates4 r;
  r.immig = immig_rate(gam, A, num_spec, K, mainland_n);
  r.ext = ext_rate(mu, x, extcutoff, num_spec, A);
  r.ana = laa * num_imm;
  r.clado = clado_rate(lac, d, num_spec, K, A);
  return r;
}

inline int pick_mainland(const std::vector<int>& mainland_spec) {
  const int n = static_cast<int>(mainland_spec.size());
  return mainland_spec[sample2_index(n)];
}

// -------------------------------------- constant rates / rate shift ----------

// DAISIE_sim_core_cr (shift_times empty) and DAISIE_sim_core_cr_shift
Rcpp::List sim_const(const Rcpp::List& cfg) {
  const double total_time = Rcpp::as<double>(cfg["total_time"]);
  const std::vector<double> pars = Rcpp::as<std::vector<double>>(cfg["pars"]);
  const double d = Rcpp::as<double>(cfg["d"]);
  const double x = Rcpp::as<double>(cfg["x"]);
  const double mainland_n = Rcpp::as<double>(cfg["mainland_n"]);
  const std::vector<int> mainland_spec =
      Rcpp::as<std::vector<int>>(cfg["mainland_spec"]);
  std::vector<double> shift_times =
      Rcpp::as<std::vector<double>>(cfg["shift_times"]);  // incl. final Inf
  const bool shift = shift_times.size() > 0;

  Island is = island_from_R(Rcpp::as<Rcpp::CharacterMatrix>(cfg["island_spec"]),
                            Rcpp::as<int>(cfg["maxspecID"]), false);
  Rcpp::NumericVector stt0 = cfg["stt_start"];
  is.stt.push_back({stt0[0], stt0[1], stt0[2], stt0[3]});

  double lac, mu, K, gam, laa;
  int rate_set = 1;
  size_t shift_pos = 0;
  auto set_rates = [&](int set) {
    const int o = (set == 1) ? 0 : 5;
    lac = pars[o];
    mu = pars[o + 1];
    K = pars[o + 2];
    gam = pars[o + 3];
    laa = pars[o + 4];
  };
  set_rates(1);
  if (shift) {
    if (shift_times[0] == 0) {
      set_rates(2);
      rate_set = 2;
      shift_pos = 1;
    }
  }
  const double A = 1.0;
  const double extcutoff = kInf;

  int num_spec = is.n();
  int num_imm = is.count_type('I');
  double timeval = 0.0;

  while (timeval < total_time) {
    Rates4 r = rates4(gam, laa, lac, mu, d, x, extcutoff, K, num_spec, num_imm,
                      mainland_n, A);
    // calc_next_timeval(_shift)
    const double totalrate = r.immig + r.ext + r.ana + r.clado;
    bool rate_shift = false;
    if (totalrate != 0.0) {
      timeval = timeval + rexp1(totalrate);
    } else {
      timeval = total_time;
    }
    if (shift && timeval >= shift_times[shift_pos]) {
      timeval = shift_times[shift_pos];
      ++shift_pos;
      rate_shift = true;
    }
    if (shift) {
      if (rate_shift) {
        if (rate_set == 2) {
          set_rates(1);
          rate_set = 1;
        } else {
          set_rates(2);
          rate_set = 2;
        }
        r = rates4(gam, laa, lac, mu, d, x, extcutoff, K, num_spec, num_imm,
                   mainland_n, A);
      }
    } else if (timeval < total_time) {
      r = rates4(gam, laa, lac, mu, d, x, extcutoff, K, num_spec, num_imm,
                 mainland_n, A);
    }
    if (!shift && !(timeval < total_time)) {
      break;
    }
    const int ev = sample_prob({r.immig, r.ext, r.ana, r.clado});

    if (shift ? (timeval <= total_time && !rate_shift) : true) {
      if (ev == 1) {
        is.immigrate(pick_mainland(mainland_spec), timeval, 0);
      } else if (ev == 2) {
        is.extinct(sample2_index(is.n()));
      } else if (ev == 3) {
        std::vector<int> imm;
        for (int i = 0; i < is.n(); ++i) {
          if (is.sp[i].type == 'I') imm.push_back(i);
        }
        int a = imm[0];
        if (imm.size() > 1) a = imm[sample2_index(static_cast<int>(imm.size()))];
        is.anagenesis(a);
      } else {
        is.cladogenesis(sample2_index(is.n()), timeval);
      }
      is.push_stt(total_time - timeval);
      num_spec = is.n();
      num_imm = is.count_type('I');
    }
  }
  return Rcpp::List::create(Rcpp::Named("island_spec") = island_to_R(is),
                            Rcpp::Named("stt_table") = stt_to_R(is));
}

// ------------------------------------------------ time-dependent (thinning) --

Rcpp::List sim_time_dep(const Rcpp::List& cfg) {
  const double total_time = Rcpp::as<double>(cfg["total_time"]);
  const std::vector<double> pars = Rcpp::as<std::vector<double>>(cfg["pars"]);
  const double d = Rcpp::as<double>(cfg["d"]);
  const double x = Rcpp::as<double>(cfg["x"]);
  const double mainland_n = Rcpp::as<double>(cfg["mainland_n"]);
  const std::vector<int> mainland_spec =
      Rcpp::as<std::vector<int>>(cfg["mainland_spec"]);
  const AreaPars ap = area_from_R(Rcpp::as<Rcpp::List>(cfg["area_pars"]));
  const int ontogeny = Rcpp::as<int>(cfg["island_ontogeny"]);
  const int sea_level = Rcpp::as<int>(cfg["sea_level"]);
  const double peak = Rcpp::as<double>(cfg["peak"]);
  const double Amax = Rcpp::as<double>(cfg["Amax"]);
  const double Amin = Rcpp::as<double>(cfg["Amin"]);
  const double extcutoff = Rcpp::as<double>(cfg["extcutoff"]);
  Rcpp::List win = cfg["area_windows"];
  const std::vector<double> w_end = Rcpp::as<std::vector<double>>(win["end"]);
  const std::vector<double> w_Amin = Rcpp::as<std::vector<double>>(win["Amin"]);
  const std::vector<double> w_Amax = Rcpp::as<std::vector<double>>(win["Amax"]);
  const size_t n_win = w_end.size();

  Island is = island_from_R(Rcpp::as<Rcpp::CharacterMatrix>(cfg["island_spec"]),
                            Rcpp::as<int>(cfg["maxspecID"]), false);
  Rcpp::NumericVector stt0 = cfg["stt_start"];
  is.stt.push_back({stt0[0], stt0[1], stt0[2], stt0[3]});

  const double lac = pars[0], mu = pars[1], K = pars[2], gam = pars[3],
               laa = pars[4];
  int num_spec = is.n();
  int num_imm = is.count_type('I');
  double timeval = 0.0;
  size_t w = 0;

  while (timeval < total_time) {
    while (w + 1 < n_win && timeval >= w_end[w]) ++w;
    double win_end, win_Amin, win_Amax;
    if (num_spec == 0) {
      win_end = total_time;
      win_Amin = Amin;
      win_Amax = Amax;
    } else {
      win_end = w_end[w];
      win_Amin = w_Amin[w];
      win_Amax = w_Amax[w];
    }
    // update_max_rates(): list order ext, immig, ana, clado
    const double ext_max = ext_rate(mu, x, extcutoff, num_spec, win_Amin);
    const double immig_max = immig_rate(gam, win_Amax, num_spec, K, mainland_n);
    const double ana_max = laa * num_imm;
    const double clado_max = clado_rate(lac, d, num_spec, K, win_Amax);

    const double totalrate = ext_max + immig_max + ana_max + clado_max;
    double next;
    if (totalrate != 0.0) {
      next = timeval + rexp1(totalrate);
    } else {
      next = total_time;
    }
    if (next >= win_end) {
      timeval = win_end;
      continue;
    }
    timeval = next;

    const double A = island_area(timeval, total_time, ap, peak, ontogeny,
                                 sea_level);
    const Rates4 r = rates4(gam, laa, lac, mu, d, x, extcutoff, K, num_spec,
                            num_imm, mainland_n, A);
    if (r.ext > ext_max * (1 + 1e-8) || r.immig > immig_max * (1 + 1e-8) ||
        r.clado > clado_max * (1 + 1e-8)) {
      Rcpp::stop("Rate exceeds its upper bound in the time-dependent simulation.");
    }
    const int ev = sample_prob({immig_max, ext_max, ana_max, clado_max});
    bool real = false;
    if (ev == 1) {
      if (r.immig / immig_max >= runif1()) {
        real = true;
        is.immigrate(pick_mainland(mainland_spec), timeval, 0);
      }
    } else if (ev == 2) {
      if (r.ext / ext_max >= runif1()) {
        real = true;
        is.extinct(sample2_index(is.n()));
      }
    } else if (ev == 3) {
      real = true;
      std::vector<int> imm;
      for (int i = 0; i < is.n(); ++i) {
        if (is.sp[i].type == 'I') imm.push_back(i);
      }
      int a = imm[0];
      if (imm.size() > 1) a = imm[sample2_index(static_cast<int>(imm.size()))];
      is.anagenesis(a);
    } else {
      if (r.clado / clado_max >= runif1()) {
        real = true;
        is.cladogenesis(sample2_index(is.n()), timeval);
      }
    }
    if (real) {
      is.push_stt(total_time - timeval);
      num_spec = is.n();
      num_imm = is.count_type('I');
    }
  }
  return Rcpp::List::create(Rcpp::Named("island_spec") = island_to_R(is),
                            Rcpp::Named("stt_table") = stt_to_R(is));
}

// ------------------------------------------------------- trait dependent ----

Rcpp::List sim_trait(const Rcpp::List& cfg) {
  const double total_time = Rcpp::as<double>(cfg["total_time"]);
  const std::vector<double> pars = Rcpp::as<std::vector<double>>(cfg["pars"]);
  const double mainland_n = Rcpp::as<double>(cfg["mainland_n"]);
  const std::vector<int> mainland_spec =
      Rcpp::as<std::vector<int>>(cfg["mainland_spec"]);
  Rcpp::List tp = cfg["trait_pars"];
  const double M2 = Rcpp::as<double>(tp["M2"]);
  const double gam2 = Rcpp::as<double>(tp["immig_rate2"]);
  const double mu2 = Rcpp::as<double>(tp["ext_rate2"]);
  const double laa2 = Rcpp::as<double>(tp["ana_rate2"]);
  const double lac2 = Rcpp::as<double>(tp["clado_rate2"]);
  const double trans1 = Rcpp::as<double>(tp["trans_rate"]);
  const double trans2 = Rcpp::as<double>(tp["trans_rate2"]);
  const bool two_K = Rcpp::as<bool>(cfg["two_K"]);
  const double K2 = two_K ? Rcpp::as<double>(tp["K2"]) : 0.0;

  Island is;
  is.maxspecID = Rcpp::as<int>(cfg["maxspecID"]);
  is.trait_mode = true;
  is.push_stt(total_time);
  const double lac = pars[0], mu = pars[1], K = pars[2], gam = pars[3],
               laa = pars[4];
  const int mainland1 = static_cast<int>(mainland_spec.size());
  double timeval = 0.0;

  while (timeval < total_time) {
    const int num_spec = is.n();
    const int n1 = is.count_trait(1);
    const int n2 = is.count_trait(2);
    const int imm1 = is.count_type('I', 1);
    const int imm2 = is.count_type('I', 2);
    double immig1, immig2, clado1, clado2, trans_r1, trans_r2;
    if (two_K) {
      immig1 = std::max(mainland_n * gam * (1.0 - (n1 / K)), 0.0);
      immig2 = std::max(M2 * gam2 * (1.0 - (n2 / K2)), 0.0);
      clado1 = std::max(0.0, lac * n1 * (1.0 - n1 / K));
      clado2 = std::max(0.0, lac2 * n2 * (1.0 - n2 / K2));
      trans_r1 = trans1 * n1 * (1.0 - (n2 / K2));
      // NB: the R code (get_trans_rate) reads trait_pars$K, which does not
      // exist and, through partial matching of $ on lists, resolves to
      // trait_pars$K2. We reproduce this so that R and C++ agree; use K here
      // if get_trans_rate() is changed to use the carrying capacity K.
      trans_r2 = trans2 * n2 * (1.0 - (n1 / K2));
    } else {
      immig1 = std::max(mainland_n * gam * (1.0 - (num_spec / K)), 0.0);
      immig2 = std::max(M2 * gam2 * (1.0 - (num_spec / K)), 0.0);
      clado1 = std::max(0.0, lac * n1 * (1.0 - num_spec / K));
      clado2 = std::max(0.0, lac2 * n2 * (1.0 - num_spec / K));
      trans_r1 = trans1 * n1;
      trans_r2 = trans2 * n2;
    }
    const std::vector<double> r = {immig1,   mu * n1,   laa * imm1, clado1,
                                   trans_r1, immig2,    mu2 * n2,   laa2 * imm2,
                                   clado2,   trans_r2};
    double totalrate = 0.0;
    for (double v : r) totalrate += v;
    if (totalrate != 0.0) {
      timeval = timeval + rexp1(totalrate);
    } else {
      timeval = total_time;
    }
    if (timeval < total_time) {
      const int ev = sample_prob(r);
      auto pick_anag = [&](int state) {
        std::vector<int> imm;
        for (int i = 0; i < is.n(); ++i) {
          if (is.sp[i].type == 'I' && is.sp[i].trait == state) imm.push_back(i);
        }
        int a = imm[0];
        if (imm.size() > 1) a = imm[sample2_index(static_cast<int>(imm.size()))];
        return a;
      };
      if (ev == 1) {
        is.immigrate(pick_mainland(mainland_spec), timeval, 1);
      } else if (ev == 6) {
        const int total = mainland1 + static_cast<int>(M2);
        // sample2((mainland1 + 1):mainland_total, 1)
        const int cnt = total - mainland1;
        const int colonist = mainland1 + 1 + sample2_index(cnt);
        is.immigrate(colonist, timeval, 2);
      } else if (ev == 2 || ev == 7) {
        const std::vector<int> v = is.which_trait(ev == 2 ? 1 : 2);
        is.extinct(v[sample2_index(static_cast<int>(v.size()))]);
      } else if (ev == 3 || ev == 8) {
        is.anagenesis(pick_anag(ev == 3 ? 1 : 2));
      } else if (ev == 4 || ev == 9) {
        const std::vector<int> v = is.which_trait(ev == 4 ? 1 : 2);
        is.cladogenesis(v[sample2_index(static_cast<int>(v.size()))], timeval);
      } else {  // 5 and 10: transition
        const std::vector<int> v = is.which_trait(ev == 5 ? 1 : 2);
        const int t = v[sample2_index(static_cast<int>(v.size()))];
        is.sp[t].trait = (ev == 5) ? 2 : 1;
      }
      is.push_stt(total_time - timeval);
    }
  }
  return Rcpp::List::create(Rcpp::Named("island_spec") = island_to_R(is),
                            Rcpp::Named("stt_table") = stt_to_R(is));
}

// ------------------------------------------------------ SR (Hauffe et al.) --

Rcpp::List sim_sr(const Rcpp::List& cfg) {
  const double total_time = Rcpp::as<double>(cfg["total_time"]);
  const std::vector<double> pars = Rcpp::as<std::vector<double>>(cfg["pars"]);
  const int mainland_n = Rcpp::as<int>(cfg["mainland_n"]);
  std::vector<int> mainland_spec;
  for (int i = 1; i <= mainland_n; ++i) mainland_spec.push_back(i);
  const double shift_t = pars[10];

  Island is;
  is.maxspecID = mainland_n;
  is.trait_mode = false;
  is.stt.push_back({total_time, 0.0, 0.0, 0.0});
  double timeval = 0.0;

  struct R4 {
    double ext, ana, clado, immig;
  };
  auto rates = [&](double lac, double mu, double K, double gam, double laa) {
    R4 r;
    const int n = is.n();
    r.ext = mu * n;
    r.ana = laa * is.count_type('I');
    r.clado = std::max(n * (lac * (1.0 - n / K)), 0.0);
    r.immig = std::max(mainland_n * gam * (1.0 - n / K), 0.0);
    return r;
  };

  while (timeval < total_time) {
    double lac, mu, K, gam, laa;
    if (timeval < shift_t) {
      lac = pars[0]; mu = pars[1]; K = pars[2]; gam = pars[3]; laa = pars[4];
    } else {
      lac = pars[5]; mu = pars[6]; K = pars[7]; gam = pars[8]; laa = pars[9];
    }
    R4 r = rates(lac, mu, K, gam, laa);
    double totalrate = r.ext + r.clado + r.ana + r.immig;
    double dt = rexp1(totalrate);
    if (std::isnan(dt)) dt = total_time - timeval;
    if (timeval < shift_t && (timeval + dt) >= shift_t) {
      lac = pars[5]; mu = pars[6]; K = pars[7]; gam = pars[8]; laa = pars[9];
      r = rates(lac, mu, K, gam, laa);
      totalrate = r.ext + r.clado + r.ana + r.immig;
      dt = rexp1(totalrate);
      if (std::isnan(dt)) dt = total_time - shift_t;
      timeval = shift_t + dt;
    } else {
      timeval = timeval + dt;
    }
    const int ev = sample_prob({r.immig, r.ext, r.ana, r.clado});
    if (timeval < total_time) {
      if (ev == 1) {
        is.immigrate(pick_mainland(mainland_spec), timeval, 0);
      } else if (ev == 2) {
        is.extinct(sample2_index(is.n()));
      } else if (ev == 3) {
        std::vector<int> imm;
        for (int i = 0; i < is.n(); ++i) {
          if (is.sp[i].type == 'I') imm.push_back(i);
        }
        int a = imm[0];
        if (imm.size() > 1) a = imm[sample2_index(static_cast<int>(imm.size()))];
        is.anagenesis(a);
      } else {
        is.cladogenesis(sample2_index(is.n()), timeval);
      }
    }
    // DAISIE_SR_sim_core adds a row after every iteration
    is.push_stt(total_time - timeval);
  }
  is.stt.back()[0] = 0.0;
  return Rcpp::List::create(Rcpp::Named("island_spec") = island_to_R(is),
                            Rcpp::Named("stt_table") = stt_to_R(is));
}

}  // namespace

// Entry point: cfg$mode selects the simulation core.
extern "C" SEXP daisie_sim_engine(SEXP cfgSEXP) {
  BEGIN_RCPP
  Rcpp::RNGScope rngScope;
  Rcpp::List cfg(cfgSEXP);
  const std::string mode = Rcpp::as<std::string>(cfg["mode"]);
  if (mode == "const") return sim_const(cfg);
  if (mode == "time_dep") return sim_time_dep(cfg);
  if (mode == "trait") return sim_trait(cfg);
  if (mode == "sr") return sim_sr(cfg);
  Rcpp::stop("unknown simulation mode");
  END_RCPP
}
