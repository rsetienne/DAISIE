// Core of the archipelago simulation DAISIE_A (version 10, divergence states).
//
// This is a line-by-line translation of the R reference implementation in
// R/DAISIE_A_sim_core.R (DAISIE_A_sim_core_R() and the da_* helpers): the same
// events in the same order, the same arithmetic and the same calls to R's
// random number generator (R::rexp, R::runif), so that both implementations
// give identical output for the same seed. Any change here must be made there
// too.
//
// Indices are 0-based here: islands 0..n-1, the mainland is node n. Genealogy
// node IDs are 1-based, as in R.

#include <Rcpp.h>
#include <R_ext/Random.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace {

const double NEG_INF = -std::numeric_limits<double>::infinity();
const double POS_INF = std::numeric_limits<double>::infinity();
const double BIG = 2147483647.0;  // .Machine$integer.max

struct Lineage {
  double id, anc, colt, type, ana, node, side, clade;
  std::vector<int> occ;
  std::vector<double> S;    // (n + 1)^2, -1 = not tracked
  std::vector<double> Tau;  // (n + 1)^2
  std::vector<int> s;
  std::vector<double> stau;
};

struct Cfg {
  int n = 0, n1 = 0, k = 1, mainland_n = 1;
  std::vector<double> lac, mu, K, gam, beta, D, emerge, submerge;
  std::vector<int> A1;
  double alpha = 0.0, total_time = 0.0;
  bool split_single = false, protracted = false;
  bool adj(int a, int b) const { return A1[a * n1 + b] != 0; }
  double d(int i, int j) const { return D[i * n + j]; }
};

struct Event {
  int type, r, a, b;
};

struct Sim {
  const Cfg& c;
  std::vector<Lineage> lin;
  std::vector<std::array<double, 8>> gen;
  int n_clades = 0;
  double maxspecID = 0.0;
  // species-through-time tables
  std::vector<std::array<double, 4>> stt_table;
  std::vector<std::vector<double>> stt_islands;
  double n_events = 0.0;

  explicit Sim(const Cfg& cfg) : c(cfg) {}

  // ---- lineages and pairs -------------------------------------------------
  Lineage new_lineage(double id, double anc, double colt, double type,
                      double ana, double node, double side, double clade) const {
    Lineage L;
    L.id = id; L.anc = anc; L.colt = colt; L.type = type; L.ana = ana;
    L.node = node; L.side = side; L.clade = clade;
    L.occ.assign(c.n, 0);
    L.S.assign(c.n1 * c.n1, -1.0);
    L.Tau.assign(c.n1 * c.n1, NEG_INF);
    L.s.assign(c.n, 0);
    L.stau.assign(c.n, NEG_INF);
    return L;
  }
  double& S(Lineage& L, int a, int b) const { return L.S[a * c.n1 + b]; }
  double Sv(const Lineage& L, int a, int b) const { return L.S[a * c.n1 + b]; }
  double Tv(const Lineage& L, int a, int b) const { return L.Tau[a * c.n1 + b]; }
  void set_pair(Lineage& L, int a, int b, double s, double tau) const {
    L.S[a * c.n1 + b] = s;
    L.S[b * c.n1 + a] = s;
    L.Tau[a * c.n1 + b] = tau;
    L.Tau[b * c.n1 + a] = tau;
  }
  int range_size(const Lineage& L) const {
    int m = 0;
    for (int i = 0; i < c.n; ++i) m += L.occ[i];
    return m;
  }
  std::vector<int> nodes(const Lineage& L) const {
    std::vector<int> v;
    for (int i = 0; i < c.n; ++i) if (L.occ[i] == 1) v.push_back(i);
    if (L.type == 1.0) v.push_back(c.n);
    return v;
  }
  // breadth-first components; edge(u, w) decides adjacency
  template <typename F>
  std::vector<int> comp(const std::vector<int>& nd, F edge) const {
    std::vector<int> cp(c.n1, 0);
    int nc = 0;
    for (int v : nd) {
      if (cp[v] != 0) continue;
      ++nc;
      cp[v] = nc;
      std::vector<int> queue(1, v);
      for (size_t h = 0; h < queue.size(); ++h) {
        int u = queue[h];
        for (int w : nd) {
          if (cp[w] == 0 && edge(u, w)) {
            cp[w] = nc;
            queue.push_back(w);
          }
        }
      }
    }
    return cp;
  }
  static bool less(double s1, double t1, double s2, double t2) {
    return s1 < s2 || (s1 == s2 && t1 < t2);
  }
  void bottleneck(const Lineage& L, const std::vector<int>& nd,
                  std::vector<double>& bs, std::vector<double>& bt) const {
    const int n1 = c.n1;
    bs.assign(n1 * n1, BIG);
    bt.assign(n1 * n1, POS_INF);
    for (int a : nd) for (int b : nd) {
      if (a == b) {
        bs[a * n1 + b] = 0.0;
        bt[a * n1 + b] = NEG_INF;
      } else if (Sv(L, a, b) >= 0) {
        bs[a * n1 + b] = Sv(L, a, b);
        bt[a * n1 + b] = Tv(L, a, b);
      }
    }
    for (int m : nd) for (int a : nd) {
      if (bs[a * n1 + m] == BIG) continue;
      for (int b : nd) {
        if (bs[m * n1 + b] == BIG) continue;
        double cs, ct;
        if (less(bs[a * n1 + m], bt[a * n1 + m], bs[m * n1 + b], bt[m * n1 + b])) {
          cs = bs[m * n1 + b];
          ct = bt[m * n1 + b];
        } else {
          cs = bs[a * n1 + m];
          ct = bt[a * n1 + m];
        }
        if (less(cs, ct, bs[a * n1 + b], bt[a * n1 + b])) {
          bs[a * n1 + b] = cs;
          bt[a * n1 + b] = ct;
        }
      }
    }
  }
  void drop_gaps(Lineage& L) const {
    std::vector<int> nd = nodes(L);
    std::vector<int> cp = comp(nd, [&](int u, int w) {
      return c.adj(u, w) && Sv(L, u, w) >= 0 && Sv(L, u, w) < c.k;
    });
    for (int a : nd) for (int b : nd) {
      if (a < b && Sv(L, a, b) >= 0 && !c.adj(a, b) && cp[a] == cp[b]) {
        set_pair(L, a, b, -1.0, NEG_INF);
      }
    }
  }
  void add_island(Lineage& L, int j, int src) const {
    std::vector<int> nd = nodes(L);
    std::vector<double> bs, bt;
    bottleneck(L, nd, bs, bt);
    for (int l : nd) {
      if (!c.adj(j, l)) continue;
      if (l == src) {
        set_pair(L, j, l, 0.0, NEG_INF);
      } else if (Sv(L, src, l) >= 0) {
        set_pair(L, j, l, Sv(L, src, l), Tv(L, src, l));
      } else {
        set_pair(L, j, l, bs[src * c.n1 + l], bt[src * c.n1 + l]);
      }
    }
    L.occ[j] = 1;
    L.s[j] = 0;
    L.stau[j] = NEG_INF;
    drop_gaps(L);
  }
  void remove_island(Lineage& L, int l) const {
    std::vector<int> nd = nodes(L);
    std::vector<int> nb;
    for (int v : nd) if (v != l && Sv(L, l, v) >= 0) nb.push_back(v);
    std::vector<double> ls(c.n1), lt(c.n1);
    for (int v = 0; v < c.n1; ++v) {
      ls[v] = Sv(L, l, v);
      lt[v] = Tv(L, l, v);
    }
    for (int v = 0; v < c.n1; ++v) set_pair(L, l, v, -1.0, NEG_INF);
    L.occ[l] = 0;
    L.s[l] = 0;
    L.stau[l] = NEG_INF;
    std::vector<int> nd2;
    for (int v : nd) if (v != l) nd2.push_back(v);
    std::vector<int> cp = comp(nd2, [&](int u, int w) { return Sv(L, u, w) >= 0; });
    const int m = static_cast<int>(nb.size());
    for (int x = 0; x + 1 < m; ++x) for (int y = x + 1; y < m; ++y) {
      int a = nb[x], b = nb[y];
      if (cp[a] == cp[b]) continue;
      if (less(ls[a], lt[a], ls[b], lt[b])) {
        set_pair(L, a, b, ls[b], lt[b]);
      } else {
        set_pair(L, a, b, ls[a], lt[a]);
      }
    }
  }
  Lineage restrict(const Lineage& L0, const std::vector<int>& P) const {
    Lineage L = L0;
    std::vector<int> keep(c.n1, 0);
    for (int v : P) keep[v] = 1;
    for (int a = 0; a < c.n1; ++a) for (int b = 0; b < c.n1; ++b) {
      if (!keep[a] || !keep[b]) {
        L.S[a * c.n1 + b] = -1.0;
        L.Tau[a * c.n1 + b] = NEG_INF;
      }
    }
    for (int i = 0; i < c.n; ++i) {
      if (!keep[i]) {
        L.occ[i] = 0;
        L.s[i] = 0;
        L.stau[i] = NEG_INF;
      }
    }
    return L;
  }
  double cut_tau(const Lineage& L, const std::vector<int>& P,
                 const std::vector<int>& Q) const {
    double m = NEG_INF;
    for (int a : P) for (int b : Q) {
      if (Sv(L, a, b) >= 0 && Tv(L, a, b) > m) m = Tv(L, a, b);
    }
    return m;
  }
  std::vector<int> perm(int m) const {
    std::vector<int> x(m);
    for (int i = 0; i < m; ++i) x[i] = i;
    for (int j = m; j >= 2; --j) {
      int idx = static_cast<int>(std::floor(R::runif(0.0, 1.0) * j));
      std::swap(x[j - 1], x[idx]);
    }
    return x;
  }

  // ---- genealogy ----------------------------------------------------------
  int add_node(double parent, double time, double kind, double side,
               double clade, double species, double time_split) {
    std::array<double, 8> row = {static_cast<double>(gen.size() + 1), parent,
                                 time, kind, side, clade, species, time_split};
    gen.push_back(row);
    return static_cast<int>(gen.size());
  }
  double gen_time(double node) const {
    return gen[static_cast<size_t>(node) - 1][2];
  }
  void clade_singleton(double clade) {
    int count = 0, row = -1;
    for (size_t q = 0; q < lin.size(); ++q) {
      if (lin[q].clade == clade) { ++count; row = static_cast<int>(q); }
    }
    if (count == 1 && lin[row].type == 3.0) {
      lin[row].type = 2.0;
      lin[row].ana = 2.0;
    }
  }
  void set_clade_colt(double clade, double tau) {
    for (auto& L : lin) if (L.clade == clade) L.colt = tau;
    for (auto& g : gen) if (g[3] == 1.0 && g[5] == clade) g[2] = tau;
  }

  // ---- splits -------------------------------------------------------------
  static std::vector<int> concat(const std::vector<std::vector<int>>& ord,
                                 size_t from, const std::vector<int>& extra) {
    std::vector<int> out;
    for (size_t q = from; q < ord.size(); ++q) {
      out.insert(out.end(), ord[q].begin(), ord[q].end());
    }
    out.insert(out.end(), extra.begin(), extra.end());
    return out;
  }
  void peel(int r, const std::vector<std::vector<int>>& ord, int n_peel,
            const std::vector<int>& extra, double t) {
    for (int j = 0; j < n_peel; ++j) {
      const std::vector<int>& P = ord[j];
      std::vector<int> Q = concat(ord, j + 1, extra);
      Lineage L = lin[r];
      double t_node = c.protracted ? std::max(cut_tau(L, P, Q), gen_time(L.node)) : t;
      int v = add_node(L.node, t_node, 2.0, L.side, L.clade, NA_REAL, t);
      Lineage B = restrict(L, P);
      maxspecID += 1.0;
      B.id = maxspecID;
      B.type = 3.0;
      B.ana = 0.0;
      B.node = v;
      B.side = 2.0;
      L = restrict(L, Q);
      L.node = v;
      L.side = 1.0;
      lin[r] = L;
      lin.push_back(B);
    }
  }
  void split(int r, std::vector<std::vector<int>> pieces, double t) {
    const int M = c.n;
    Lineage L = lin[r];
    if (L.type == 1.0) {
      size_t pm = 0;
      for (size_t q = 0; q < pieces.size(); ++q) {
        if (std::find(pieces[q].begin(), pieces[q].end(), M) != pieces[q].end()) pm = q;
      }
      std::vector<std::vector<int>> others;
      for (size_t q = 0; q < pieces.size(); ++q) if (q != pm) others.push_back(pieces[q]);
      if (pieces[pm].size() > 1) {
        std::vector<std::vector<int>> ord = others;
        if (others.size() >= 2) {
          std::vector<int> p = perm(static_cast<int>(others.size()));
          for (size_t q = 0; q < p.size(); ++q) ord[q] = others[p[q]];
        }
        peel(r, ord, static_cast<int>(ord.size()), pieces[pm], t);
        return;
      }
      std::vector<int> isl = concat(others, 0, std::vector<int>());
      double tau_m = cut_tau(L, std::vector<int>(1, M), isl);
      lin[r] = restrict(L, isl);
      if (c.protracted) set_clade_colt(L.clade, tau_m);
      if (others.size() == 1) {
        maxspecID += 1.0;
        lin[r].id = maxspecID;
        lin[r].type = 2.0;
        lin[r].ana = 1.0;
        return;
      }
      pieces = others;
    }
    const int m = static_cast<int>(pieces.size());
    std::vector<std::vector<int>> ord = pieces;
    if (m >= 3) {
      std::vector<int> p = perm(m);
      for (int q = 0; q < m; ++q) ord[q] = pieces[p[q]];
    }
    peel(r, ord, m - 1, std::vector<int>(), t);
    maxspecID += 1.0;
    lin[r].id = maxspecID;
    lin[r].type = 3.0;
    lin[r].ana = 0.0;
  }
  void check_split(int r, double t) {
    const Lineage& L = lin[r];
    std::vector<int> nd = nodes(L);
    std::vector<int> cp = comp(nd, [&](int u, int w) {
      return Sv(L, u, w) >= 0 && Sv(L, u, w) < c.k;
    });
    int nc = 0;
    for (int v : nd) nc = std::max(nc, cp[v]);
    if (nc == 1) return;
    std::vector<std::vector<int>> pieces(nc);
    for (int x = 1; x <= nc; ++x) {
      for (int v : nd) if (cp[v] == x) pieces[x - 1].push_back(v);
    }
    split(r, pieces, t);
  }

  // ---- events ---------------------------------------------------------------
  std::vector<double> rates;
  std::vector<Event> events;
  void add_event(double x, int type, int r, int a, int b) {
    rates.push_back(x);
    events.push_back(Event{type, r, a, b});
  }
  void list_events(double timeval) {
    rates.clear();
    events.clear();
    const int n = c.n, k = c.k;
    std::vector<int> exists(n);
    std::vector<double> f(n);
    for (int i = 0; i < n; ++i) {
      exists[i] = c.emerge[i] <= timeval && c.submerge[i] > timeval;
      double ni = 0.0;
      for (const auto& L : lin) ni += L.occ[i];
      f[i] = std::max(0.0, 1.0 - ni / c.K[i]);
    }
    std::vector<int> ne_row(c.mainland_n, -1);
    for (size_t r = 0; r < lin.size(); ++r) {
      if (lin[r].type == 1.0) ne_row[static_cast<int>(lin[r].id) - 1] = static_cast<int>(r);
    }
    for (int m = 0; m < c.mainland_n; ++m) {
      for (int i = 0; i < n; ++i) {
        if (!exists[i] || c.gam[i] == 0.0) continue;
        int r = ne_row[m];
        double x;
        if (r < 0) {
          x = c.gam[i] * f[i];
        } else {
          const Lineage& L = lin[r];
          if (L.occ[i] == 0) {
            x = c.gam[i] * f[i];
          } else if (range_size(L) == 1) {
            x = c.gam[i];
          } else if (Sv(L, i, n) > 0 || (c.protracted && L.s[i] > 0)) {
            x = c.gam[i];
          } else {
            x = 0.0;
          }
        }
        if (x > 0) add_event(x, 1, m, i, -1);
      }
    }
    for (size_t rr = 0; rr < lin.size(); ++rr) {
      const int r = static_cast<int>(rr);
      const Lineage& L = lin[rr];
      std::vector<int> R;
      for (int i = 0; i < n; ++i) if (L.occ[i] == 1) R.push_back(i);
      for (int i : R) if (c.mu[i] > 0) add_event(c.mu[i], 2, r, i, -1);
      for (int i : R) {
        double x;
        if (c.protracted) {
          x = (L.s[i] == k - 1) ? c.beta[i] * f[i] : c.beta[i];
        } else {
          x = c.lac[i] * f[i];
        }
        if (x > 0) add_event(x, 3, r, i, -1);
      }
      for (int i : R) {
        for (int j = 0; j < n; ++j) {
          if (j == i || c.d(i, j) == 0.0 || !exists[j]) continue;
          double x;
          if (L.occ[j] == 0) {
            x = c.d(i, j) * f[j];
          } else if (Sv(L, i, j) > 0 || (c.protracted && L.s[j] > 0)) {
            x = c.d(i, j);
          } else {
            x = 0.0;
          }
          if (x > 0) add_event(x, 4, r, i, j);
        }
      }
      if (c.alpha > 0) {
        std::vector<int> nd = nodes(L);
        for (int a : nd) for (int b : nd) {
          if (a < b && Sv(L, a, b) >= 0 && Sv(L, a, b) < k) add_event(c.alpha, 5, r, a, b);
        }
      }
    }
  }

  void two_new_species(const Lineage& L, int i, double v, Lineage& A, Lineage& B) {
    for (int side = 1; side <= 2; ++side) {
      maxspecID += 1.0;
      Lineage x = new_lineage(maxspecID, L.anc, L.colt, 3.0, 0.0, v, side, L.clade);
      x.occ[i] = 1;
      if (side == 1) A = x; else B = x;
    }
  }
  void within_speciation(int r, int i, double t) {
    Lineage L = lin[r];
    int v = add_node(L.node, t, 2.0, L.side, L.clade, NA_REAL, t);
    if (range_size(L) == 1 && c.split_single) {
      Lineage A, B;
      two_new_species(L, i, v, A, B);
      lin[r] = A;
      lin.push_back(B);
    } else {
      L.node = v;
      L.side = 1.0;
      lin[r] = L;
      maxspecID += 1.0;
      Lineage B = new_lineage(maxspecID, L.anc, L.colt, 3.0, 0.0, v, 2.0, L.clade);
      B.occ[i] = 1;
      lin.push_back(B);
    }
  }
  void within_step(int r, int i, double t) {
    Lineage L = lin[r];
    L.s[i] += 1;
    if (L.s[i] == 1) L.stau[i] = t;
    if (L.s[i] < c.k) {
      lin[r] = L;
      return;
    }
    if (range_size(L) == 1) {
      int v = add_node(L.node, t, 2.0, L.side, L.clade, NA_REAL, t);
      Lineage A, B;
      two_new_species(L, i, v, A, B);
      lin[r] = A;
      lin.push_back(B);
    } else {
      double t_start = std::max(L.stau[i], gen_time(L.node));
      int v1 = add_node(L.node, t_start, 2.0, L.side, L.clade, NA_REAL, t);
      int v2 = add_node(v1, t, 2.0, 2.0, L.clade, NA_REAL, t);
      Lineage A, B;
      two_new_species(L, i, v2, A, B);
      L.node = v1;
      L.side = 1.0;
      remove_island(L, i);
      lin[r] = L;
      lin.push_back(A);
      lin.push_back(B);
      check_split(r, t);
    }
  }
  void apply(const Event& ev, double t) {
    const int n = c.n;
    if (ev.type == 1) {
      const double m = ev.r + 1;  // mainland species ID
      const int i = ev.a;
      int r = -1;
      for (size_t q = 0; q < lin.size(); ++q) {
        if (lin[q].type == 1.0 && lin[q].id == m) r = static_cast<int>(q);
      }
      if (r < 0 || (range_size(lin[r]) == 1 && lin[r].occ[i] == 1)) {
        double old_clade = 0.0;
        if (r >= 0) {
          const Lineage& old = lin[r];
          old_clade = old.clade;
          add_node(old.node, t, 4.0, old.side, old.clade, old.id, t);
        }
        n_clades += 1;
        int v = add_node(0.0, t, 1.0, 0.0, n_clades, m, t);
        Lineage L = new_lineage(m, m, t, 1.0, 0.0, v, 0.0, n_clades);
        L.occ[i] = 1;
        S(L, i, n) = 0.0;
        S(L, n, i) = 0.0;
        if (r < 0) {
          lin.push_back(L);
        } else {
          lin[r] = L;
          clade_singleton(old_clade);
        }
      } else {
        Lineage L = lin[r];
        if (L.occ[i] == 0) {
          add_island(L, i, n);
        } else {
          set_pair(L, i, n, 0.0, NEG_INF);
          if (c.protracted) {
            L.s[i] = 0;
            L.stau[i] = NEG_INF;
          }
          drop_gaps(L);
        }
        lin[r] = L;
      }
    } else if (ev.type == 2) {
      Lineage L = lin[ev.r];
      if (range_size(L) == 1) {
        add_node(L.node, t, 3.0, L.side, L.clade, L.id, t);
        lin.erase(lin.begin() + ev.r);
        clade_singleton(L.clade);
      } else {
        remove_island(L, ev.a);
        lin[ev.r] = L;
        check_split(ev.r, t);
      }
    } else if (ev.type == 3) {
      if (c.protracted) {
        within_step(ev.r, ev.a, t);
      } else {
        within_speciation(ev.r, ev.a, t);
      }
    } else if (ev.type == 4) {
      Lineage L = lin[ev.r];
      if (L.occ[ev.b] == 0) {
        add_island(L, ev.b, ev.a);
      } else {
        set_pair(L, ev.a, ev.b, 0.0, NEG_INF);
        if (c.protracted) {
          L.s[ev.b] = 0;
          L.stau[ev.b] = NEG_INF;
        }
        drop_gaps(L);
      }
      lin[ev.r] = L;
    } else {
      Lineage L = lin[ev.r];
      double new_s = Sv(L, ev.a, ev.b) + 1.0;
      set_pair(L, ev.a, ev.b, new_s, new_s == 1.0 ? t : Tv(L, ev.a, ev.b));
      lin[ev.r] = L;
      if (new_s == c.k) check_split(ev.r, t);
    }
  }
  void submerge(int i, double t) {
    size_t r = 0;
    while (r < lin.size()) {
      if (lin[r].occ[i] == 0) {
        ++r;
        continue;
      }
      Lineage L = lin[r];
      if (range_size(L) == 1) {
        add_node(L.node, t, 3.0, L.side, L.clade, L.id, t);
        lin.erase(lin.begin() + r);
        clade_singleton(L.clade);
        continue;
      }
      remove_island(L, i);
      lin[r] = L;
      check_split(static_cast<int>(r), t);
      ++r;
    }
  }

  // ---- species-through-time ---------------------------------------------
  void init_stt() {
    stt_table.push_back({c.total_time, 0.0, 0.0, 0.0});
    std::vector<double> row(1 + 4 * c.n, 0.0);
    row[0] = c.total_time;
    stt_islands.push_back(row);
  }
  void stt_row(double time_bp) {
    std::array<double, 4> arch = {time_bp, 0.0, 0.0, 0.0};
    std::vector<double> isl(1 + 4 * c.n, 0.0);
    isl[0] = time_bp;
    for (const auto& L : lin) {
      const int ty = static_cast<int>(L.type);
      arch[ty] += 1.0;
      const bool single = L.type != 1.0 && range_size(L) == 1;
      for (int i = 0; i < c.n; ++i) {
        if (L.occ[i] != 1) continue;
        isl[1 + 4 * i + ty - 1] += 1.0;
        if (single) isl[1 + 4 * i + 3] += 1.0;
      }
    }
    const auto& last_a = stt_table.back();
    const auto& last_i = stt_islands.back();
    bool same = true;
    for (int q = 1; q < 4; ++q) if (last_a[q] != arch[q]) same = false;
    for (size_t q = 1; q < isl.size(); ++q) if (last_i[q] != isl[q]) same = false;
    if (same) return;
    stt_table.push_back(arch);
    stt_islands.push_back(isl);
  }

  // ---- main loop ------------------------------------------------------------
  void run() {
    const double total_time = c.total_time;
    maxspecID = c.mainland_n;
    init_stt();
    double timeval = 0.0;
    std::vector<double> cs;
    while (timeval < total_time) {
      Rcpp::checkUserInterrupt();
      list_events(timeval);
      cs.assign(rates.size(), 0.0);
      double acc = 0.0;
      for (size_t e = 0; e < rates.size(); ++e) {
        acc = acc + rates[e];
        cs[e] = acc;
      }
      const double total = rates.empty() ? 0.0 : cs.back();
      double new_time;
      if (total > 0) {
        new_time = timeval + R::rexp(1.0 / total);
      } else {
        new_time = total_time;
      }
      double next_island_event = POS_INF;
      for (int i = 0; i < c.n; ++i) {
        if (c.emerge[i] > timeval) next_island_event = std::min(next_island_event, c.emerge[i]);
        if (c.submerge[i] > timeval) next_island_event = std::min(next_island_event, c.submerge[i]);
      }
      if (next_island_event < total_time && new_time >= next_island_event) {
        timeval = next_island_event;
        bool any = false;
        for (int i = 0; i < c.n; ++i) {
          if (c.submerge[i] == timeval) {
            submerge(i, timeval);
            any = true;
          }
        }
        if (any) stt_row(total_time - timeval);
        continue;
      }
      timeval = new_time;
      if (timeval < total_time) {
        const double target = R::runif(0.0, 1.0) * total;
        size_t e = 0;
        while (cs[e] <= target) ++e;
        const Event ev = events[e];
        apply(ev, timeval);
        n_events += 1.0;
        stt_row(total_time - timeval);
      }
    }
  }
};

std::vector<double> num_vec(const Rcpp::List& x, const char* name) {
  return Rcpp::as<std::vector<double>>(x[name]);
}

}  // namespace

extern "C" SEXP DAISIE_A_sim_cpp(SEXP cfgSEXP) {
  BEGIN_RCPP
  Rcpp::RNGScope rngScope;
  Rcpp::List x(cfgSEXP);
  Cfg c;
  c.n = Rcpp::as<int>(x["n"]);
  c.n1 = c.n + 1;
  c.k = Rcpp::as<int>(x["k"]);
  c.mainland_n = Rcpp::as<int>(x["mainland_n"]);
  c.lac = num_vec(x, "lac");
  c.mu = num_vec(x, "mu");
  c.K = num_vec(x, "K");
  c.gam = num_vec(x, "gam");
  c.beta = num_vec(x, "beta");
  c.emerge = num_vec(x, "emerge");
  c.submerge = num_vec(x, "submerge");
  c.alpha = Rcpp::as<double>(x["alpha"]);
  c.total_time = Rcpp::as<double>(x["total_time"]);
  c.split_single = Rcpp::as<bool>(x["split_single"]);
  c.protracted = Rcpp::as<bool>(x["protracted"]);
  Rcpp::NumericMatrix D = x["D"];
  Rcpp::LogicalMatrix A1 = x["A1"];
  c.D.assign(c.n * c.n, 0.0);
  for (int i = 0; i < c.n; ++i) for (int j = 0; j < c.n; ++j) c.D[i * c.n + j] = D(i, j);
  c.A1.assign(c.n1 * c.n1, 0);
  for (int i = 0; i < c.n1; ++i) for (int j = 0; j < c.n1; ++j) c.A1[i * c.n1 + j] = A1(i, j) ? 1 : 0;

  Sim sim(c);
  sim.run();

  const int L = static_cast<int>(sim.lin.size());
  Rcpp::NumericMatrix spec(L, 8 + c.n);
  for (int r = 0; r < L; ++r) {
    const Lineage& l = sim.lin[r];
    const double fixed[8] = {l.id, l.anc, l.colt, l.type, l.ana, l.node, l.side, l.clade};
    for (int q = 0; q < 8; ++q) spec(r, q) = fixed[q];
    for (int i = 0; i < c.n; ++i) spec(r, 8 + i) = l.occ[i];
  }
  const int G = static_cast<int>(sim.gen.size());
  Rcpp::NumericMatrix gen(G, 8);
  for (int g = 0; g < G; ++g) for (int q = 0; q < 8; ++q) gen(g, q) = sim.gen[g][q];
  // append the final rows (time 0) of the species-through-time tables
  sim.stt_table.push_back(sim.stt_table.back());
  sim.stt_table.back()[0] = 0.0;
  sim.stt_islands.push_back(sim.stt_islands.back());
  sim.stt_islands.back()[0] = 0.0;
  const int T = static_cast<int>(sim.stt_table.size());
  Rcpp::NumericMatrix stt_table(T, 4);
  Rcpp::NumericMatrix stt_islands(T, 1 + 4 * c.n);
  for (int q = 0; q < T; ++q) {
    for (int w = 0; w < 4; ++w) stt_table(q, w) = sim.stt_table[q][w];
    for (int w = 0; w < 1 + 4 * c.n; ++w) stt_islands(q, w) = sim.stt_islands[q][w];
  }
  return Rcpp::List::create(Rcpp::Named("spec") = spec, Rcpp::Named("gen") = gen,
                            Rcpp::Named("stt_table") = stt_table,
                            Rcpp::Named("stt_islands") = stt_islands,
                            Rcpp::Named("n_events") = sim.n_events);
  END_RCPP
}
