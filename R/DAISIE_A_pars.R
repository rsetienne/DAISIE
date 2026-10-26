#' Create the parameter object of the archipelago DAISIE model (DAISIE_A)
#'
#' @description
#' Collects, expands and checks all rates and settings of the archipelago
#' extension of DAISIE (DAISIE_A). In DAISIE_A every lineage on the
#' archipelago has a range: the non-empty set of islands on which it occurs.
#' For \eqn{n} islands there are \eqn{2^n - 1} possible ranges. Lineages
#' expand their range by dispersal, contract it by local extinction (losing the
#' last island is global extinction), and speciate either within an island
#' (in-situ cladogenesis; budding for widespread lineages) or between islands
#' (vicariance). Anagenesis acts either on the whole lineage or on the
#' population of a single island of a widespread lineage (island-level
#' anagenesis). Colonisation from the mainland is always to a single island.
#'
#' Founder-event (jump) speciation is not a separate process: a founder event
#' is dispersal to a new island followed by island-level anagenesis on that
#' island.
#'
#' All island-specific rates (\code{lac}, \code{mu}, \code{K}, \code{gam},
#' \code{laa_isl}) can be given as a single number, which is then used for every island, or as
#' a vector with one value per island.
#'
#' @param n_islands Integer, number of islands in the archipelago.
#' @param lac Within-island (in-situ) cladogenesis rate per island population
#'   (\eqn{\lambda^c_i}). Scalar or vector of length \code{n_islands}.
#' @param mu Local extinction rate per island population (\eqn{\mu_i}).
#'   For a lineage on a single island this is the global extinction rate.
#'   Scalar or vector of length \code{n_islands}.
#' @param K Clade-specific carrying capacity per island (\eqn{K_i}). Arrivals
#'   on island \eqn{i} (colonisation, dispersal) and in-situ
#'   speciation on island \eqn{i} are multiplied by
#'   \eqn{\max(0, 1 - n_i / K_i)}, where \eqn{n_i} is the number of species of
#'   the clade on island \eqn{i}. Use \code{Inf} for no diversity dependence.
#'   Scalar or vector of length \code{n_islands}.
#' @param gam Colonisation rate from the mainland to island \eqn{i}, per
#'   mainland species (\eqn{\gamma_i}). Scalar or vector of length
#'   \code{n_islands}. Note that a scalar is the rate per island, so the total
#'   colonisation rate of an empty archipelago is \code{n_islands * gam}.
#' @param laa Anagenesis rate per non-endemic lineage (\eqn{\lambda^a}). The
#'   whole (possibly multi-island) lineage becomes endemic, for example
#'   because the mainland sister lineage changes or speciates.
#' @param laa_isl Island-level anagenesis rate per island population
#'   (\eqn{\lambda^{ia}_i}; default 0). The population on island \eqn{i} of a
#'   lineage (non-endemic or endemic) with range \eqn{R}, \eqn{|R| \ge 2},
#'   becomes a new endemic species with range \eqn{\{i\}}; the parent keeps
#'   its identity and type with range \eqn{R \setminus \{i\}}. For a lineage
#'   on a single island this would be whole-lineage anagenesis, so it does not
#'   apply there. Scalar or vector of length \code{n_islands}.
#' @param lav Between-island (vicariant) speciation rate
#'   (\eqn{\lambda^v}) of lineages that occur on two or more islands, either
#'   per lineage or per admissible split, see \code{vic_per}.
#' @param disp Dispersal (range expansion) rate multiplier. The rate at which a
#'   lineage with range \eqn{R} expands to island \eqn{j} is
#'   \code{disp} \eqn{\times \sum_{i \in R} C_{ij}}, with \eqn{C} the
#'   connectivity matrix.
#' @param connectivity Either a character string (\code{"full"}, \code{"line"},
#'   \code{"ring"} or \code{"none"}), passed to
#'   \code{\link{DAISIE_A_connectivity}()}, or a non-negative
#'   \code{n_islands} by \code{n_islands} matrix \eqn{C} of relative
#'   dispersal weights from row island to column island (diagonal ignored).
#'   The non-zero pattern of \eqn{C + C^T} also defines which islands are
#'   adjacent (used by \code{vic_mode = "connected"} and
#'   \code{connected_ranges = TRUE}).
#' @param vic_mode How the range of a lineage is split in a vicariance event.
#'   All admissible splits are equally likely:
#'   \itemize{
#'     \item \code{"connected"} (default): if the range forms one connected
#'       block in the adjacency graph, both daughter ranges must be connected;
#'       if the range consists of several disconnected blocks, the split
#'       follows the gaps (blocks are never broken up).
#'     \item \code{"all"}: all \eqn{2^{|R|-1} - 1} bipartitions.
#'     \item \code{"single"}: one island versus the rest.
#'   }
#' @param vic_per To what the vicariance rate \code{lav} applies:
#'   \itemize{
#'     \item \code{"lineage"} (default): every lineage with two or more
#'       islands splits at rate \code{lav}, whatever its range; the split is
#'       chosen uniformly among the admissible splits.
#'     \item \code{"split"}: every admissible split of the range occurs at rate
#'       \code{lav}, so a lineage with range \eqn{R} splits at rate
#'       \code{lav} \eqn{\times s(R)}, with \eqn{s(R)} the number of
#'       admissible splits (e.g. \eqn{|R| - 1} for a contiguous range on a line
#'       with \code{vic_mode = "connected"}: a rate per link between
#'       neighbouring islands; \eqn{2^{|R|-1} - 1} with
#'       \code{vic_mode = "all"}).
#'   }
#' @param clado_single What happens when a lineage confined to a single island
#'   speciates in situ: \code{"split"} (default, as in DAISIE: two new endemic
#'   species replace the parent) or \code{"bud"} (the parent keeps its
#'   identity, so a non-endemic remains non-endemic, and one new endemic
#'   species buds off). Widespread lineages always bud.
#' @param max_range_size Maximum number of islands a lineage may occupy.
#'   Events that would create a larger range (dispersal, recolonisation of a
#'   non-endemic on a new island) are not allowed. Reduces the number of range
#'   states.
#' @param connected_ranges Logical. If \code{TRUE}, only ranges that are
#'   connected in the adjacency graph are allowed; dispersal,
#'   recolonisation, local extinction, island-level anagenesis and vicariance
#'   events that would create a disconnected range are not allowed. Reduces the number of range states
#'   (e.g. to \eqn{n(n+1)/2} for islands on a line).
#' @param island_ages Optional vector of island ages (time before present at
#'   which each island emerges). \code{NULL} (default) means that all islands
#'   are as old as the simulation. Before it emerges an island cannot receive
#'   colonists or dispersers.
#' @param island_names Optional character vector of island names.
#'
#' @return An object of class \code{DAISIE_A_pars}: a named list with the
#'   expanded per-island vectors (\code{lac}, \code{mu}, \code{K},
#'   \code{gam}, \code{laa_isl}), the scalars (\code{laa}, \code{lav},
#'   \code{disp}), the connectivity matrix \code{C}, the dispersal rate matrix
#'   \code{D = disp * C}, the logical \code{adjacency} matrix and the settings.
#' @author Rampal S. Etienne
#' @family DAISIE_A
#' @export
#' @examples
#' # Three islands on a line, equal rates on all islands
#' pars <- DAISIE_A_create_pars(
#'   n_islands = 3, lac = 0.5, mu = 0.2, K = Inf, gam = 0.02, laa = 0.5,
#'   lav = 0.1, disp = 0.3, connectivity = "line"
#' )
#' # Island-specific extinction rates and a general dispersal matrix
#' C <- matrix(c(0, 1, 0.2,
#'               1, 0, 0.5,
#'               0.2, 0.5, 0), 3, 3, byrow = TRUE)
#' pars2 <- DAISIE_A_create_pars(
#'   n_islands = 3, lac = 0.5, mu = c(0.1, 0.2, 0.4), K = 20, gam = 0.02,
#'   laa = 0.5, lav = 0.1, disp = 0.3, connectivity = C
#' )
DAISIE_A_create_pars <- function(n_islands,
                                 lac,
                                 mu,
                                 K = Inf,
                                 gam,
                                 laa,
                                 laa_isl = 0,
                                 lav = 0,
                                 disp = 0,
                                 connectivity = "full",
                                 vic_mode = "connected",
                                 vic_per = "lineage",
                                 clado_single = "split",
                                 max_range_size = n_islands,
                                 connected_ranges = FALSE,
                                 island_ages = NULL,
                                 island_names = NULL) {
  if (length(n_islands) != 1 || !is.numeric(n_islands) ||
      n_islands < 1 || n_islands != round(n_islands)) {
    stop("'n_islands' must be a single positive integer.")
  }
  n_islands <- as.integer(n_islands)
  if (n_islands > 30) {
    stop("'n_islands' larger than 30 is not supported.")
  }
  expand <- function(x, name) {
    if (!is.numeric(x) || any(is.na(x))) {
      stop("'", name, "' must be numeric.")
    }
    if (length(x) == 1) x <- rep(x, n_islands)
    if (length(x) != n_islands) {
      stop("'", name, "' must have length 1 or n_islands (", n_islands, ").")
    }
    if (any(x < 0)) stop("'", name, "' must be non-negative.")
    as.numeric(x)
  }
  lac <- expand(lac, "lac")
  mu <- expand(mu, "mu")
  K <- expand(K, "K")
  if (any(K == 0)) stop("'K' must be positive.")
  gam <- expand(gam, "gam")
  laa_isl <- expand(laa_isl, "laa_isl")
  for (nm in c("laa", "lav", "disp")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1 || is.na(v) || v < 0) {
      stop("'", nm, "' must be a single non-negative number.")
    }
  }
  if (is.character(connectivity)) {
    C <- DAISIE_A_connectivity(n_islands = n_islands, type = connectivity)
  } else {
    C <- as.matrix(connectivity)
    if (!is.numeric(C) || nrow(C) != n_islands || ncol(C) != n_islands) {
      stop("'connectivity' must be a character string or a ",
           "n_islands x n_islands numeric matrix.")
    }
    if (any(is.na(C)) || any(C < 0)) {
      stop("'connectivity' must be non-negative.")
    }
    diag(C) <- 0
  }
  dimnames(C) <- NULL
  vic_mode <- match.arg(vic_mode, c("connected", "all", "single"))
  vic_per <- match.arg(vic_per, c("lineage", "split"))
  clado_single <- match.arg(clado_single, c("split", "bud"))
  if (length(max_range_size) != 1 || max_range_size < 1 ||
      max_range_size != round(max_range_size)) {
    stop("'max_range_size' must be a single positive integer.")
  }
  max_range_size <- as.integer(min(max_range_size, n_islands))
  if (!is.logical(connected_ranges) || length(connected_ranges) != 1 ||
      is.na(connected_ranges)) {
    stop("'connected_ranges' must be TRUE or FALSE.")
  }
  if (!is.null(island_ages)) {
    if (!is.numeric(island_ages) || length(island_ages) != n_islands ||
        any(is.na(island_ages)) || any(island_ages <= 0)) {
      stop("'island_ages' must be NULL or a positive vector of length ",
           "n_islands.")
    }
  }
  if (is.null(island_names)) {
    island_names <- paste0("island_", seq_len(n_islands))
  }
  if (length(island_names) != n_islands) {
    stop("'island_names' must have length n_islands.")
  }
  adjacency <- (C + t(C)) > 0
  pars <- list(
    n_islands = n_islands,
    island_names = as.character(island_names),
    lac = lac,
    mu = mu,
    K = K,
    gam = gam,
    laa = laa,
    laa_isl = laa_isl,
    lav = lav,
    disp = disp,
    C = C,
    D = disp * C,
    adjacency = adjacency,
    vic_mode = vic_mode,
    vic_per = vic_per,
    clado_single = clado_single,
    max_range_size = max_range_size,
    connected_ranges = connected_ranges,
    island_ages = island_ages
  )
  class(pars) <- "DAISIE_A_pars"
  return(pars)
}

#' Check whether an object is a valid DAISIE_A parameter object
#'
#' @param pars Object to check.
#' @return \code{TRUE} or \code{FALSE}
#' @keywords internal
DAISIE_A_are_pars <- function(pars) {
  if (!inherits(pars, "DAISIE_A_pars")) return(FALSE)
  needed <- c("n_islands", "lac", "mu", "K", "gam", "laa", "lav", "D",
              "C", "adjacency", "vic_mode", "clado_single", "max_range_size",
              "connected_ranges")
  if (!all(needed %in% names(pars))) return(FALSE)
  n <- pars$n_islands
  all(length(pars$lac) == n, length(pars$mu) == n, length(pars$K) == n,
      length(pars$gam) == n, all(dim(pars$D) == c(n, n)),
      length(DAISIE_A_laa_isl(pars)) == n)
}

#' @export
print.DAISIE_A_pars <- function(x, ...) {
  cat("DAISIE_A parameters for an archipelago of", x$n_islands, "islands\n")
  tab <- rbind(lac = x$lac, mu = x$mu, K = x$K, gam = x$gam,
               laa_isl = DAISIE_A_laa_isl(x))
  colnames(tab) <- x$island_names
  print(tab)
  cat("laa =", x$laa, " lav =", x$lav, " disp =", x$disp, "\n")
  cat("connectivity matrix C:\n")
  print(x$C)
  cat("vic_mode =", x$vic_mode, " vic_per =", DAISIE_A_vic_per(x),
      " clado_single =", x$clado_single,
      " max_range_size =", x$max_range_size,
      " connected_ranges =", x$connected_ranges, "\n")
  if (!is.null(x$island_ages)) cat("island ages:", x$island_ages, "\n")
  invisible(x)
}

#' Standard connectivity matrices for DAISIE_A
#'
#' @description Creates a connectivity matrix \eqn{C} for an archipelago,
#'   with \eqn{C_{ij}} the relative dispersal weight from island \eqn{i} to
#'   island \eqn{j}.
#'
#' @param n_islands Number of islands.
#' @param type One of \code{"full"} (all pairs connected with weight 1),
#'   \code{"line"} (stepping stone: only adjacent islands \eqn{i} and
#'   \eqn{i \pm 1} are connected), \code{"ring"} (as line, but the first and
#'   last island are also connected), or \code{"none"} (no dispersal).
#'
#' @return An \code{n_islands} by \code{n_islands} matrix with zero diagonal.
#' @family DAISIE_A
#' @export
#' @examples
#' DAISIE_A_connectivity(4, "line")
DAISIE_A_connectivity <- function(n_islands, type = "full") {
  type <- match.arg(type, c("full", "line", "ring", "none"))
  C <- matrix(0, n_islands, n_islands)
  if (type == "full") {
    C[] <- 1
  } else if (type %in% c("line", "ring") && n_islands > 1) {
    for (i in seq_len(n_islands - 1)) {
      C[i, i + 1] <- 1
      C[i + 1, i] <- 1
    }
    if (type == "ring" && n_islands > 2) {
      C[1, n_islands] <- 1
      C[n_islands, 1] <- 1
    }
  }
  diag(C) <- 0
  return(C)
}

#' Distance-based connectivity matrix for DAISIE_A
#'
#' @description Creates a connectivity matrix from a matrix of inter-island
#'   distances, so that the whole dispersal matrix is described by one or two
#'   parameters instead of \eqn{n(n-1)}.
#'
#' @param distances Symmetric (or not) matrix of distances between islands.
#' @param scale Scale parameter of the kernel.
#' @param kernel \code{"exponential"}: \eqn{C_{ij} = \exp(-d_{ij}/s)};
#'   \code{"power"}: \eqn{C_{ij} = d_{ij}^{-s}}.
#' @param threshold Distances larger than this are set to have zero
#'   connectivity (default \code{Inf}).
#'
#' @return Connectivity matrix with zero diagonal.
#' @family DAISIE_A
#' @export
#' @examples
#' xy <- cbind(c(0, 1, 3), c(0, 0, 1))
#' DAISIE_A_connectivity_distance(as.matrix(dist(xy)), scale = 1)
DAISIE_A_connectivity_distance <- function(distances,
                                           scale = 1,
                                           kernel = "exponential",
                                           threshold = Inf) {
  kernel <- match.arg(kernel, c("exponential", "power"))
  distances <- as.matrix(distances)
  if (nrow(distances) != ncol(distances)) {
    stop("'distances' must be a square matrix.")
  }
  if (kernel == "exponential") {
    C <- exp(-distances / scale)
  } else {
    C <- distances^(-scale)
  }
  C[distances > threshold] <- 0
  diag(C) <- 0
  dimnames(C) <- NULL
  return(C)
}

#' Unit to which the vicariance rate applies
#'
#' @param pars A \code{DAISIE_A_pars} object.
#' @return \code{"lineage"} or \code{"split"} (\code{"lineage"} for
#'   parameter objects created before this setting existed).
#' @keywords internal
DAISIE_A_vic_per <- function(pars) {
  if (is.null(pars$vic_per)) "lineage" else pars$vic_per
}

#' Island-level anagenesis rates
#'
#' @param pars A \code{DAISIE_A_pars} object.
#' @return Vector of island-level anagenesis rates (zeros for parameter objects
#'   created before this process existed).
#' @keywords internal
DAISIE_A_laa_isl <- function(pars) {
  if (is.null(pars$laa_isl)) rep(0, pars$n_islands) else pars$laa_isl
}
