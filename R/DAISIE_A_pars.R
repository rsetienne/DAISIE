#' Create the parameter object of the archipelago DAISIE model (DAISIE_A)
#'
#' @description
#' Collects, expands and checks all rates and settings of the archipelago
#' extension of DAISIE (DAISIE_A, version 10: the divergence-state model). In
#' DAISIE_A every lineage on the archipelago has a range: the non-empty set of
#' islands on which it occurs. Lineages expand their range by colonisation from
#' the mainland and by dispersal between islands, contract it by local
#' extinction (losing the last island is global extinction) and speciate within
#' an island (in-situ cladogenesis, by budding).
#'
#' Speciation between islands (vicariance) and anagenesis are not separate
#' rates: they result from divergence states. Every pair of adjacent islands in
#' the range of a lineage has a state \eqn{d_{ij} \in \{0, \ldots, k\}}, and
#' so has every pair of the mainland and an island of the range of a
#' non-endemic lineage (\eqn{e_i}, for islands with \eqn{\gamma_i > 0}). A
#' state below \eqn{k} advances by one at rate \eqn{\alpha}; gene flow between
#' the two populations (dispersal between the islands, colonisation of the
#' island from the mainland) resets it to 0. When the network of pairs below
#' \eqn{k} falls apart, the lineage splits at once into one species per part:
#' for a non-endemic lineage the part with the mainland stays non-endemic, and
#' if that part has no islands the lineage becomes endemic (anagenesis). With
#' one island and \eqn{k = 1} the model is DAISIE with anagenesis rate
#' \eqn{\alpha} (with \code{clado_single = "split"}). See the report for the
#' complete description (gap pairs of fragmented ranges, the copy rule for
#' founders, multifurcations).
#'
#' All island-specific rates (\code{lac}, \code{mu}, \code{K}, \code{gam},
#' \code{beta}) can be given as a single number, which is then used for every
#' island, or as a vector with one value per island.
#'
#' @param n_islands Integer, number of islands in the archipelago.
#' @param lac Within-island (in-situ) cladogenesis rate per island population
#'   (\eqn{\lambda^c_i}). Scalar or vector of length \code{n_islands}.
#' @param mu Local extinction rate per island population (\eqn{\mu_i}).
#'   For a lineage on a single island this is the global extinction rate.
#'   Scalar or vector of length \code{n_islands}.
#' @param K Clade-specific carrying capacity per island (\eqn{K_i}). Events
#'   that add a species to island \eqn{i} (colonisation of the island by a
#'   lineage that does not occur there, dispersal to it, in-situ speciation on
#'   it) are multiplied by \eqn{\max(0, 1 - n_i / K_i)}, where \eqn{n_i} is the
#'   number of species of the clade on island \eqn{i}. Recolonisation (by a
#'   mainland species already present on the island), gene flow and divergence
#'   are not diversity dependent. Use \code{Inf} for no diversity dependence.
#' @param gam Colonisation rate from the mainland to island \eqn{i}, per
#'   mainland species (\eqn{\gamma_i}). Scalar or vector of length
#'   \code{n_islands}; a scalar is the rate per island.
#' @param alpha Rate at which a divergence state below \code{k} advances by one
#'   (\eqn{\alpha}).
#' @param k Number of divergence steps after which two populations are
#'   separate species (integer, at least 1).
#' @param disp Dispersal rate multiplier. Dispersal from island \eqn{i} to
#'   island \eqn{j} occurs at rate \code{disp} \eqn{\times C_{ij}} per
#'   population on \eqn{i}; to an island outside the range it expands the
#'   range, to an island in the range it is gene flow.
#' @param connectivity Either a character string (\code{"full"}, \code{"line"},
#'   \code{"ring"} or \code{"none"}), passed to
#'   \code{\link{DAISIE_A_connectivity}()}, or a non-negative
#'   \code{n_islands} by \code{n_islands} matrix \eqn{C} of relative
#'   dispersal weights from row island to column island (diagonal ignored).
#'   The non-zero pattern of \eqn{C + C^T} defines which islands are adjacent,
#'   i.e. which pairs of islands have a divergence state.
#' @param clado_single What happens when a lineage confined to a single island
#'   speciates in situ: \code{"bud"} (default: the parent keeps its identity,
#'   so a non-endemic remains non-endemic, and one new endemic species buds
#'   off) or \code{"split"} (as in DAISIE: two new endemic species replace the
#'   parent). Widespread lineages always bud.
#' @param beta Rate per island population at which the within-island
#'   divergence state of the protracted simulation option advances
#'   (\eqn{\beta_i}; only used by \code{\link{DAISIE_A_sim}(protracted = TRUE)},
#'   where it replaces \code{lac}). \code{NULL} (default) means equal to
#'   \code{lac}.
#' @param island_ages Optional vector of island ages (time before present at
#'   which each island emerges). \code{NULL} (default) means that all islands
#'   are as old as the simulation. Before it emerges an island cannot receive
#'   colonists or dispersers.
#' @param submergence_ages Optional vector with, for each island, the time
#'   before present at which it submerged, or \code{NA} for an island that
#'   still exists. \code{NULL} (default) means that no island has submerged.
#'   When an island submerges, all its populations are lost, as with local
#'   extinction.
#' @param island_names Optional character vector of island names.
#'
#' @return An object of class \code{DAISIE_A_pars}: a named list with the
#'   expanded per-island vectors (\code{lac}, \code{mu}, \code{K},
#'   \code{gam}, \code{beta}), the scalars (\code{alpha}, \code{k},
#'   \code{disp}), the connectivity matrix \code{C}, the dispersal rate matrix
#'   \code{D = disp * C}, the logical \code{adjacency} matrix and the settings.
#' @author Rampal S. Etienne
#' @family DAISIE_A
#' @export
#' @examples
#' # Three islands on a line, equal rates on all islands
#' pars <- DAISIE_A_create_pars(
#'   n_islands = 3, lac = 0.5, mu = 0.2, K = Inf, gam = 0.02, alpha = 0.5,
#'   k = 2, disp = 0.3, connectivity = "line"
#' )
#' # Island-specific extinction rates and a general dispersal matrix
#' C <- matrix(c(0, 1, 0.2,
#'               1, 0, 0.5,
#'               0.2, 0.5, 0), 3, 3, byrow = TRUE)
#' pars2 <- DAISIE_A_create_pars(
#'   n_islands = 3, lac = 0.5, mu = c(0.1, 0.2, 0.4), K = 20, gam = 0.02,
#'   alpha = 0.5, k = 1, disp = 0.3, connectivity = C
#' )
DAISIE_A_create_pars <- function(n_islands,
                                 lac,
                                 mu,
                                 K = Inf,
                                 gam,
                                 alpha,
                                 k = 1,
                                 disp = 0,
                                 connectivity = "full",
                                 clado_single = "bud",
                                 beta = NULL,
                                 island_ages = NULL,
                                 submergence_ages = NULL,
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
  beta <- if (is.null(beta)) lac else expand(beta, "beta")
  for (nm in c("alpha", "disp")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1 || is.na(v) || v < 0 ||
        is.infinite(v)) {
      stop("'", nm, "' must be a single non-negative number.")
    }
  }
  if (!is.numeric(k) || length(k) != 1 || is.na(k) || k < 1 ||
      k != round(k)) {
    stop("'k' must be a single integer of at least 1.")
  }
  k <- as.integer(k)
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
  clado_single <- match.arg(clado_single, c("bud", "split"))
  if (!is.null(island_ages)) {
    if (!is.numeric(island_ages) || length(island_ages) != n_islands ||
        any(is.na(island_ages)) || any(island_ages <= 0)) {
      stop("'island_ages' must be NULL or a positive vector of length ",
           "n_islands.")
    }
  }
  adjacency <- (C + t(C)) > 0
  if (!is.null(submergence_ages)) {
    if (!is.numeric(submergence_ages) && !all(is.na(submergence_ages))) {
      stop("'submergence_ages' must be NULL or a numeric vector.")
    }
    submergence_ages <- as.numeric(submergence_ages)
    if (length(submergence_ages) != n_islands ||
        any(submergence_ages[!is.na(submergence_ages)] <= 0)) {
      stop("'submergence_ages' must be NULL or a vector of length n_islands ",
           "with positive values or NA (not submerged).")
    }
    if (!is.null(island_ages) &&
        any(submergence_ages >= island_ages, na.rm = TRUE)) {
      stop("An island must emerge before it submerges: 'submergence_ages' ",
           "must be smaller than 'island_ages'.")
    }
    if (all(is.na(submergence_ages))) submergence_ages <- NULL
  }
  if (is.null(island_names)) {
    island_names <- paste0("island_", seq_len(n_islands))
  }
  if (length(island_names) != n_islands) {
    stop("'island_names' must have length n_islands.")
  }
  pars <- list(
    n_islands = n_islands,
    island_names = as.character(island_names),
    lac = lac,
    mu = mu,
    K = K,
    gam = gam,
    beta = beta,
    alpha = alpha,
    k = k,
    disp = disp,
    C = C,
    D = disp * C,
    adjacency = adjacency,
    clado_single = clado_single,
    island_ages = island_ages,
    submergence_ages = submergence_ages,
    version = 10L
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
  needed <- c("n_islands", "lac", "mu", "K", "gam", "beta", "alpha", "k",
              "D", "C", "adjacency", "clado_single")
  if (!all(needed %in% names(pars))) return(FALSE)
  n <- pars$n_islands
  all(length(pars$lac) == n, length(pars$mu) == n, length(pars$K) == n,
      length(pars$gam) == n, length(pars$beta) == n,
      all(dim(pars$D) == c(n, n)))
}

#' @export
print.DAISIE_A_pars <- function(x, ...) {
  cat("DAISIE_A parameters for an archipelago of", x$n_islands, "islands\n")
  tab <- rbind(lac = x$lac, mu = x$mu, K = x$K, gam = x$gam, beta = x$beta)
  colnames(tab) <- x$island_names
  print(tab)
  cat("alpha =", x$alpha, " k =", x$k, " disp =", x$disp,
      " clado_single =", x$clado_single, "\n")
  cat("connectivity matrix C:\n")
  print(x$C)
  if (!is.null(x$island_ages)) cat("island ages:", x$island_ages, "\n")
  if (!is.null(x$submergence_ages)) {
    cat("submergence ages:", x$submergence_ages, "\n")
  }
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
