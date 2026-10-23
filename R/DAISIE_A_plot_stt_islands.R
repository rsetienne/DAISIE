#' Plot per-island species-through-time (STT) plots of archipelago simulations
#'
#' @description Produces one panel per island with the median (line) and the
#'   2.5--97.5th (light grey) and 25--75th (dark grey) percentiles across
#'   replicates of the number of species on that island, as in
#'   \code{\link{DAISIE_plot_sims}()}.
#'
#' @param island_replicates Output of \code{\link{DAISIE_A_sim}()}.
#' @param what Which species to count: \code{"total"} (default),
#'   \code{"endemic"} (\code{nA + nC}), \code{"nonendemic"} (\code{nI}) or
#'   \code{"SIE"} (single-island endemics).
#' @param log_scale Logical; plot \eqn{\log(N + 1)} on the y-axis.
#'
#' @return Invisibly, a list with the per-island quantiles.
#' @family DAISIE_A
#' @export
#' @examples
#' pars <- DAISIE_A_create_pars(n_islands = 3, lac = 0.5, mu = 0.2,
#'   gam = 0.02, laa = 0.5, lav = 0.2, disp = 0.3, connectivity = "line")
#' set.seed(1)
#' sims <- DAISIE_A_sim(time = 2, M = 50, pars = pars, replicates = 5,
#'   verbose = FALSE)
#' DAISIE_A_plot_stt_islands(sims)
DAISIE_A_plot_stt_islands <- function(island_replicates,
                                      what = "total",
                                      log_scale = FALSE) {
  what <- match.arg(what, c("total", "endemic", "nonendemic", "SIE"))
  first <- island_replicates[[1]][[1]]
  n <- first$n_islands
  if (is.null(n)) stop("'island_replicates' must be the output of DAISIE_A_sim().")
  pick <- function(m) {
    switch(what,
           total = m[, "nI"] + m[, "nA"] + m[, "nC"],
           endemic = m[, "nA"] + m[, "nC"],
           nonendemic = m[, "nI"],
           SIE = m[, "nSIE"])
  }
  times <- first$stt_islands[[1]][, "Time"]
  res <- list()
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  nc <- ceiling(sqrt(n))
  graphics::par(mfrow = c(ceiling(n / nc), nc), mar = c(4, 4, 2, 1))
  for (i in seq_len(n)) {
    vals <- sapply(island_replicates, function(r) pick(r[[1]]$stt_islands[[i]]))
    vals <- matrix(vals, nrow = length(times))
    if (log_scale) vals <- log(vals + 1)
    q <- apply(vals, 1, stats::quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    res[[first$island_names[i]]] <- q
    graphics::plot(NULL, xlim = rev(range(times)), ylim = c(0, max(q, 1)),
                   xlab = "Time before present",
                   ylab = if (log_scale) "log(N + 1)" else "N",
                   main = first$island_names[i])
    graphics::polygon(c(times, rev(times)), c(q[1, ], rev(q[5, ])),
                      col = "lightgrey", border = NA)
    graphics::polygon(c(times, rev(times)), c(q[2, ], rev(q[4, ])),
                      col = "darkgrey", border = NA)
    graphics::lines(times, q[3, ], lwd = 2)
  }
  invisible(res)
}
