TRAISIE_n_series <- function(
    datalist, clade, parameter,
    num_observed_states, num_hidden_states,
    atol = 1e-10, rtol = 1e-10,
    max_coefficients = 256L,
    max_assignments = 4096L) {

  S <- num_observed_states
  H <- num_hidden_states
  n <- S * H

  missing <- clade$missing_species
  traits <- as.numeric(clade$traits)
  stac <- clade$stac
  brts <- clade$branching_times
  phy <- clade$phylogeny




  if (any(diag(parameter[[5]]) != 0)) {
    stop("q must have zero diagonal.")
  }

  if (length(brts) < 2 ||
      any(!is.finite(brts)) ||
      any(brts < 0) ||
      brts[1] < brts[2]) {
    stop("Invalid branching times.")
  }

  # One coefficient for each combination of missing-species degrees.
  K <- prod(missing + 1)

  if (K > max_coefficients) {
    stop("Too many series coefficients.")
  }

  degrees <- as.matrix(
    expand.grid(lapply(missing, function(m) 0:m))
  )

  strides <- cumprod(c(1, head(missing + 1, -1)))

  # Precompute the coefficient pairs used in each convolution.
  pairs <- lapply(seq_len(K), function(k) {
    keep <- which(
      rowSums(sweep(degrees, 2, degrees[k, ], ">")) == 0
    )

    remainder <- matrix(
      degrees[k, ],
      nrow = length(keep),
      ncol = S,
      byrow = TRUE
    ) - degrees[keep, , drop = FALSE]

    list(
      left = keep,
      right = as.integer(1 + remainder %*% strides)
    )
  })

  # Truncated multivariate polynomial multiplication.
  multiply <- function(a, b) {
    columns <- max(ncol(a), ncol(b))

    if (ncol(a) == 1L && columns > 1L) {
      a <- a[, rep(1L, columns), drop = FALSE]
    }

    if (ncol(b) == 1L && columns > 1L) {
      b <- b[, rep(1L, columns), drop = FALSE]
    }

    out <- matrix(0, K, columns)

    for (k in seq_len(K)) {
      pair <- pairs[[k]]

      out[k, ] <- colSums(
        a[pair$left, , drop = FALSE] *
          b[pair$right, , drop = FALSE]
      )
    }

    out
  }

  # A numerical constant has only a degree-zero coefficient.
  constant <- function(value) {
    out <- matrix(0, K, length(value))
    out[1, ] <- value
    out
  }

  scale_columns <- function(x, rates) {
    sweep(x, 2, rates, "*")
  }

  # Use one common layout throughout the series calculation.
  indices <- setNames(
    lapply(0:4, function(i) i * n + seq_len(n)),
    c("DE", "DM1", "DM2", "DM3", "E")
  )

  indices$DA2 <- 5 * n + 1
  indices$DA3 <- 5 * n + 2

  empty <- function() {
    matrix(0, K, 5 * n + 2)
  }

  lc <- parameter[[1]]
  mu <- parameter[[2]]
  la <- parameter[[4]]
  q <- parameter[[5]]

  pq <- parameter[[6]] * q
  opq <- (1 - parameter[[6]]) * q
  loss <- lc + mu + la + rowSums(q)

  # Reuse TRAISIE's existing mainland-ancestor weights.
  weights_function <- get0(
    "TRAISIE_weights",
    mode = "function",
    inherits = TRUE
  )

  if (is.null(weights_function)) {
    weights_function <- getFromNamespace(
      "TRAISIE_weights", "DAISIE"
    )
  }

  weights <- weights_function(
    clade$root_state, S, H, datalist
  )

  if (length(weights) != n ||
      any(!is.finite(weights)) ||
      any(weights < 0) ||
      abs(sum(weights) - 1) > 1e-8) {
    stop("Invalid mainland ancestor weights.")
  }

  # E_s(0) = z_observed_state(s).
  # If no species are missing in a state, its z is set to zero.
  z <- matrix(0, K, n)

  for (s in seq_len(n)) {
    observed <- (s - 1) %/% H + 1

    if (missing[observed] > 0) {
      z[1 + strides[observed], s] <- 1
    }
  }

  ancestor_coefficient <- function(assigned, ancestor) {
    g <- numeric(n)
    g[ancestor] <- parameter[[3]][ancestor]
    G <- sum(g)

    rhs <- function(t, y, interval) {
      x <- matrix(y, K)

      block <- function(name) {
        x[, indices[[name]], drop = FALSE]
      }

      dx <- empty()

      E <- block("E")
      DE <- block("DE")
      EE <- multiply(E, E)

      f <- constant(mu) +
        scale_columns(E, la) +
        scale_columns(EE, lc) +
        E %*% t(pq)

      h <- scale_columns(DE, la) +
        2 * scale_columns(multiply(DE, E), lc) +
        DE %*% t(pq)

      reset <- function(M) {
        -G * M +
          (M %*% g)[, rep(1L, n), drop = FALSE]
      }

      nonendemic <- function(M, A) {
        -scale_columns(M, loss) +
          multiply(f, A) +
          M %*% t(opq)
      }

      dx[, indices$E] <- constant(mu) -
        scale_columns(E, lc + mu + rowSums(q)) +
        scale_columns(EE, lc) +
        E %*% t(q)

      if (interval != 4) {
        dx[, indices$DE] <-
          -scale_columns(DE, lc + mu + rowSums(q)) +
          2 * scale_columns(multiply(DE, E), lc) +
          DE %*% t(q)

        M3 <- block("DM3")
        A3 <- block("DA3")

        dx[, indices$DM3] <-
          nonendemic(M3, A3) + reset(M3)

        dx[, indices$DA3] <- -G * A3 + M3 %*% g
      }

      if (interval == 2) {
        M2 <- block("DM2")

        dx[, indices$DM2] <-
          -scale_columns(M2, loss + G) +
          multiply(h, block("DA3")) +
          M2 %*% t(opq)
      }

      if (interval == 3) {
        M1 <- block("DM1")
        M2 <- block("DM2")
        A2 <- block("DA2")

        dx[, indices$DM1] <-
          nonendemic(M1, A2) - G * M1 +
          (M2 %*% g)[, rep(1L, n), drop = FALSE]

        dx[, indices$DM2] <-
          nonendemic(M2, A2) +
          multiply(h, block("DA3")) +
          reset(M2)

        dx[, indices$DA2] <- -G * A2 + M2 %*% g
      }

      if (interval == 4) {
        M1 <- block("DM1")

        # During interval4, the DA2 slot holds DA1.
        A1 <- block("DA2")

        dx[, indices$DM1] <-
          nonendemic(M1, A1) + reset(M1)

        dx[, indices$DA2] <- -G * A1 + M1 %*% g
      }

      list(as.vector(dx))
    }

    solve <- function(x, duration, interval) {
      if (!is.finite(duration) || duration < -1e-8) {
        stop("Invalid branch duration.")
      }

      if (duration <= 0) {
        return(x)
      }

      final_coefficients <- .Call(
        "TRAISIE_series_cpp",
        as.vector(x),
        parameter,
        as.numeric(g),
        pairs,
        as.integer(interval),
        as.numeric(duration),
        as.numeric(atol),
        as.numeric(rtol),
        PACKAGE = "DAISIE"
      )

      matrix(final_coefficients, nrow = K)
    }

    tip <- function(trait) {
      x <- empty()
      x[, indices$E] <- z

      x[1, indices$DE] <- rep(
        as.numeric(seq_len(S) - 1 == trait),
        each = H
      )

      x[1, indices$DA3] <- 1
      x
    }

    if (length(brts) == 2L) {
      if (length(assigned) != 1L || stac == 6) {
        stop("Singleton series supports stac 2 or 5 and one tip.")
      }

      x <- tip(assigned)
      interval <- if (stac == 2) 2 else 3
      x <- solve(x, brts[2], interval)

    } else {
      if (!stac %in% c(2, 6)) {
        stop("Tree series currently supports stac 2 or 6.")
      }

      if (!inherits(phy, "phylo") ||
          length(phy$tip.label) != length(assigned)) {
        stop("Supply a phylo tree with one trait per tip.")
      }

      if (is.null(phy$edge.length) ||
          any(!is.finite(phy$edge.length)) ||
          any(phy$edge.length < 0)) {
        stop("Supply finite nonnegative tree edge lengths.")
      }

      roots <- setdiff(
        unique(phy$edge[, 1]),
        phy$edge[, 2]
      )

      if (length(roots) != 1L) {
        stop("Tree must have one root.")
      }

      visit <- function(node) {
        if (node <= length(assigned)) {
          return(list(
            x = tip(assigned[node]),
            height = 0
          ))
        }

        edges <- which(phy$edge[, 1] == node)

        if (length(edges) != 2L) {
          stop("Tree must be bifurcating.")
        }

        children <- lapply(edges, function(e) {
          child <- visit(phy$edge[e, 2])

          list(
            x = solve(child$x, phy$edge.length[e], 1),
            height = child$height + phy$edge.length[e]
          )
        })

        if (abs(children[[1]]$height - children[[2]]$height) >
            1e-7 * max(1, children[[1]]$height)) {
          stop("Tree must be ultrametric.")
        }

        x <- children[[1]]$x

        x[, indices$DE] <- scale_columns(
          multiply(
            x[, indices$DE, drop = FALSE],
            children[[2]]$x[, indices$DE, drop = FALSE]
          ),
          lc
        )

        list(x = x, height = children[[1]]$height)
      }

      root <- visit(roots)

      if (abs(root$height - brts[3]) >
          1e-3 * max(1, brts[3])) {
        stop("Tree crown age must match brts[3].")
      }

      x <- root$x

      x[, indices$DM2] <- multiply(
        x[, indices$DE, drop = FALSE],
        x[, indices$DA3, drop = FALSE]
      )

      interval <- if (stac == 2) 2 else 3
      x <- solve(x, brts[2] - brts[3], interval)
    }

    if (stac == 2) {
      event <- x[, indices$DM2, drop = FALSE] %*% g

      x[, indices$DM1] <-
        event[, rep(1L, n), drop = FALSE]

      x[, indices$DA2] <- event
    }

    x <- solve(x, brts[1] - brts[2], 4)

    # The last row is the requested multivariate coefficient.
    x[K, indices$DA2]
  }

  evaluate_assignment <- function(assigned) {
    coefficients <- vapply(seq_len(n), function(a) {
      if (weights[a] == 0) {
        return(0)
      }

      ancestor_coefficient(assigned, a)
    }, numeric(1))

    coefficient <- sum(weights * coefficients)

    if (!is.finite(coefficient) || coefficient < 0) {
      stop("Invalid coefficient; check tolerances and model equations.")
    }

    sampled <- tabulate(assigned + 1L, nbins = S)

    log(coefficient) -
      sum(lchoose(sampled + missing, missing))
  }

  # Unknown tips require assignment-specific normalization.
  unknown <- which(is.na(traits))

  if (S ^ length(unknown) > max_assignments) {
    stop("Too many unknown-tip assignments.")
  }

  add_logs <- function(a, b) {
    if (a == -Inf) return(b)
    if (b == -Inf) return(a)

    scale <- max(a, b)
    scale + log(exp(a - scale) + exp(b - scale))
  }

  enumerate <- function(assigned, k = 1L) {
    if (k > length(unknown)) {
      return(evaluate_assignment(assigned))
    }

    total <- -Inf

    for (state in seq_len(S) - 1L) {
      assigned[unknown[k]] <- state
      total <- add_logs(total, enumerate(assigned, k + 1L))
    }

    total
  }

  enumerate(traits)
}
