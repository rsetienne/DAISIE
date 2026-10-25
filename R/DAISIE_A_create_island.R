#' Convert the final state of an archipelago simulation into output
#'
#' @description Internal function. Closes the genealogy, reconstructs the
#'   phylogeny of every colonisation event with extant descendants, and
#'   creates (i) the standard DAISIE output for the archipelago treated as a
#'   single island (via \code{\link{DAISIE_create_island}()}) and (ii) the
#'   archipelago-specific output (ranges, trees, per-island STT).
#'
#' @param st Final state list.
#' @param stt_table Archipelago-level STT table (DAISIE format).
#' @param stt_islands Per-island STT table.
#' @param total_time Simulation time.
#' @param mainland_n Number of mainland species.
#' @param pars A \code{DAISIE_A_pars} object.
#'
#' @return List, see \code{\link{DAISIE_A_sim_core}()}.
#' @keywords internal
DAISIE_A_create_island <- function(st,
                                   stt_table,
                                   stt_islands,
                                   total_time,
                                   mainland_n,
                                   pars) {
  n <- pars$n_islands
  spec <- st$spec
  L <- nrow(spec)
  # close the genealogy: one tip node per extant species
  tip_node <- rep(NA_integer_, L)
  for (r in seq_len(L)) {
    st <- DAISIE_A_add_node(st, parent = spec[r, .DA_NODE],
                            time = total_time, kind = 5,
                            side = spec[r, .DA_SIDE],
                            clade = spec[r, .DA_CLADE],
                            species = spec[r, .DA_ID])
    tip_node[r] <- st$n_gen
  }
  gen <- st$gen[seq_len(st$n_gen), , drop = FALSE]

  type_char <- c("I", "A", "C")
  ana_char <- c(NA, "Immig_parent", "Clado_extinct")
  occ <- spec[, .DA_NFIX + seq_len(n), drop = FALSE]

  # reconstructed trees of all colonisation events with extant species
  clade_ids <- unique(spec[, .DA_CLADE])
  clade_ids <- clade_ids[order(gen[match(clade_ids, gen[, "clade"]), "time"])]
  daisie_code <- rep(NA_character_, L)
  daisie_btime <- rep(NA_real_, L)
  daisie_type <- type_char[spec[, .DA_TYPE]]
  clades <- list()
  for (k in seq_along(clade_ids)) {
    rows <- which(spec[, .DA_CLADE] == clade_ids[k])
    rec <- DAISIE_A_reconstruct_clade(
      gen = gen,
      clade = clade_ids[k],
      total_time = total_time
    )
    tip_pos <- match(tip_node[rows], rec$tips)
    mixed <- any(spec[rows, .DA_TYPE] == 1) && any(spec[rows, .DA_TYPE] != 1)
    if (length(rows) > 1) {
      # DAISIE view: all members of a multi-species clade are cladogenetic
      daisie_type[rows] <- "C"
      daisie_code[rows] <- rec$code[tip_pos]
      daisie_btime[rows] <- rec$daisie_btime[tip_pos]
    }
    tips <- data.frame(
      species = spec[rows, .DA_ID],
      type = type_char[spec[rows, .DA_TYPE]],
      range = DAISIE_A_range_string(occ[rows, , drop = FALSE]),
      range_code = DAISIE_A_range_code(occ[rows, , drop = FALSE]),
      n_islands = rowSums(occ[rows, , drop = FALSE]),
      stringsAsFactors = FALSE
    )
    tips <- tips[order(tip_pos), , drop = FALSE]
    rownames(tips) <- NULL
    clades[[k]] <- list(
      colonisation_time = total_time - rec$col_time,
      crown_age = if (length(rows) > 1) total_time - rec$crown_time else NA,
      branching_times = c(total_time,
                          total_time - rec$col_time,
                          sort(total_time - rec$node_times,
                               decreasing = TRUE)),
      newick = rec$newick,
      nodes = rec$nodes,
      tips = tips,
      mixed = mixed
    )
  }

  # standard DAISIE output, treating the archipelago as one island
  island_spec <- matrix(nrow = 0, ncol = 7)
  if (L > 0) {
    island_spec <- cbind(
      as.character(spec[, .DA_ID]),
      as.character(spec[, .DA_ANC]),
      as.character(spec[, .DA_COLT]),
      daisie_type,
      daisie_code,
      as.character(daisie_btime),
      ana_char[spec[, .DA_ANA] + 1]
    )
  }
  island <- DAISIE_create_island(
    stt_table = stt_table,
    total_time = total_time,
    island_spec = island_spec,
    mainland_n = mainland_n
  )

  # archipelago-specific output
  clade_number <- match(spec[, .DA_CLADE], clade_ids)
  species <- data.frame(
    species = spec[, .DA_ID],
    mainland_ancestor = spec[, .DA_ANC],
    colonisation_time = total_time - spec[, .DA_COLT],
    type = type_char[spec[, .DA_TYPE]],
    anagenetic_origin = ana_char[spec[, .DA_ANA] + 1],
    range = DAISIE_A_range_string(occ),
    range_code = DAISIE_A_range_code(occ),
    n_islands = rowSums(occ),
    clade = clade_number,
    stringsAsFactors = FALSE
  )
  if (L > 0) {
    species <- species[order(species$clade, species$species), , drop = FALSE]
  }
  rownames(species) <- NULL
  island$stt_islands <- stt_islands
  island$species <- species
  island$clades <- clades
  if (length(clades) > 0 && any(vapply(clades, function(x) x$mixed,
                                       logical(1)))) {
    island$mixed_clade <- TRUE
  }
  return(island)
}

#' Reconstruct the phylogeny of one colonisation event
#'
#' @description Prunes the extinct lineages from the full genealogy of one
#'   colonisation event and returns the reconstructed tree, together with the
#'   DAISIE branch codes and branching times of its tips.
#'
#' @param gen Genealogy matrix (closed, i.e. with extant tips).
#' @param clade Clade (colonisation event) ID.
#' @param total_time Simulation time.
#'
#' @return A list with \code{col_time} (forward colonisation time),
#'   \code{crown_time}, \code{node_times} (forward times of the internal
#'   nodes of the reconstructed tree), \code{tips} (genealogy node IDs of the
#'   extant tips), \code{code} and \code{daisie_btime} (forward) per tip,
#'   and \code{newick}.
#' @keywords internal
DAISIE_A_reconstruct_clade <- function(gen, clade, total_time) {
  nn <- nrow(gen)
  parent <- gen[, "parent"]
  time <- gen[, "time"]
  kind <- gen[, "kind"]
  side <- gen[, "side"]
  in_clade <- gen[, "clade"] == clade
  root <- which(in_clade & kind == 1)
  testit::assert(length(root) == 1)
  tips <- which(in_clade & kind == 5)
  marked <- rep(FALSE, nn)
  for (u in tips) {
    w <- u
    while (w > 0 && !marked[w]) {
      marked[w] <- TRUE
      w <- parent[w]
    }
  }
  n_marked_children <- rep(0, nn)
  for (u in which(marked & kind != 1)) {
    n_marked_children[parent[u]] <- n_marked_children[parent[u]] + 1
  }
  kept <- rep(FALSE, nn)
  kept[tips] <- TRUE
  kept[which(kind == 2 & n_marked_children == 2)] <- TRUE
  kept_nodes <- which(kept)
  rec_parent <- rep(NA_integer_, nn)
  rec_side <- rep(NA_real_, nn)
  for (u in kept_nodes) {
    w <- u
    repeat {
      p <- parent[w]
      if (kind[p] == 1) {
        rec_parent[u] <- 0L
        break
      }
      if (kept[p]) {
        rec_parent[u] <- p
        rec_side[u] <- side[w]
        break
      }
      w <- p
    }
  }
  crown <- kept_nodes[rec_parent[kept_nodes] == 0]
  testit::assert(length(crown) == 1)
  # DAISIE branch codes and branching times (the time of the most recent node
  # at which the path to the tip took the B branch; the colonisation time if
  # there is none)
  code <- rep(NA_character_, nn)
  last_b <- rep(NA_real_, nn)
  code[crown] <- ""
  last_b[crown] <- time[root]
  for (u in kept_nodes) {
    if (u == crown) next
    p <- rec_parent[u]
    code[u] <- paste0(code[p], ifelse(rec_side[u] == 1, "A", "B"))
    last_b[u] <- if (rec_side[u] == 2) time[p] else last_b[p]
  }
  internal <- kept_nodes[kind[kept_nodes] == 2]
  # Newick representation (tip labels are species IDs, root edge is the stem)
  children <- split(kept_nodes[kept_nodes != crown],
                    rec_parent[kept_nodes[kept_nodes != crown]])
  newick_rec <- function(u) {
    bl <- if (u == crown) time[u] - time[root] else time[u] - time[rec_parent[u]]
    if (kind[u] == 5) {
      return(paste0(gen[u, "species"], ":", format(bl, digits = 15)))
    }
    ch <- children[[as.character(u)]]
    ch <- ch[order(rec_side[ch])]
    paste0("(", newick_rec(ch[1]), ",", newick_rec(ch[2]), "):",
           format(bl, digits = 15))
  }
  newick <- paste0(newick_rec(crown), ";")
  if (length(tips) == 1) {
    newick <- paste0("(", newick_rec(crown), ");")
  }
  # node table of the reconstructed tree (crown first, parents before
  # children)
  new_id <- match(kept_nodes, kept_nodes)
  nodes <- data.frame(
    node = new_id,
    parent = ifelse(rec_parent[kept_nodes] == 0, 0L,
                    match(rec_parent[kept_nodes], kept_nodes)),
    age = total_time - time[kept_nodes],
    species = ifelse(kind[kept_nodes] == 5, gen[kept_nodes, "species"], NA),
    side = ifelse(rec_parent[kept_nodes] == 0, NA, rec_side[kept_nodes])
  )
  return(list(
    col_time = time[root],
    crown_time = time[crown],
    node_times = time[internal],
    tips = tips,
    code = code[tips],
    daisie_btime = last_b[tips],
    newick = newick,
    nodes = nodes
  ))
}
