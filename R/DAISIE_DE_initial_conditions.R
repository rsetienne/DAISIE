# initial conditions system of equation interval 2
#' @keywords internal
DAISIE_DE_get_initial_conditions2 <- function(status,
                                    brts,
                                    sampling_fraction) {

  if (status == 6 || status == 2 && length(brts) > 2) {
    initial_conditions2 <- c(
      DE = sampling_fraction,
      DM3 = 0,
      E = 1 - sampling_fraction,
      DA3 = 1
    )
  } else if (status == 3 && length(brts) > 2) {
    initial_conditions2 <- c(
      DE = sampling_fraction,
      DM3 = 1,
      E = 1 - sampling_fraction,
      DA3 = 0
    )
  } else if (status == 9 || status == 2 && length(brts) == 2) {
    initial_conditions2 <- c(
      DE = sampling_fraction,
      DM2 = 0,
      DM3 = 0,
      E = 1 - sampling_fraction,
      DA3 = 1
    )

  } else if (status == 3 && length(brts) == 2) {
    initial_conditions2 <- c(
      DE = sampling_fraction,
      DM2 = 0,
      DM3 = 1,
      E = 1 - sampling_fraction,
      DA3 = 0
    )

  } else if (status %in% c(4, 8)) {
    initial_conditions2 <- c(
      DM2 = 1,
      E = 0
    )

  } else {
    stop("Unsupported status or brts length in get_initial_conditions2()")
  }

  return(initial_conditions2)
}



# initial conditions system of equation interval 3
#' @keywords internal
DAISIE_DE_get_initial_conditions3 <- function(status,
                                    solution = NULL,
                                    sampling_fraction = NA) {

  if (status == 1) {
    initial_conditions3 <- c(
      DM1 = 0,
      DM2 = 1,
      E = 0,
      DA2 = 0
    )

  } else if (status == 6) {

    initial_conditions3 <- c(
      DE  = solution[, "DE"][[2]],
      DM1 = solution[, "DM1"][[2]],
      DM2 = solution[, "DM2"][[2]],
      DM3 = solution[, "DM3"][[2]],
      DA2 = solution[, "DA2"][[2]],
      DA3 = solution[, "DA3"][[2]],
      E   = solution[, "E"][[2]]
    )

  } else if (status == 5) {
    initial_conditions3 <- c(
      DE = sampling_fraction,
      DM1 = 0,
      DM2 = 0,
      DM3 = 0,
      E = 1 - sampling_fraction,
      DA2 = 0,
      DA3 = 1
    )

  } else if (status == 8 || status == 9) {
    initial_conditions3 <- c(
      DE = solution[, "DE"][[2]],
      DM1 = 0,
      DM2 = solution[, "DM2"][[2]],
      DM3 = solution[, "DM3"][[2]],
      E   = solution[, "E"][[2]],
      DA2 = 0,
      DA3 = solution[, "DA3"][[2]]
    )

  } else {
    stop("Unsupported status in get_initial_conditions3()")
  }

  return(initial_conditions3)
}



# initial conditions system of equation interval 4
#' @keywords internal
DAISIE_DE_get_initial_conditions4 <- function(status,
                                    solution,
                                    parameter) {

  if (status %in% c(2, 3, 4)) {
    initial_conditions4 <- c(
      DA1 = parameter[4] * solution[, "DM2"][[2]],
      DM1 = parameter[4] * solution[, "DM2"][[2]],
      E   = solution[, "E"][[2]]
    )

  } else if (status %in% c(1, 5, 6, 8, 9)) {
    initial_conditions4 <- c(
      DA1 = solution[, "DA2"][[2]],
      DM1 = solution[, "DM1"][[2]],
      E   = solution[, "E"][[2]]
    )

  } else {
    stop("Unsupported status in get_initial_conditions4()")
  }

  return(initial_conditions4)
}

