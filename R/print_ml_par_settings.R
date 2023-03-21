#' Print optimisation settings
#'
#' @inheritParams default_params_doc
#'
#' @return Nothing. Prints a `message()` to the console with the parameters
#'   that are to be optimized, fixed, and shifted.
#' @noRd
print_ml_par_settings <- function(namepars,
                                  idparsopt,
                                  idparsfix,
                                  idparsnoshift,
                                  all_no_shift) {
  if (length(namepars[idparsopt]) == 0) {
    optstr <- "nothing"
  } else {
    optstr <- namepars[idparsopt]
  }

  message("You are optimizing ", paste(optstr, collapse = " "))
  if (length(namepars[idparsfix]) == 0) {
    fixstr <- "nothing"
  } else {
    fixstr <- namepars[idparsfix]
  }
  message("You are fixing ", paste(fixstr, collapse = " "))

  if (sum(idparsnoshift %in% (all_no_shift)) != 5) {
    noshiftstring <- namepars[idparsnoshift]
    message("You are not shifting", paste(noshiftstring, collapse = " "))
  }
}
