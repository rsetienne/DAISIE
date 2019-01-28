#' Determines if object x are rates
#' @param x object to be determined to be rates
#' @return TRUE if object x is a list of rates
are_rates <- function(x) {
  if (!"immig_rate" %in% names(x)) return(FALSE)
  if (!"ext_rate" %in% names(x)) return(FALSE)
  if (!"ana_rate" %in% names(x)) return(FALSE)
  if (!"clado_rate" %in% names(x)) return(FALSE)
  if (!"ext_rate_max" %in% names(x)) return(FALSE)
  if (!"immig_rate_max" %in% names(x)) return(FALSE)
  if (!"clado_rate_max" %in% names(x)) return(FALSE)
  if (x$immig_rate < 0.0) return(FALSE)
  if (x$ext_rate < 0.0) return(FALSE)
  if (x$ana_rate < 0.0) return(FALSE)
  if (x$clado_rate < 0.0) return(FALSE)
  if (x$ext_rate_max < 0.0 || x$ext_rate_max < x$ext_rate) return(FALSE)
  if (x$immig_rate_max < 0.0 || x$immig_rate_max < x$immig_rate) return(FALSE)
  if (x$clado_rate_max < 0.0 || x$clado_rate_max < x$clado_rate) return(FALSE)
  TRUE
}