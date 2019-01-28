#' Extract the sorted branching times, in million years ago.
#' from a data table
#' @param data_table data table
#' @return the sorted branching times, in million years ago
#' @examples 
#'   data(Galapagos_datatable)
#'   brts_mya <- DAISIE_get_brts_mya(data_table = Galapagos_datatable)
#'   testit::assert(length(brts_mya) > 1)
#'   testit::assert(all(brts_mya > 0))
#' @author Richel J.C. Bilderbeek
DAISIE_get_brts_mya <- function(data_table) {
  testit::assert("Branching_times" %in% names(data_table))
  brts_mya <- c()
  for (t in data_table$Branching_times) {
    brts_mya <- as.numeric(c(brts_mya, strsplit(as.character(t),split = ",")[[1]]))
  }
  sort(brts_mya)
}