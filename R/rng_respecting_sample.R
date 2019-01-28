#' Sampling in which zero probabilities are removed
#' @inheritParams base::sample
#' @return a vector of length size with elements drawn 
#'   from either x or from the integers 1:x.
#' @examples 
#'   # Number of draws
#'   n <- 1000
#'   
#'   # Do normal sampling
#'   set.seed(42)
#'   draws_1 <- rng_respecting_sample(
#'     1:3, size = n, replace = TRUE, prob = c(1.0, 1.0, 1.0)
#'   )
#'   
#'   # Do a sampling with one element of probabily zero
#'   set.seed(42)
#'   draws_2 <- rng_respecting_sample(
#'     1:4, size = n, replace = TRUE, prob = c(1.0, 1.0, 1.0, 0.0)
#'   )
#'   testit::assert(sum(draws_2 == 4) == 0)
#'   testit::assert(draws_1 == draws_2)
#'   
#'   # Use base sampling will give different results,
#'   # as it results in different RNG values
#'   set.seed(42)
#'   draws_3 <- sample(
#'     1:4, size = n, replace = TRUE, prob = c(1.0, 1.0, 1.0, 0.0)
#'   )
#'   testit::assert(sum(draws_3 == 4) == 0)
#'   testit::assert(!all(draws_1 == draws_3))
#'   
#' @author Richel J.C. Bilderbeek
#' @seealso See \code{\link[base]{sample}} for more details
#' @note thanks to Pedro Neves for finding this feature in base::sample
rng_respecting_sample <- function(x, size, replace, prob) {
  which_non_zero <- prob > 0.0
  non_zero_prob <- prob[which_non_zero]
  non_zero_x <- prob[which_non_zero]
  sample(x = non_zero_x, size = size, replace = replace, prob = non_zero_prob)
}