context("update_rates")

test_that("update rates use", {

  #testit::assert(is.matrix(c()))
  # Does not give errors. One day, it can be checked to be silent  
  set.seed(42)
  update_rates(
    timeval = 0, 
    totaltime = 1, 
    gam = 0.009, 
    mu = 2.0, 
    laa = 1.0, 
    lac = 2.5, 
    Apars = create_area_params(
      max_area = 1.0, 
      proportional_peak_t = 0.5, 
      peak_sharpness = 1.0, 
      total_island_age = 1.0
    ), 
    Epars = c(0.5, 10.0),
    island_ontogeny = "quadratic", 
    extcutoff = 1000.0, 
    K = 3, 
    island_spec = c(), 
    mainland_n = 1, 
    t_hor = 0.5
    )
  are_rates
})

test_that("update_rates classic behavior", {
  
  #testit::assert(is.matrix(c()))
  # Does not give errors. One day, it can be checked to be silent  
  set.seed(42)
  update_rates(
    timeval = 0, 
    totaltime = 1, 
    gam = 0.009, 
    mu = 2.0, 
    laa = 1.0, 
    lac = 2.5, 
    Apars = NULL, 
    Epars = NULL,
    island_ontogeny = NULL, 
    extcutoff = 1000.0, 
    K = 3, 
    island_spec = c(), 
    mainland_n = 1, 
    t_hor = 0.5
  )
  
})
