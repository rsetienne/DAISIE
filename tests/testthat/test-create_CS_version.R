context("create_CS_version")

test_that("create_CS_version produces correct output", {
  CS_version <- create_CS_version(model = 1,
                                  pick_parameter = NULL,
                                  distribution = NULL,
                                  sd = NULL,
                                  multi_rate_optim_method = NULL)
  expect_equal(CS_version, list(model = 1,
                                pick_parameter = NULL,
                                distribution = NULL,
                                sd = NULL,
                                optimmethod = NULL))

  CS_version <- create_CS_version(model = 2,
                                  pick_parameter = "cladogenesis",
                                  distribution = "gamma",
                                  sd = 1,
                                  multi_rate_optim_method = "optimize")
  expect_equal(CS_version, list(model = 2,
                                pick_parameter = "cladogenesis",
                                distribution = "gamma",
                                sd = 1,
                                optimmethod = "optimize"))

  CS_version <- create_CS_version(model = 3,
                                  pick_parameter = NULL,
                                  distribution = NULL,
                                  sd = NULL,
                                  multi_rate_optim_method = NULL)
  expect_equal(CS_version, list(model = 3,
                                pick_parameter = NULL,
                                distribution = NULL,
                                sd = NULL,
                                optimmethod = NULL))

})


test_that("abuse create_CS_version", {
  expect_error(create_CS_version(model = 4,
                                  pick_parameter = NULL,
                                  distribution = NULL,
                                  sd = NULL,
                                  multi_rate_optim_method = NULL))

  expect_error(create_CS_version(model = 2,
                                  pick_parameter = NULL,
                                  distribution = "gamma",
                                  sd = 1,
                                  multi_rate_optim_method = "optimize"))

  expect_error(create_CS_version(model = 2,
                                 pick_parameter = "cladogenesis",
                                 distribution = NULL,
                                 sd = 1,
                                 multi_rate_optim_method = "optimize"))

  expect_error(create_CS_version(model = 2,
                                 pick_parameter = "cladogenesis",
                                 distribution = "gamma",
                                 sd = NULL,
                                 multi_rate_optim_method = "optimize"))

  expect_error(create_CS_version(model = 2,
                                 pick_parameter = "cladogenesis",
                                 distribution = "gamma",
                                 sd = 1,
                                 multi_rate_optim_method = NULL))
})
