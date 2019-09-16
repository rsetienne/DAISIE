test_that("use", {
  
  expect_false(is_simulation_outputs("nonsense"))
  expect_true(is_simulation_outputs(create_test_simulation_outputs()))

})
