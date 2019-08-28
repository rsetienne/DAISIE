context("create_trait_state_params")

test_that("basic use", {
  expect_silent(
    create_trait_state_params(
         trans_rate = 0.5,
         immig_rate2 = 0.1,
         ext_rate2 = 0.2,
         ana_rate2 = 0.3,
         clado_rate2 = 0.4,
         trans_rate2 = 0.5,
         M2 = 1000
    )
  )

})

test_that("abuse", {
  expect_error(
    create_trait_state_params(
      trans_rate = "nonsense",
      immig_rate2 = 0.1,
      ext_rate2 = 0.2,
      ana_rate2 = 0.3,
      clado_rate2 = 0.4,
      trans_rate2 = 0.5,
      M2 = 1000
    ),
    "trans_rate is not of class 'double'"
  )
  
  expect_error(
    create_trait_state_params(
      trans_rate = 0.2,
      immig_rate2 = "nonsense",
      ext_rate2 = 0.2,
      ana_rate2 = 0.3,
      clado_rate2 = 0.4,
      trans_rate2 = 0.5,
      M2 = 1000
    ),
    "immig_rate2 is not of class 'double'"
  )
  
  expect_error(
    create_trait_state_params(
      trans_rate = 0.2,
      immig_rate2 = 0.1,
      ext_rate2 = "nonsense",
      ana_rate2 = 0.3,
      clado_rate2 = 0.4,
      trans_rate2 = 0.5,
      M2 = 1000
    ),
    "ext_rate2 is not of class 'double'"
  )
  
  expect_error(
    create_trait_state_params(
      trans_rate = 0.2,
      immig_rate2 = 0.1,
      ext_rate2 = 0.2,
      ana_rate2 = "nonsense",
      clado_rate2 = 0.4,
      trans_rate2 = 0.5,
      M2 = 1000
    ),
    "ana_rate2 is not of class 'double'"
  )
  
  expect_error(
    create_trait_state_params(
      trans_rate = 0.2,
      immig_rate2 = 0.1,
      ext_rate2 = 0.2,
      ana_rate2 = 0.3,
      clado_rate2 = "nonsense",
      trans_rate2 = 0.5,
      M2 = 1000
    ),
    "clado_rate2 is not of class 'double'"
  )
  
  expect_error(
    create_trait_state_params(
      trans_rate = 0.2,
      immig_rate2 = 0.1,
      ext_rate2 = 0.2,
      ana_rate2 = 0.3,
      clado_rate2 = 0.4,
      trans_rate2 = "nonsense",
      M2 = 1000
    ),
    "trans_rate2 is not of class 'double'"
  )
  
  expect_error(
    create_trait_state_params(
      trans_rate = 0.2,
      immig_rate2 = 0.1,
      ext_rate2 = 0.2,
      ana_rate2 = 0.3,
      clado_rate2 = 0.4,
      trans_rate2 = 0.5,
      M2 = "nonsense"
    ),
    "M2 is not of class 'numeric'"
  )
})
