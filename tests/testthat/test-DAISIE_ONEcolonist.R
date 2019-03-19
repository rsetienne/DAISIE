context("DAISIE_ONEcolonist")

test_that("use", {

  sim_time <- 10

  # We need to create 'stt_table' and 'island_spec' 
  if (1 == 2) {
    # Run the code if you need to recreate 'stt_table' and 'island_spec', 
    # Add a print in DAISIE_sim_core_1_4 right before calling 
    # 'DAISIE_ONEcolonist'
    set.seed(42)
    n_mainland_species <- 1
    clado_rate <- 1.0 
    ext_rate <- 0.1
    carr_cap <- 4
    imm_rate <- 1.0
    ana_rate <- 1.0
    DAISIE::DAISIE_sim_core_1_4(
      time = sim_time,
      mainland_n = n_mainland_species,
      pars = c(clado_rate, ext_rate, carr_cap, imm_rate, ana_rate)
    )
  }
  
  #       Time nI nA nC
  #  [1,] 10.0000000  0  0  0
  #  [2,]  9.8016632  1  0  0
  #  [3,]  9.7869740  0  0  2
  #  [4,]  9.6022772  1  0  2
  #  [5,]  9.4117448  1  0  2
  #  [6,]  8.8270858  0  1  2
  #  [7,]  7.9951318  0  0  4
  #  [8,]  7.8522784  0  1  2
  #  [9,]  6.9258096  1  1  2
  # [10,]  6.6612814  0  2  2
  # [11,]  4.9967010  0  3  0
  # [12,]  4.7943237  0  2  2
  # [13,]  4.0461558  0  3  0
  # [14,]  3.8633278  0  2  2
  # [15,]  0.6456057  0  1  2
  # [16,]  0.2628436  0  2  0
  # [17,]  0.0000000  0  2  0  
  stt_table <- data.frame(
    Time = c(
      10,
      9.8016632,
      9.786974,
      9.6022772,
      9.4117448,
      8.8270858,
      7.9951318,
      7.8522784,
      6.9258096,
      6.6612814,
      4.996701,
      4.7943237,
      4.0461558,
      3.8633278,
      0.6456057,
      0.2628436,
      0
    ),
    nI = c(0,1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
    nA = c(0,0, 0, 0, 0, 1, 0, 1, 1, 2, 3, 2, 3, 2, 1, 2, 2),
    nC = c(0,0, 2, 2, 2, 2, 4, 2, 2, 2, 0, 2, 0, 2, 2, 0, 0)
  )
  
  # island spec is a matrix of strings converted to a data frame.
  # Obtained by using code above
  #
  #     Species Mainland Ancestor Colonisation time (BP) Species type branch_code branching time (BP)
  # [1,] "7"     "1"               "6.92580955162582"     "A"          NA          NA                 
  # [2,] "11"    "1"               "9.41174479159888"     "A"          NA          NA                 
  #      Anagenetic_origin
  # [1,] "Immig_parent"   
  # [2,] "Clado_extinct" 
  island_spec <- matrix(nrow = 2, ncol = 7, data = "x")
  island_spec[,1] <- c("7", "11")
  island_spec[,2] <- c("1", "1")
  island_spec[,3] <- c("6.92580955162582", "9.41174479159888")
  island_spec[,4] <- c("A", "A")
  island_spec[,5] <- c(NA, NA)
  island_spec[,6] <- c(NA, NA)
  island_spec[,7] <- c("Immig_parent", "Clado_extinct")
  colnames(island_spec) <- c(
    "Species",
    "Mainland Ancestor",
    "Colonisation time (BP)",
    "Species type",
    "branch_code",
    "branching time (BP)",
    "Anagenetic_origin"
  )
  result <- DAISIE:::DAISIE_ONEcolonist(
    time = sim_time,
    island_spec = island_spec,
    stt_table = stt_table
  )
  # result:
  #
  # $stt_table [RJCB: same as put in]
  #          Time nI nA nC
  # 1  10.0000000  0  0  0
  # 2   9.8016632  1  0  0
  # 3   9.7869740  0  0  2
  # 4   9.6022772  1  0  2
  # 5   9.4117448  1  0  2
  # 6   8.8270858  0  1  2
  # 7   7.9951318  0  0  4
  # 8   7.8522784  0  1  2
  # 9   6.9258096  1  1  2
  # 10  6.6612814  0  2  2
  # 11  4.9967010  0  3  0
  # 12  4.7943237  0  2  2
  # 13  4.0461558  0  3  0
  # 14  3.8633278  0  2  2
  # 15  0.6456057  0  1  2
  # 16  0.2628436  0  2  0
  # 17  0.0000000  0  2  0
  # 
  # $branching_times
  # [1] 10.000000  9.411745
  # 
  # $stac
  # [1] 3
  # 
  # $missing_species
  # [1] 0
  # 
  # $other_clades_same_ancestor
  # $other_clades_same_ancestor[[1]]
  # $other_clades_same_ancestor[[1]]$brts_miss
  # [1] 6.92581
  # 
  # $other_clades_same_ancestor[[1]]$species_type
  # [1] "A"  
  
  expect_equal(result$stt_table, stt_table)
  expect_true(
    all.equal(
      result$branching_times, 
      c(10.000000, 9.411745), 
      tolerance = 1.0e-7 # OK, caused by string conversion
    )
  )
  expect_equal(result$stac, 3)
  expect_equal(result$missing_species, 0)
  expect_equal(length(result$other_clades_same_ancestor), 1)
  expect_true(
    all.equal(
      result$other_clades_same_ancestor[[1]]$brts_miss, 
      6.92581,
      tolerance = 1.0e-7 # OK, caused by string conversion
    )
  )
  expect_equal(
    result$other_clades_same_ancestor[[1]]$species_type,
    "A"
  )
  
})
