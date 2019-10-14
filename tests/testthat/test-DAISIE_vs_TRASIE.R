context("DAISIE vs TRASEI")
test_that("TRASIE can get the same result as DAISIE, for IW model", {

  n_mainland_state1 = 1000
  n_mainland_state2 = 1000
  island_age <- 1
  clado_rate <- 0.1 # cladogenesis rate
  ext_rate <- 0.2 # extinction rate (not used)
  island_wide_cap <- 100  # island_wide carrying capacity
  imm_rate <- 0.01 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  trans_rate12 <- 0 # transition rate from state1 to state2
  trans_rate21 <- 0 # transition rate from state2 to state1
  pars = c(clado_rate, ext_rate, island_wide_cap, imm_rate, ana_rate)
  pars_empty = c(0, 0, island_wide_cap, 0, 0)
  Tpars_empty = create_trait_state_params(trans_rate = 0,
                                    immig_rate2 = 0,
                                    ext_rate2 = 0,
                                    ana_rate2 = 0,
                                    clado_rate2 = 0,
                                    trans_rate2 = 0,
                                    M2 = 0)
  Tpars = create_trait_state_params(trans_rate = trans_rate12,
                                          immig_rate2 = imm_rate,
                                          ext_rate2 = ext_rate,
                                          ana_rate2 = ana_rate,
                                          clado_rate2 = clado_rate,
                                          trans_rate2 = trans_rate21,
                                          M2 = n_mainland_state2)
  set.seed(1)
  DAISIE_output <- DAISIE_sim(time = island_age,
                              M = n_mainland_state1,
                              pars = pars,
                              Tpars = NULL,
                              replicates=10,
                              divdepmodel = "IW")
  set.seed(1)
  TRASIE_output1 <- DAISIE_sim(time = island_age,
                               M = n_mainland_state1,
                               pars = pars,
                               Tpars = Tpars_empty,
                               replicates=10,
                               divdepmodel = "IW")
  set.seed(1)
  TRASIE_output2 <- DAISIE_sim(time = island_age,
                               M = 0,
                               pars = pars_empty,
                               Tpars = Tpars,
                               replicates=10,
                               divdepmodel = "IW")
  for(i in 1:length(DAISIE_output)){
    for(j in 2:length(DAISIE_output[[i]])){
      expect_identical(object = TRASIE_output1[[i]][[j]],expected = TRASIE_output2[[i]][[j]])
    }
    expect_identical(object = TRASIE_output1[[i]][[1]]$island_age,expected = TRASIE_output2[[i]][[1]]$island_age)
    expect_identical(object = TRASIE_output1[[i]][[1]]$not_present,expected = TRASIE_output2[[i]][[1]]$not_present)
    expect_identical(object = TRASIE_output1[[i]][[1]]$brts_table,expected = TRASIE_output2[[i]][[1]]$brts_table)
    expect_identical(object = TRASIE_output1[[i]][[1]]$stt_all[,1],expected = TRASIE_output2[[i]][[1]]$stt_all[,1])
    # expect_identical(object = TRASIE_output1[[i]][[1]]$stt_all[1:26,2:4],
    #                  expected = TRASIE_output2[[i]][[1]]$stt_all[1:26,5:7])
  }
  for(i in 1:length(DAISIE_output)){
    for(j in 2:length(DAISIE_output[[i]])){
      expect_identical(object = DAISIE_output[[i]][[j]],expected = TRASIE_output1[[i]][[j]])
    }
    expect_identical(object = DAISIE_output[[i]][[1]]$island_age,expected = TRASIE_output1[[i]][[1]]$island_age)
    expect_identical(object = DAISIE_output[[i]][[1]]$not_present,expected = TRASIE_output1[[i]][[1]]$not_present)
    expect_identical(object = DAISIE_output[[i]][[1]]$brts_table,expected = TRASIE_output1[[i]][[1]]$brts_table)
    expect_identical(object = DAISIE_output[[i]][[1]]$stt_all[,1],expected = TRASIE_output1[[i]][[1]]$stt_all[,1])
    expect_identical(object = DAISIE_output[[i]][[1]]$stt_all,
                     expected = TRASIE_output1[[i]][[1]]$stt_all[,1:4])
  }
})
test_that("TRASIE can get the same result as DAISIE, for CS model", {

  n_mainland_state1 = 50
  n_mainland_state2 = 50
  island_age <- 1
  clado_rate <- 0.1 # cladogenesis rate
  ext_rate <- 0.2 # extinction rate (not used)
  clade_spec_cap <- 10  # clade-level carrying capacity
  imm_rate <- 1 # immigration rate
  ana_rate <- 0.1 # anagenesis rate
  trans_rate12 <- 0 # transition rate from state1 to state2
  trans_rate21 <- 0 # transition rate from state2 to state1
  pars = c(clado_rate, ext_rate, clade_spec_cap, imm_rate, ana_rate)
  pars_empty = c(0, 0, clade_spec_cap, 0, 0)
  Tpars_check = create_trait_state_params(trans_rate = trans_rate12,
                                    immig_rate2 = imm_rate,
                                    ext_rate2 = ext_rate,
                                    ana_rate2 = ana_rate,
                                    clado_rate2 = clado_rate,
                                    trans_rate2 = trans_rate21,
                                    M2 = n_mainland_state2)
  Tpars_balance = create_trait_state_params(trans_rate = trans_rate12,
                                            immig_rate2 = imm_rate,
                                            ext_rate2 = ext_rate,
                                            ana_rate2 = ana_rate,
                                            clado_rate2 = clado_rate,
                                            trans_rate2 = trans_rate21,
                                            M2 = 25)
  set.seed(1)
  DAISIE_output <- DAISIE_sim(time = island_age,
                              M = n_mainland_state1,
                              pars = pars,
                              Tpars = NULL,
                              replicates=5,
                              divdepmodel = "CS")
  set.seed(1)
  TRASIE_output1 <- DAISIE_sim(time = island_age,
                               M = 0,
                               pars = pars_empty,
                               Tpars = Tpars_check,
                               replicates=5,
                               divdepmodel = "CS")
  set.seed(1)
  TRASIE_output2 <- DAISIE_sim(time = island_age,
                               M = 25,
                               pars = pars,
                               Tpars = Tpars_balance,
                               replicates=5,
                               divdepmodel = "CS")
  for(i in 1:length(DAISIE_output)){
    for(j in 2:length(DAISIE_output[[i]])){
      expect_identical(object = TRASIE_output1[[i]][[j]],expected = TRASIE_output2[[i]][[j]])
    }
    expect_identical(object = TRASIE_output1[[i]][[1]]$island_age,expected = TRASIE_output2[[i]][[1]]$island_age)
    expect_identical(object = TRASIE_output1[[i]][[1]]$not_present,expected = TRASIE_output2[[i]][[1]]$not_present)
    expect_identical(object = TRASIE_output1[[i]][[1]]$brts_table,expected = TRASIE_output2[[i]][[1]]$brts_table)
    expect_identical(object = TRASIE_output1[[i]][[1]]$stt_all[,1],expected = TRASIE_output2[[i]][[1]]$stt_all[,1])
    expect_identical(object = TRASIE_output1[[i]][[1]]$stt_all[,"present"],
                     expected = TRASIE_output2[[i]][[1]]$stt_all[,"present"])
  }
  for(i in 1:length(DAISIE_output)){
    for(j in 2:length(DAISIE_output[[i]])){
      expect_identical(object = DAISIE_output[[i]][[j]],expected = TRASIE_output1[[i]][[j]])
    }
    expect_identical(object = DAISIE_output[[i]][[1]]$island_age,expected = TRASIE_output1[[i]][[1]]$island_age)
    expect_identical(object = DAISIE_output[[i]][[1]]$not_present,expected = TRASIE_output1[[i]][[1]]$not_present)
    expect_identical(object = DAISIE_output[[i]][[1]]$brts_table,expected = TRASIE_output1[[i]][[1]]$brts_table)
    expect_identical(object = DAISIE_output[[i]][[1]]$stt_all[,1],expected = TRASIE_output1[[i]][[1]]$stt_all[,1])
    expect_identical(object = DAISIE_output[[i]][[1]]$stt_all[,"present"],
                     expected = TRASIE_output1[[i]][[1]]$stt_all[,"present"])
  }
})
