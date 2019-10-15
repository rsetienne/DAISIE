if (1 == 2) {
  # Issue #81   estimate params using simulated data from DAISIE( Tpars = NULL), use ML_IW and ML_CS separately. Work!
  test_that("general DAISIE_IW", {
    clado_rate <- 0.5 # cladogenesis rate
    ext_rate <- 0.2 # extinction rate
    island_wide_cap <- 20 # clade-level carrying capacity
    imm_rate <- 0.002 # immigration rate
    ana_rate <- 0.1 # anagenesis rate
    pars1 <- c(clado_rate, ext_rate, island_wide_cap, imm_rate, ana_rate)
    set.seed(100);islands_replicates_IW = DAISIE_sim(time=4,
                                                     M=500,
                                                     pars=pars1,
                                                     Tpars = NULL,
                                                     replicates=1,
                                                     divdepmodel = "IW")
    ###  estimate params using ML_CS
    lambda_c <- factor()
    mu <- factor()
    K <- factor()
    gamma <- factor()
    lambda_a <- factor()
    loglik <- factor()
    df <- factor()
    conv <- factor()
    DAISIEML_CS_table <- data.frame(lambda_c,mu,K,gamma,lambda_a
                                    ,loglik,df,conv)
    for (i in 1:length(islands_replicates_IW)){
      DAISIEML_CS <- DAISIE::DAISIE_ML_CS(datalist=islands_replicates_IW[[i]],
                                          initparsopt = c(0.1,0.1,10,0.1,0.1),
                                          ddmodel=0,
                                          idparsopt = 1:5,
                                          parsfix = NULL,
                                          idparsfix = NULL)
      DAISIEML_CS_table <- rbind(DAISIEML_CS_table, DAISIEML_CS)
    }
    # Calculating the likelihood for the initial parameters.
    # The loglikelihood for the initial parameter values is -174.3363
    # Optimizing the likelihood - this may take a while.
    #
    # Maximum likelihood parameter estimates: lambda_c: 0.190015, mu: 0.000000, K: 10.564105, gamma: 0.001000, lambda_a: 1650.628953
    # Maximum loglikelihood: -18.476255

    ###  estimate params using ML_IW
    DAISIEML_IW_table <- data.frame(lambda_c,mu,K,gamma,lambda_a
                                    ,loglik,df,conv)

    for (i in 1:length(islands_replicates_IW)){
      DAISIEML_IW <- DAISIE::DAISIE_ML_IW(datalist=islands_replicates_IW[[i]],
                                          initparsopt = c(0.5,0.2,10,0.001,0.1),
                                          ddmodel=0,
                                          idparsopt = 1:5,
                                          parsfix = NULL,
                                          idparsfix = NULL)
      DAISIEML_IW_table <- rbind(DAISIEML_IW_table, DAISIEML_IW)
    }

    # The loglikelihood for the initial parameter values is -19.5761 ,too lang time to get estimated values, so stop
  })
}
####  test when some of the values are zero    The crash probably caused by terrible initial value with 0s.
if (1 == 2) {
  # Issue #81  not crash now, while cause numerical problem.
  test_that("TRASIE IW crash", {
    pars1 <- c(0,0,Inf,0.01,0)
    Tpars1 <- create_trait_state_params(trans_rate = 0,
                                        immig_rate2 = 0.01,
                                        ext_rate2 = 0,
                                        ana_rate2 = 0,
                                        clado_rate2 = 0,
                                        trans_rate2 = 0,
                                        M2 = 20)
    set.seed(100);islands_replicates_IW = DAISIE_sim(time=10,
                                                     M=20,
                                                     pars=pars1,
                                                     Tpars = Tpars1,
                                                     replicates=1,
                                                     divdepmodel = "IW")
    lambda_c <- factor()
    mu <- factor()
    K <- factor()
    gamma <- factor()
    lambda_a <- factor()
    loglik <- factor()
    df <- factor()
    conv <- factor()
    DAISIEML_IW_table <- data.frame(lambda_c,mu,K,gamma,lambda_a
                                    ,loglik,df,conv)

    for (i in 1:length(islands_replicates_IW)){
      DAISIEML_IW <- DAISIE::DAISIE_ML_IW(datalist=islands_replicates_IW[[i]],
                                          initparsopt = c(0,0,20,0.1,0),
                                          ddmodel=0,
                                          idparsopt = 1:5,
                                          parsfix = NULL,
                                          idparsfix = NULL)
      DAISIEML_IW_table <- rbind(DAISIEML_IW_table, DAISIEML_IW)
    }

    # Calculating the likelihood for the initial parameters.
    # The loglikelihood for the initial parameter values is -42.30258
    # Optimizing the likelihood - this may take a while.
    # There are parameter values used which cause numerical problems.
  })
}

##### test what will happen when estimating crashed simulated data using ML_CS
if (1 == 2) {
  # Issue #81
  test_that("TRASIE IW crash", {
    pars1 <- c(0,0,Inf,0.001,0)
    Tpars1 <- create_trait_state_params(trans_rate = 0,
                                        immig_rate2 = 0.001,
                                        ext_rate2 = 0,
                                        ana_rate2 = 0,
                                        clado_rate2 = 0,
                                        trans_rate2 = 0,
                                        M2 = 500)
    set.seed(100);islands_replicates_IW = DAISIE_sim(time=4,
                                                     M=500,
                                                     pars=pars1,
                                                     Tpars = Tpars1,
                                                     replicates=2,
                                                     divdepmodel = "IW")
    lambda_c <- factor()
    mu <- factor()
    K <- factor()
    gamma <- factor()
    lambda_a <- factor()
    loglik <- factor()
    df <- factor()
    conv <- factor()
    DAISIEML_CS_table <- data.frame(lambda_c,mu,K,gamma,lambda_a
                                    ,loglik,df,conv)

    for (i in 1:length(islands_replicates_IW)){
      DAISIEML_CS <- DAISIE::DAISIE_ML_CS(datalist=islands_replicates_IW[[i]],
                                          initparsopt = c(0.1,0.1,10,0.1,0.1),
                                          ddmodel=0,
                                          idparsopt = 1:5,
                                          parsfix = NULL,
                                          idparsfix = NULL)
      DAISIEML_CS_table <- rbind(DAISIEML_CS_table, DAISIEML_CS)

      ###The loglikelihood for the initial parameter values is -333.8213
      ###Optimizing the likelihood - this may take a while.

      ##Maximum likelihood parameter estimates: lambda_c: 0.000000, mu: 0.000009, K: 24.068738, gamma: 0.000250, lambda_a: 0.000007
      ##Maximum loglikelihood: -9.294080
    }
  })
}
#### test what will happen when no zero values
if (1 == 2) {
  # Issue #81
  test_that("TRASIE IW crash", {
    pars1 <- c(0.5,0.2,20,0.01,0.1)
    Tpars1 <- create_trait_state_params(trans_rate = 0,
                                        immig_rate2 = 0.01,
                                        ext_rate2 = 0.2,
                                        ana_rate2 = 0.1,
                                        clado_rate2 = 0.5,
                                        trans_rate2 = 0,
                                        M2 = 20)
    set.seed(100);islands_replicates_IW = DAISIE_sim(time=4,
                                                     M=20,
                                                     pars=pars1,
                                                     Tpars = Tpars1,
                                                     replicates=2,
                                                     divdepmodel = "IW")
    lambda_c <- factor()
    mu <- factor()
    K <- factor()
    gamma <- factor()
    lambda_a <- factor()
    loglik <- factor()
    df <- factor()
    conv <- factor()
    DAISIEML_IW_table <- data.frame(lambda_c,mu,K,gamma,lambda_a
                                    ,loglik,df,conv)

    for (i in 1:length(islands_replicates_IW)){
      DAISIEML_IW <- DAISIE::DAISIE_ML_IW(datalist=islands_replicates_IW[[i]],
                                          initparsopt = c(0.1,0.1,10,0.1,0.1),
                                          ddmodel=0,
                                          idparsopt = 1:5,
                                          parsfix = NULL,
                                          idparsfix = NULL)
      DAISIEML_IW_table <- rbind(DAISIEML_IW_table, DAISIEML_IW)
    }
  })
}
if (1 == 2) {
  # Issue #81
  test_that("check CS", {
    pars1 <- c(0.5,0.2,20,1,0.2)
    Tpars1 <- create_trait_state_params(trans_rate = 0,
                                        immig_rate2 = 1,
                                        ext_rate2 = 0.2,
                                        ana_rate2 = 0.2,
                                        clado_rate2 = 0.5,
                                        trans_rate2 = 0,
                                        M2 = 10)
    set.seed(100);islands_replicates_CS = DAISIE_sim(time=4,
                                                     M=10,
                                                     pars=pars1,
                                                     Tpars = Tpars1,
                                                     replicates=1,
                                                     divdepmodel = "CS")
    lambda_c <- factor()
    mu <- factor()
    K <- factor()
    gamma <- factor()
    lambda_a <- factor()
    loglik <- factor()
    df <- factor()
    conv <- factor()
    DAISIEML_CS_table <- data.frame(lambda_c,mu,K,gamma,lambda_a
                                    ,loglik,df,conv)

    for (i in 1:length(islands_replicates_CS)){
      DAISIEML_CS <- DAISIE::DAISIE_ML_CS(datalist=islands_replicates_CS[[i]],
                                          initparsopt = c(0.1,0.1,10,0.1,0.1),
                                          ddmodel=0,
                                          idparsopt = 1:5,
                                          parsfix = NULL,
                                          idparsfix = NULL)
      DAISIEML_CS_table <- rbind(DAISIEML_CS_table, DAISIEML_CS)
    }
  })
}
