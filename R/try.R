totaltime=4
M=1000
ddmodel = c(1,0,1)
island_type = "oceanic"
nonoceanic = NULL
Apars = NULL
Epars = NULL
Tpars = create_trait_state_params(trans_rate = 0,immig_rate2 = 0,ext_rate2 = 0,ana_rate2 = 0,clado_rate2 = 0,trans_rate2 = 0,M2 = 1000)
keep_final_state = FALSE
island_spec = NULL
replicates = 25
mainland_params = NULL
divdepmodel = 'IW'
nonoceanic = NULL
prop_type2_pool = NA
replicates_apply_type2 = TRUE
sample_freq = 25
plot_sims = FALSE
island_ontogeny = "const" 
Tpars = Tpars
Apars = NULL
Epars = NULL
keep_final_state = FALSE
stored_data = NULL
verbose = TRUE
pars = c(0.1, 0.2, 10, 0.1, 0.2)

set.seed(1);DAISIE_sim_core(time = totaltime,mainland_n = M,pars = pars,island_spec=island_spec)
set.seed(1);DAISIE_sim_core(time = totaltime,mainland_n = M,pars = pars,island_spec=island_spec,Tpars = Tpars)

DAISIE_sim_core(
  time = totaltime,
  mainland_n = M,
  pars = pars,
  ddmodel = ddmodel,
  island_type = island_type,
  nonoceanic = nonoceanic,
  island_ontogeny = island_ontogeny,
  Apars = Apars,
  Epars = Epars,
  Tpars = Tpars,
  keep_final_state = keep_final_state,
  island_spec = island_spec
)

set.seed(1);DAISIE_sim(
  time = totaltime,
  M = M,
  pars = pars,
  replicates = replicates,
  mainland_params = NULL,
  divdepmodel = 'IW',
  island_type = "oceanic",
  nonoceanic = NULL,
  prop_type2_pool = NA,
  replicates_apply_type2 = TRUE,
  sample_freq = 25,
  plot_sims = FALSE,
  island_ontogeny = "const",
  Tpars = Tpars,
  Apars = NULL,
  Epars = NULL,
  keep_final_state = FALSE,
  stored_data = NULL,
  verbose = TRUE)
