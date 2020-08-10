## code to prepare `frog_M1_ML` dataset goes here

data(frogs_datalist)

frog_M1_ML <- DAISIE_ML(
  datalist = frogs_datalist,
  initparsopt = c(0.18,0.03,0.0006,2),
  idparsopt = c(1,2,4,5),
  ddmodel = 0,
  parsfix = Inf,
  idparsfix = 3
)

saveRDS(frog_M1_ML, file = "vignettes/demo_CSvsIW_files/frog_M1_ML.rds")
