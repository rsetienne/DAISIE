## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
 collapse = TRUE,
 comment = "#>"
)

## -----------------------------------------------------------------------------
library(DAISIE)

## -----------------------------------------------------------------------------
data(frogs_datatable, package = "DAISIE")

## ---- echo=FALSE, results= 'asis'---------------------------------------------
knitr::kable(frogs_datatable,caption= 'Eleutherodactylus data table')

## -----------------------------------------------------------------------------
frogs_datalist <- DAISIE_dataprep(
        datatable = frogs_datatable,
        island_age = 30,
        M = 300)

## ---- results='hide'----------------------------------------------------------
frogs_datalist

## ---- echo=TRUE---------------------------------------------------------------
DAISIE_plot_island(frogs_datalist)

## ---- echo=TRUE, fig.align='left'---------------------------------------------
# DAISIE_plot_age_diversity(frogs_datalist)

## -----------------------------------------------------------------------------
data(frogs_datalist)

## ----ML_DAISIE_M1, warning=TRUE, cache= TRUE----------------------------------
DAISIE_ML(
  datalist = frogs_datalist,
  initparsopt = c(0.18,0.03,0.0006,2),
  idparsopt = c(1,2,4,5),
  ddmodel = 0,
  parsfix = Inf,
  idparsfix = 3
)

## ----ML_DAISIE_M2, eval=FALSE-------------------------------------------------
#  DAISIE_ML(
#    datalist = frogs_datalist,
#    initparsopt = c(0.44,0.11,36.44,0.0007,2),
#    idparsopt = c(1,2,3,4,5),
#    ddmodel = 11,
#    parsfix = NULL,
#    idparsfix = NULL
#  )

## ----ML_DAISIE_M3, eval=FALSE-------------------------------------------------
#  DAISIE_ML(
#    datalist = frogs_datalist,
#    initparsopt = c(0.44,0.11,36.44,0.0007),
#    idparsopt = c(1,2,3,4),
#    ddmodel = 11,
#    parsfix = 0,
#    idparsfix = 5
#  )

## ----ML_DAISIE_M4, eval=FALSE-------------------------------------------------
#  DAISIE_ML_IW(
#    datalist = frogs_datalist,
#    initparsopt = c(0.41, 0.17, 131.7, 0.0012, 2),
#    idparsopt = c(1,2,3,4,5),
#    ddmodel = 11,
#    parsfix = NULL,
#    idparsfix = NULL
#  )

## ----ML_DAISIE_M5, eval=FALSE-------------------------------------------------
#  DAISIE_ML_IW(
#    datalist = frogs_datalist,
#    initparsopt = c(0.40, 0.17, 131.83, 0.0012),
#    idparsopt = c(1,2,3,4),
#    ddmodel = 11,
#    parsfix = 0,
#    idparsfix = 5
#  )

## ----DAISIE_sim_CS, results='hide', cache= TRUE-------------------------------
frog_sims_CS<-DAISIE_sim_constant_rate(
  time=30,
  M=300,
  pars=c(0.44,0.11,36.44,0.0007,0),
  divdepmodel = "CS",
  replicates= 100,
  plot_sims = FALSE)

## ----DAISIE_sim_IW, results='hide', cache= TRUE-------------------------------
frog_sims_CS<-DAISIE_sim_constant_rate(
  time=30,
  M=300,
  pars=c(0.40,0.17,131.83,0.0012,0),
  divdepmodel = "IW",
  replicates= 100,
  plot_sims = FALSE)

## ---- echo=TRUE---------------------------------------------------------------
DAISIE_plot_sims(frog_sims_CS)

