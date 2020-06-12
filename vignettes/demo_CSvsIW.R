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

