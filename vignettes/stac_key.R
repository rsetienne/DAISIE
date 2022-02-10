## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
library(readr)

## ----load table, echo = FALSE, message = FALSE--------------------------------
data(stac_table, package = "DAISIE")

## ----table, echo=FALSE--------------------------------------------------------
knitr::kable(stac_table)

