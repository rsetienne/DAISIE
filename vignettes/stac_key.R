## ----setup, include=FALSE-----------------------------------------------------
library(knitr)

## ----load table, echo = FALSE, message = FALSE--------------------------------
data(stac_table, package = "DAISIE")

## ----table, echo=FALSE--------------------------------------------------------
knitr::kable(stac_table)

