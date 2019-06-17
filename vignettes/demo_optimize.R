## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
 collapse = TRUE,
 comment = "#>"
)

## ------------------------------------------------------------------------
library(DAISIE)

## ------------------------------------------------------------------------
data(Galapagos_datatable)
knitr::kable(Galapagos_datatable)

## ----fig.width=7, fig.height=7-------------------------------------------
DAISIE::DAISIE_plot_island(Galapagos_datatable, island_age = 4)

## ------------------------------------------------------------------------
data(Galapagos_datatable) 
      
Galapagos_datalist <- DAISIE_dataprep( 
  datatable = Galapagos_datatable, 
  island_age = 4, 
  M = 1000
)

## ------------------------------------------------------------------------
data(Galapagos_datatable) 

Galapagos_datalist_2types <- DAISIE_dataprep( 
  datatable = Galapagos_datatable, 
  island_age = 4, 
  M = 1000, 
  number_clade_types = 2, 
  list_type2_clades = "Finches", 
  prop_type2_pool = 0.163
)

