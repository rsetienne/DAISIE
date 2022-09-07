#' @title Count number of species in DAISIE datalist or simulated data.
#'
#' @description Calculates various island diversity metrics from island datasets.
#'
#' @inheritParams default_params_doc
#' @author Luis Valente
#' @seealso \code{\link{DAISIE_dataprep}},
#' \code{\link{DAISIE_plot_island}}
#'
#' @return The output is a list containing the following items:
#' \item{clade_sizes_sorted}{ List showing the total number of species in each
#' island clade (including missing species). Each item [[i]] on the list
#' gives the sizes of all clades for a single island. If option
#' sort_clade_sizes = T,
#' the clade sizes for are
#' sorted by increasing number of species. If option sort_clade_sizes = F
#' the clade sizes are given in the same order as in the input datalist.}
#' \item{size_largest_clade}{ The total number of species in the largest
#' island clade
#' for each island.}
#' \item{mean_clade_size}{ Mean clade size (average of all island clades)}
#' \item{number_colonisations}{ The total number of colonisations (clades) on
#' each island.}
#' \item{total_number_species}{ The total number of species on each island. These
#' are the extant species at present, including missing species; in case of
#' simulations, this is the number of species present on the island at the
#' end of the
#' simulation.}
#' @examples
#' # Run function with clade sizes in ascending order
#' data("NewZealand_birds_datalist")
#' DAISIE_count_species(NewZealand_birds_datalist)
#'
#' # Run function with clade sizes in the order they appear in the input data
#' DAISIE_count_species(NewZealand_birds_datalist,
#'                        sort_clade_sizes = TRUE)
#' @export
DAISIE_count_species<- function(islands, sort_clade_sizes = TRUE){

if(length(grep("not_present",islands))==1) {islands<-list(islands)}

replicates<-length(islands)
time<-islands[[1]][[1]]$island_age


###### Calculate overall species richness and
## colonization statistics across all islands
number_colonists<-c()
number_species<-c()
size_largest_clade<-c()
mean_clade_size<-c()
clade_sizes<-c()

for(i in 1:replicates){
  the_island<-islands[[i]]
  ncols<-length(the_island)-1
  number_colonists<-append(number_colonists,ncols)

  if(ncols==0){
    number_species<-append(number_species,0)
    size_largest_clade<-append(size_largest_clade,0)
    mean_clade_size<-append(mean_clade_size,0)
    clade_sizes<-append(clade_sizes,0)
     }

  if(ncols>0){
    btimes<-c()
    miss_specs<-c()

    for(o in 2:length(the_island)){
      btimes<-append(btimes,length(the_island[[o]]$branching_times)-1)
      miss_specs<-append(miss_specs,the_island[[o]]$missing_species)
    }
    number_species<-append(number_species,sum(btimes,miss_specs))}
    clade_sizes_dist<-btimes+miss_specs
    size_largest_clade<-append(size_largest_clade,max(clade_sizes_dist))
    mean_clade_size<-append(mean_clade_size,round(mean(clade_sizes_dist),2))
    if(sort_clade_sizes==T) {clade_sizes[i]<-list(sort(clade_sizes_dist))} else{
    clade_sizes[i]<-list(clade_sizes_dist)}

}

overall_results<-list(clade_sizes_sorted=clade_sizes,
                      size_largest_clade=size_largest_clade,
                      mean_clade_size=mean_clade_size,
                      number_colonisations=number_colonists,
                      total_number_species=number_species
                      )

if(sort_clade_sizes==F) {names(overall_results)[1]<-"clade_sizes"}

return(overall_results)
}
