#' Plot clade age against clade diversity.
#' 
#' Plots clade age against clade diversity for all clades for which
#' colonisation time is known.
#' 
#' R plot showing for each clade in the island object the time of colonisation
#' and the diversity of the clade. Only clades for which colonisation time is
#' known are plotted. Blue - endemic; black - non-endemic. C = number of
#' independent colonisations on island, N = number of species on the island
#' 
#' @param island Island data object. Can be in DAISIE list format (see
#' Galapagos_datalist and DAISIE_data_prep for examples) or in table format
#' (see Galapagos_datatable for an example)
#' @param title Title of the plot
#' @param island_age If island input is in table format, the age of the island
#' must be specified. If island input is in DAISIE list format, this option
#' will override the island age specified in the island list.
#' @return R plot.
#' @author Luis Valente
#' @seealso \code{\link{DAISIE_sim}} \code{\link{DAISIE_plot_island}}
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#' @keywords models
#' @examples
#'  
#' 
#' ### Plot Galapagos age-diversity for Galapagos dataset
#' utils::data(Galapagos_datalist)
#' DAISIE_plot_age_diversity(Galapagos_datalist)
#' 
#' 
#' @export DAISIE_plot_age_diversity
DAISIE_plot_age_diversity = function(
  island,
  title = 'Clade age vs clade diversity',
  island_age = NA) 
{

  if(class(island)=='data.frame'){
    if(is.na(island_age)) {stop('Island age must be specified if the input file is a table')}  
    
    island<-DAISIE_dataprep(island,island_age,1000)}
  
  if(is.na(island_age)){island_age<-island[[1]]$island_age}
  

if(island[[2]]$stac==0){ 
  graphics::plot(NULL,NULL,xlim=c(-island_age,0),ylim=c(0,2))
  graphics::legend("center", paste('Island has no species'),bty="n") 
  }else{
  
island[[1]]<-NULL

stac_vec = unlist(island)[which(names(unlist(island)) == "stac")]
stac_age_known<-which(stac_vec == 2 | stac_vec == 3 | stac_vec == 4 )
stacs_used<-stac_vec[stac_age_known]

btimes<-list()
colonisation_time<-c()
diversity<-c()
  
for (i in 1:length(stac_age_known)) {
  btimes[[i]]<-island[[stac_age_known[i]]]$branching_times[-1]
  colonisation_time[i]<-max(btimes[[i]])
  diversity[i]<-length(btimes[[i]])+island[[stac_age_known[i]]]$missing_species
  }
graphics::plot(-colonisation_time,diversity,xlim=c(-island_age,0),ylab='Diversity',xlab='Colonisation time',col=(stacs_used==2 | stacs_used==3)*4 + (stacs_used==4)*1,pch=(stacs_used==2 | stacs_used==3)*16 + (stacs_used==4)*15,cex=1.2,main=title)
graphics::legend("topleft", paste('C=',length(stac_age_known),' ','N=',sum(diversity),sep=''),bty="n") 

}
}
