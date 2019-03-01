#' Plot colonisation and branching time of species found on an island dataset.
#' 
#' Produces an image with the times of colonisation, branching times, clade
#' name and clade status for all clades found in a given dataset.
#' 
#' R plot showing for each clade in the island object: time of colonisation,
#' branching times, species status (endemic or non-endemic), total number of
#' species in clade (n), number of species not sampled (m).
#' 
#' @param island Island data object. Can be in DAISIE list format (see
#' Galapagos_datalist and DAISIE_data_prep for examples) or in table format
#' (see Galapagos_datatable for an example)
#' @param island_age If island input is in table format, the age of the island
#' must be specified. If island input is in DAISIE list format, this option
#' will override the island age specified in the island list.
#' @return R plot.
#' @author Luis Valente
#' @seealso \code{\link{DAISIE_sim}} \code{\link{DAISIE_plot_age_diversity}}
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#' @keywords models
#' @examples
#'  
#' 
#' ### Plot Galapagos islands dataset
#' utils::data(Galapagos_datalist)
#' DAISIE_plot_island(Galapagos_datalist)
#' 
#' 
#' @export DAISIE_plot_island
DAISIE_plot_island = function (island,island_age=NA) { 
  
  if(class(island)=='data.frame'){
    if(is.na(island_age)) {stop('Island age must be specified if the input file is a table')}  
    
    island<-DAISIE_dataprep(island,island_age,1000)}
  
  if(is.na(island_age)){island_age<-island[[1]]$island_age}
  
  island[[1]]<-NULL
  
  btimes<-list()
  for (i in 1:length(island)) {
    btimes[[i]]<-island[[i]]$branching_times[-1] 
  }
  island<-island[rev(order(sapply(btimes,"[",1)))]
  
  
  colonists<-c()
  for (i in 1:length(island)) {
    
    if(is.null(island[[i]]$colonist_name)==T){
    island[[i]]$colonist_name<-paste('Col',i,sep='')
    }
    
    colonists[i]<-strsplit(island[[i]]$colonist_name,'_')[[1]][1]
  }
  
  ymax<-max(island_age,island[[1]]$branching_times[2])
  
  
  graphics::plot(NULL,NULL,xlim=c(0,2*length(island)),ylim=c(0,ymax*1.1),xaxt='n',ann=FALSE,yaxs="i")
  
  graphics::abline(h=island_age,lty=2)
  
  position<-seq(2,length(island)*2,by=2)
  
  graphics::axis(1,labels=colonists, at=position, las=2,tick=F,cex.axis=0.7,line=-0.5)
  
  graphics::legend(x=length(island)+length(island)/5,y=ymax,
         c("Non-endemic","Non-endemic MaxAge","Endemic",'Endemic MaxAge'),
         pch=c(22,22,21,21),col=c(1,1,4,4),pt.bg=c(1,0,4,0),pt.cex=1.2,cex=0.8,border=NA,bty='n')
  
  
  for (i in 1:length(island)){
    
    brts<-island[[i]]$branching_times[-1]
    stac<-island[[i]]$stac
    
    graphics::segments(position[i],brts[1],position[i],0)
    graphics::points(position[i],brts[1],
           pch=(stac == 1)*22 + (stac == 2)*21 + (stac == 3)*21 + (stac == 4)*22 + (stac == 5)*21,
           col=(stac == 1)*1 + (stac == 2)*4 + (stac == 3)*4 + (stac == 4)*1 + (stac == 5)*4,
           bg=(stac == 1)*0 + (stac == 2)*4 + (stac == 3)*4 + (stac == 4)*1 + (stac == 5)*0,
           cex=1.5)
    graphics::text(position[i],ymax*(brts[1]/ymax+0.035),labels=paste('n=',length(brts),sep = ''),cex=0.6)
    graphics::text(position[i],ymax*(brts[1]/ymax+0.02),labels=paste('m=',island[[i]]$missing_species,sep = ''),
         cex=0.6)
    
    if(length(brts)>1) {graphics::points(rep(position[i],length(brts)-1),brts[-1],pch=3,cex=1)}
  }
  
}
