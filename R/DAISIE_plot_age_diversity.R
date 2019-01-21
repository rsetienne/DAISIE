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
  plot(NULL,NULL,xlim=c(-island_age,0),ylim=c(0,2))
  legend("center", paste('Island has no species'),bty="n") 
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
plot(-colonisation_time,diversity,xlim=c(-island_age,0),ylab='Diversity',xlab='Colonisation time',col=(stacs_used==2 | stacs_used==3)*4 + (stacs_used==4)*1,pch=(stacs_used==2 | stacs_used==3)*16 + (stacs_used==4)*15,cex=1.2,main=title)
legend("topleft", paste('C=',length(stac_age_known),' ','N=',sum(diversity),sep=''),bty="n") 

}
}