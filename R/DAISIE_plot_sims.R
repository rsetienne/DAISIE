#' @title Plot island species-through-time (STT) plots
#' @description Produces STT plots. If only one type of species is present in the simulated
#' islands, STT is plotted for all species. If two types are present, three
#' plots are produced: STT for all, STT for type 1 and STT for type 2.
#' 
#' R plots with number of total, endemic and non-endemic STTs for different
#' types of species for the entire time span the islands were simulated.
#' 2.5-97.5th percentiles are plotted in light grey, 25-75th percentiles
#' plotted in dark grey.
#' 
#' @param island_replicates Island replicates in DAISIE format (Produced in
#' DAISIE_sim with format=TRUE option, or in DAISIE_format_sim)
#' @return R plot.
#' @author Luis Valente
#' @seealso \code{\link{DAISIE_sim}}
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#' @keywords models
#' @examples
#'  
#' 
#' ### Plot islands with single process (only one type of species)
#' utils::data(islands_1type_1000reps)
#' DAISIE_plot_sims(island_replicates = islands_1type_1000reps)
#' 
#' 
#' ### Plot island with type 1 and type 2
#' utils::data(islands_2types_1000reps)
#' DAISIE_plot_sims(island_replicates = islands_2types_1000reps)
#' 
#' 
#' @export DAISIE_plot_sims
DAISIE_plot_sims<-function(island_replicates)
{

replicates<-length(island_replicates)
time<-max(island_replicates[[1]][[1]]$stt_all[,1])

###  STT ALL species
s_freq<-length(island_replicates[[1]][[1]]$stt_all[,1])
complete_arr <-array(dim=c(s_freq,6,replicates))

for(x in 1:replicates)
{
sum_endemics<-island_replicates[[x]][[1]]$stt_all[,"nA"]+island_replicates[[x]][[1]]$stt_all[,"nC"]
total<-island_replicates[[x]][[1]]$stt_all[,"nA"]+island_replicates[[x]][[1]]$stt_all[,"nC"]+island_replicates[[x]][[1]]$stt_all[,"nI"]
complete_arr[,,x]<-cbind(island_replicates[[x]][[1]]$stt_all[,c('Time',"nI","nA","nC")],sum_endemics,total)
}


stt_average_all<-apply(complete_arr,c(1,2),stats::median)
stt_q0.025_all<-apply(complete_arr,c(1,2),stats::quantile,0.025)
stt_q0.25_all<-apply(complete_arr,c(1,2),stats::quantile,0.25)
stt_q0.75_all<-apply(complete_arr,c(1,2),stats::quantile,0.75)
stt_q0.975_all<-apply(complete_arr,c(1,2),stats::quantile,0.975)

colnames(stt_average_all)<-c("Time","nI","nA","nC","Endemic","Total")
colnames(stt_q0.025_all)<-c("Time","nI","nA","nC","Endemic","Total")
colnames(stt_q0.25_all)<-c("Time","nI","nA","nC","Endemic","Total")
colnames(stt_q0.75_all)<-c("Time","nI","nA","nC","Endemic","Total")
colnames(stt_q0.975_all)<-c("Time","nI","nA","nC","Endemic","Total")


if(is.null(island_replicates[[1]][[1]]$stt_type1)==FALSE)
{



###  STT TYPE1
s_freq<-length(island_replicates[[1]][[1]]$stt_type1[,1])
complete_arr <-array(dim=c(s_freq,7,replicates))

for(x in 1:replicates)
{
sum_endemics<-island_replicates[[x]][[1]]$stt_type1[,"nA"]+island_replicates[[x]][[1]]$stt_type1[,"nC"]
total<-island_replicates[[x]][[1]]$stt_type1[,"nA"]+island_replicates[[x]][[1]]$stt_type1[,"nC"]+island_replicates[[x]][[1]]$stt_type1[,"nI"]
complete_arr[,,x]<-cbind(island_replicates[[x]][[1]]$stt_type1,sum_endemics,total)
}


stt_average_type1<-apply(complete_arr,c(1,2),stats::median)
stt_q0.025_type1<-apply(complete_arr,c(1,2),stats::quantile,0.025)
stt_q0.25_type1<-apply(complete_arr,c(1,2),stats::quantile,0.25)
stt_q0.75_type1<-apply(complete_arr,c(1,2),stats::quantile,0.75)
stt_q0.975_type1<-apply(complete_arr,c(1,2),stats::quantile,0.975)

colnames(stt_average_type1)<-c("Time","nI","nA","nC","present","Endemic","Total")
colnames(stt_q0.025_type1)<-c("Time","nI","nA","nC","present","Endemic","Total")
colnames(stt_q0.25_type1)<-c("Time","nI","nA","nC","present","Endemic","Total")
colnames(stt_q0.75_type1)<-c("Time","nI","nA","nC","present","Endemic","Total")
colnames(stt_q0.975_type1)<-c("Time","nI","nA","nC","present","Endemic","Total")



###  STT TYPE2
s_freq<-length(island_replicates[[1]][[1]]$stt_type2[,1])
complete_arr <-array(dim=c(s_freq,7,replicates))

for(x in 1:replicates)
{
sum_endemics<-island_replicates[[x]][[1]]$stt_type2[,"nA"]+island_replicates[[x]][[1]]$stt_type2[,"nC"]
total<-island_replicates[[x]][[1]]$stt_type2[,"nA"]+island_replicates[[x]][[1]]$stt_type2[,"nC"]+island_replicates[[x]][[1]]$stt_type2[,"nI"]
complete_arr[,,x]<-cbind(island_replicates[[x]][[1]]$stt_type2,sum_endemics,total)
}

stt_average_type2<-apply(complete_arr,c(1,2),stats::median)
stt_q0.025_type2<-apply(complete_arr,c(1,2),stats::quantile,0.025)
stt_q0.25_type2<-apply(complete_arr,c(1,2),stats::quantile,0.25)
stt_q0.75_type2<-apply(complete_arr,c(1,2),stats::quantile,0.75)
stt_q0.975_type2<-apply(complete_arr,c(1,2),stats::quantile,0.975)

colnames(stt_average_type2)<-c("Time","nI","nA","nC","present","Endemic","Total")
colnames(stt_q0.025_type2)<-c("Time","nI","nA","nC","present","Endemic","Total")
colnames(stt_q0.25_type2)<-c("Time","nI","nA","nC","present","Endemic","Total")
colnames(stt_q0.75_type2)<-c("Time","nI","nA","nC","present","Endemic","Total")
colnames(stt_q0.975_type2)<-c("Time","nI","nA","nC","present","Endemic","Total")

grDevices::dev.new(width=12, height=4)
graphics::par(mfrow=c(1,3)) 

suppressWarnings(graphics::plot(NULL,NULL,xlim=rev(c(0, time)),ylim=c(1,max(stt_q0.975_all)),ylab="No of species + 1",bty="l", xaxs="i",xlab="Time before present",main="Species-through-time - All species",log='y',cex.lab=1.5,cex.main=1.2,cex.axis=1.2))
graphics::polygon(c(stt_average_all[,"Time"],rev(stt_average_all[,"Time"])),c(stt_q0.025_all[,"Total"]+1,rev(stt_q0.975_all[,"Total"]+1)),col="light grey",border=NA)
graphics::polygon(c(stt_average_all[,"Time"],rev(stt_average_all[,"Time"])),c(stt_q0.25_all[,"Total"]+1,rev(stt_q0.75_all[,"Total"]+1)),col="dark grey",border=NA)
graphics::lines(stt_average_all[,"Time"],stt_average_all[,"Total"]+1,lwd=2)
graphics::lines(stt_average_all[,"Time"],stt_average_all[,"nI"]+1,lwd=2,col='cyan3')
graphics::lines(stt_average_all[,"Time"],stt_average_all[,"Endemic"]+1,lwd=2,col='dodgerblue1')

graphics::legend(time,max(stt_q0.975_all),c("Total","Non-endemic","Endemic"),lty=1,lwd=2,col=c("black","cyan3","dodgerblue1"),cex=1.2,border=NA,bty='n')

suppressWarnings(graphics::plot(NULL,NULL,xlim=rev(c(0, time)),ylim=c(1,max(stt_q0.975_type1)),ylab="No of species + 1",bty="l", xaxs="i",xlab="Time before present",main="STT Type 1 species",log='y',cex.lab=1.5,cex.main=1.2,cex.axis=1.2))
graphics::polygon(c(stt_average_type1[,"Time"],rev(stt_average_type1[,"Time"])),c(stt_q0.025_type1[,"Total"]+1,rev(stt_q0.975_type1[,"Total"]+1)),col="light grey",border=NA)
graphics::polygon(c(stt_average_type1[,"Time"],rev(stt_average_type1[,"Time"])),c(stt_q0.25_type1[,"Total"]+1,rev(stt_q0.75_type1[,"Total"]+1)),col="dark grey",border=NA)
graphics::lines(stt_average_type1[,"Time"],stt_average_type1[,"Total"]+1,lwd=2)
graphics::lines(stt_average_type1[,"Time"],stt_average_type1[,"nI"]+1,lwd=2,col='cyan3')
graphics::lines(stt_average_type1[,"Time"],stt_average_type1[,"Endemic"]+1,lwd=2,col='dodgerblue1')

suppressWarnings(graphics::plot(NULL,NULL,xlim=rev(c(0, time)),ylim=c(1,max(stt_q0.975_type2)),ylab="No of species + 1",bty="l", xaxs="i",xlab="Time before present",main="STT Type 2 species",log='y',cex.lab=1.5,cex.main=1.2,cex.axis=1.2))
graphics::polygon(c(stt_average_type2[,"Time"],rev(stt_average_type2[,"Time"])),c(stt_q0.025_type2[,"Total"]+1,rev(stt_q0.975_type2[,"Total"]+1)),col="light grey",border=NA)
graphics::polygon(c(stt_average_type2[,"Time"],rev(stt_average_type2[,"Time"])),c(stt_q0.25_type2[,"Total"]+1,rev(stt_q0.75_type2[,"Total"]+1)),col="dark grey",border=NA)
graphics::lines(stt_average_type2[,"Time"],stt_average_type2[,"Total"]+1,lwd=2)
graphics::lines(stt_average_type2[,"Time"],stt_average_type2[,"nI"]+1,lwd=2,col='cyan3')
graphics::lines(stt_average_type2[,"Time"],stt_average_type2[,"Endemic"]+1,lwd=2,col='dodgerblue1')

}else{
grDevices::dev.new(width=6, height=6)
graphics::par(mfrow=c(1,1)) 
suppressWarnings(graphics::plot(NULL,NULL,xlim=rev(c(0, time)),ylim=c(1,max(stt_q0.975_all)),ylab="No of species + 1",bty="l", xaxs="i",xlab="Time before present",main="Species-through-time - All species",log='y',cex.lab=1.2,cex.main=1.2,cex.axis=1.2))
graphics::polygon(c(stt_average_all[,"Time"],rev(stt_average_all[,"Time"])),c(stt_q0.025_all[,"Total"]+1,rev(stt_q0.975_all[,"Total"]+1)),col="light grey",border=NA)
graphics::polygon(c(stt_average_all[,"Time"],rev(stt_average_all[,"Time"])),c(stt_q0.25_all[,"Total"]+1,rev(stt_q0.75_all[,"Total"]+1)),col="dark grey",border=NA)
graphics::lines(stt_average_all[,"Time"],stt_average_all[,"Total"]+1,lwd=2)
graphics::lines(stt_average_all[,"Time"],stt_average_all[,"nI"]+1,lwd=2,col='cyan3')
graphics::lines(stt_average_all[,"Time"],stt_average_all[,"Endemic"]+1,lwd=2,col='dodgerblue1')

graphics::legend(time,max(stt_q0.975_all),c("Total","Non-endemic","Endemic"),lty=1,lwd=2,col=c("black","cyan3","dodgerblue1"),cex=1.2,border=NA,bty='n')

}

}
