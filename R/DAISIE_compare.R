#' Comparing simulated results from DAISIE and TRASIE
#'
#' @param DAISIE_output List with simulation output from general DAISIE
#'   \code{DAISIE_sim}
#' @param TRASIE_output1 List with simulation output using TRASIE pars
#'   \code{DAISIE_sim}
#' @param TRASIE_output2 List with simulation output using TRASIE pars
#'   \code{DAISIE_sim}
#' @param TRASIE_output3 List with simulation output using TRASIE pars
#'   \code{DAISIE_sim}
#' @param Tpars1 A named list containing diversification rates considering two trait states:
#' \itemize{
#'   \item{[1]:A numeric with the per capita transition rate with state1}
#'   \item{[2]:A numeric with the per capita immigration rate with state2}
#'   \item{[3]:A numeric with the per capita extinction rate with state2}
#'   \item{[4]:A numeric with the per capita anagenesis rate with state2}
#'   \item{[5]:A numeric with the per capita cladogenesis rate with state2}
#'   \item{[6]:A numeric with the per capita transition rate with state2}
#'   \item{[7]:A numeric with the number of species with trait state 2 on mainland}
#' }
#' @param Tpars2 A named list containing diversification rates considering two trait states:
#' @param Tpars3 A named list containing diversification rates considering two trait states:
#' @seealso \code{\link{DAISIE_sim}},
#'   \code{\link{DAISIE_convert_to_classic_plot}}
#' @family plotting
#' @export
#'
compare_diff <- function(
  DAISIE_output,
  TRASIE_output1,
  TRASIE_output2,
  TRASIE_output3,
  Tpars1,
  Tpars2,
  Tpars3){

  DAISIE_list = DAISIE_convert_to_classic_plot(DAISIE_output,Tpars = NULL)
  TRASIE_list1 = DAISIE_convert_to_classic_plot(TRASIE_output1,Tpars = Tpars1)
  TRASIE_list2 = DAISIE_convert_to_classic_plot(TRASIE_output2,Tpars = Tpars2)
  TRASIE_list3 = DAISIE_convert_to_classic_plot(TRASIE_output3,Tpars = Tpars3)
  # TRASIE_list4 = DAISIE_convert_to_classic_plot(TRASIE_output4,Tpars = Tpars4)

  stt1 = DAISIE_list[["all_species"]]
  stt2 = TRASIE_list1[["all_species"]]
  stt3 = TRASIE_list2[["all_species"]]
  stt4 = TRASIE_list3[["all_species"]]
  # stt5 = TRASIE_list4[["all_species"]]

  ###### Compare total species number
  # DAISIE_stt = stt1$stt_average[, "Total"]
  # TRASIE_notrans_stt = stt2$stt_average[, "Total"]
  # TRASIE_trans1_stt = stt3$stt_average[, "Total"]
  # TRASIE_trans2_stt = stt4$stt_average[, "Total"]
  # TRASIE_trans3_stt = stt5$stt_average[, "Total"]

  Data_total<- data.frame(stt1 = stt1$stt_average[, "Total"],
                          stt2 = stt2$stt_average[, "Total"],
                          stt3 = stt3$stt_average[, "Total"],
                          stt4 = stt4$stt_average[, "Total"])
  Data_nonendemic <- data.frame(stt1 = stt1$stt_average[, "nI"],
                                stt2 = stt2$stt_average[, "nI"],
                                stt3 = stt3$stt_average[, "nI"],
                                stt4 = stt4$stt_average[, "nI"])
  Data_Endemic <- data.frame(stt1 = stt1$stt_average[, "Endemic"],
                             stt2 = stt2$stt_average[, "Endemic"],
                             stt3 = stt3$stt_average[, "Endemic"],
                             stt4 = stt4$stt_average[, "Endemic"])
  Data <- list(Data_total, Data_Endemic, Data_nonendemic)

  for(i in 1:length(Data)){
    Data[[i]]$diff1_2 <- Data[[i]]$stt2 - Data[[i]]$stt1
    Data[[i]]$diff1_3 <- Data[[i]]$stt3 - Data[[i]]$stt1
    Data[[i]]$diff1_4 <- Data[[i]]$stt4 - Data[[i]]$stt1
    # Data[[i]]$diff1_5 <- Data[[i]]$stt5 - Data[[i]]$stt1

    Data_differ_quantile <- apply(Data[[i]][,5:7],2,function(x) round(quantile(x,c(0.25,0.5,0.75)),digits=2))
    rownames(Data_differ_quantile)[1:3] <-c("differ_down","differ_median","differ_up")
    Data_t <- data.frame(t(Data_differ_quantile))

    pairs_name <- c("DAISIE vs. TRASIE1","DAISIE vs. TRASIE2","DAISIE vs. TRASIE3")
    pairs_median <- Data_t$differ_median
    pairs_CI <- paste("(", Data_t$differ_down, " ~ ", Data_t$differ_up, ")", sep = "")
    Data_str <- data.frame(pairs_na=pairs_name,pairs_median=pairs_median,pairs_CI=pairs_CI)
    Data_str <- as.matrix(Data_str)
    Data_str <- rbind(c(NA,"Median","Interquartile ranges"),Data_str)
    # jpeg(file = "results_Value_1.jpg",width =2000,height = 1800,units = "px",res =300)
    forestplot::forestplot(Data_str,
               c(NA,Data_t$differ_median),
               c(NA,Data_t$differ_down),
               c(NA,Data_t$differ_up),
               zero = 0,
               xlog=FALSE,
               fn.ci_norm = forestplot::fpDrawCircleCI,
               boxsize = 0.3,
               col = forestplot::fpColors(line = "#CC79A7",
                            box="#D55E00"),
               lty.ci = 7,
               lwd.ci = 3,
               ci.vertices.height = 0.15,
               txt_gp = forestplot::fpTxtGp(ticks = grid::gpar(cex = 0.5), xlab = grid::gpar(cex = 0.7), cex = 0.7),
               lineheight = "auto",
               xlab="Differences in assessment indicators between relevant pairs"
    )
    # dev.off()
  }
}
