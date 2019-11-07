# pars_equas <- c(2.550687345, 2.683454548, Inf, 0.00933207, 1.010073119)
# island_replicates_equal <- DAISIE_sim(
#   time = 4,
#   M = 1000,
#   pars = pars_equal,
#   replicates = 2
# )
# stt_list <- island_replicates_equal[[1:2]]
# unlist(stt_list[1])
# which(stt_list)
#
# small_stts <- lapply(stt_list, nrow) == 2
# second_line_stts <- lapply(stt_list, "[", 2,)
# zeros_second_line <- sapply(second_line_stts, sum) == 0
#
# comparisson <- zeros_second_line == small_stts
# testit::assert(all(comparisson))
#
# filled_stt_lists <- stt_list[!zeros_second_line]
#
# filled_stt_lists[[1]][, 1]
# filled_stt_lists[[1]][, 1]
# times <- sapply(filled_stt_lists, "[", ,1)
#
# times_without_last <- lapply(times, head, -1)
# times_without_last_first <- sapply(times_without_last, "[", -1)
#
# sort(unlist(times_without_last_first))
