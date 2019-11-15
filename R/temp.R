pars_equal <- c(2.550687345, 2.683454548, Inf, 0.00933207, 1.010073119)
time = 4
M = 1000
pars = pars_equal
replicates = 2
plot_sims = FALSE
island_replicates_equal <- DAISIE_sim(
  time = time,
  M = M,
  pars = pars_equal,
  replicates = replicates,
  plot_sims = FALSE
)
stt_list <- island_replicates_equal

small_stts <- lapply(stt_list, nrow) == 2
second_line_stts <- lapply(stt_list, "[", 2,)
zeros_second_line <- sapply(second_line_stts, sum) == 0

comparisson <- zeros_second_line == small_stts
testit::assert(all(comparisson))

filled_stt_lists <- stt_list[!zeros_second_line]

times_list <- sapply(filled_stt_lists, "[", , 1)
times_without_first <- sapply(times_list, "[", -1)

deltas_matrix <- lapply(filled_stt_lists, FUN = diff)
for (i in seq_along(deltas_matrix)) {
  deltas_matrix[[i]][, 1] <- times_without_first[[i]]
}

nI_list <- sapply(deltas_matrix, "[", , 2)
nA_list <- sapply(deltas_matrix, "[", , 3)
nC_list <- sapply(deltas_matrix, "[", , 4)

times <- unlist(times_without_first)
nI <- unlist(nI_list)
nA <- unlist(nA_list)
nC <- unlist(nC_list)

full_stt <- data.frame(times = times, nI = nI, nA = nA, nC = nC)
full_stt_matrix <- as.matrix.data.frame(full_stt)
ordered_diffs <- full_stt[order(full_stt$times, decreasing = TRUE), ]

complete_stt_table <- mapply(ordered_diffs[2:4], FUN = cumsum)
complete_stt_table <- cbind(ordered_diffs$times, complete_stt_table)
complete_stt_table
