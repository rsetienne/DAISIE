#' @description DAISIE is an island biogeography model that assumes island biota
#' assembly is governed by immigration, extinction and speciation through cladogenesis
#' and anagenesis. This likelihood-based statistical package can simulate islands and
#' estimate parameters of the DAISIE model based on phylogenetic/phylogeographic data.
#' Cladogenesis and immigration rates can be dependent on diversity.
#' @references
#' \itemize{
#' \item Valente, L. M., Etienne, R. S., & Phillimore, A. B. (2014). The effects of island ontogeny on species diversity and phylogeny. Proceedings of the Royal Society of London. Series B, Biological Sciences 281, 20133227. \doi{10.1098/rspb.2013.3227}. \cr
#' \item Valente, L., A.B. Phillimore & R.S. Etienne (2015). Equilibrium and non-equilibrium dynamics simultaneously operate in the Galapagos islands. Ecology Letters 18: 844-852, \doi{10.1111/ele.12461}. \cr
#' \item Valente, L., R.S. Etienne & L. Davalos (2017). Recent extinctions disturb path to equilibrium diversity in Caribbean bats. Nature Ecology and Evolution 1: 0026. \doi{10.1038/s41559-016-0026}.\cr
#' \item Valente, L., Illera, J. C., Havenstein, K., Pallien, T., Etienne, R. S., & Tiedemann, R. (2017). Equilibrium Bird Species Diversity in Atlantic Islands. Current Biology 27: 1660-1666. \doi{10.1016/j.cub.2017.04.053}.\cr
#' \item Valente, L., Phillimore, A. B., & Etienne, R. (2018). Using molecular phylogenies in island biogeography: It's about time. Ecography 41: 1684-1686. \doi{10.1111/ecog.03503}.\cr
#' \item Valente, L., Etienne, R. S., & Garcia-R, J. C. (2019). Deep macroevolutionary impact of humans on New Zealand's unique avifauna. Current Biology 29: 2563-2569. \doi{10.1016/j.cub.2019.06.058}.\cr
#' \item Valente, L., Phillimore, A. B., Melo, M., Warren, B. H., Clegg, S. M., Havenstein, K., & Etienne, R. S. (2020). A simple dynamic model explains the diversity of island birds worldwide. Nature 579: 92-96. \doi{10.1038/s41586-020-2022-5}.\cr
#' \item Hauffe, T., Delicado, D., Etienne, R.S., & Valente, L. (2020). Lake expansion elevates equilibrium diversity via increasing colonization. Journal of Biogeography 47: 1849–1860. \doi{10.1111/jbi.13914}.\cr
#' \item Valente, L., Kristensen, N., Phillimore, A. B., & Etienne, R. S. (2021). Report of programming bugs in the DAISIE R package: consequences and correction. EcoEvoRxiv. \doi{10.32942/osf.io/w5ntf}.\cr
#' \item Santos Neves, P., Lambert, J. W., Valente, L., & Etienne, R. S. (2021).The robustness of a simple dynamic model of island biodiversity to geological and eustatic change. bioRxiv 2021.07.26.453064. \doi{10.1101/2021.07.26.453064}.\cr
#' }
#' @keywords internal
#' @import Rcpp
#' @useDynLib DAISIE, .registration = TRUE
"_PACKAGE"
