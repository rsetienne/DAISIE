#' Opens DAISIE tutorial
#' 
#' This function opens a PDF file that contains a step-by-step tutorial on how
#' to run DAISIE functions, using the Galapagos avian dataset from Valente,
#' Phillimore & Etienne 2015 as an example
#' 
#' Opens a PDF file containing the tutorial
#' 
#' @author Luis Valente and Rampal S. Etienne
#' @references Valente, L.M., A.B. Phillimore and R.S. Etienne (2015).
#' Equilibrium and non-equilibrium dynamics simultaneously operate in the
#' Galapagos islands. Ecology Letters 18: 844-852.
#' @keywords models
#' @export DAISIE_tutorial
DAISIE_tutorial = function()
{
   filename = system.file("DAISIE_tutorial.pdf",package = "DAISIE")
   os = .Platform$OS.type
   if(os == "windows")
   {
       shell.exec(filename)
   }
   if(os == "unix")
   {
       system(paste("open",filename,sep = " "))
   }
}
