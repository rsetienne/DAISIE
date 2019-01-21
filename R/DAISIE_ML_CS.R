DAISIE_ML_CS <- DAISIE_ML <- function(
     datalist,
     datatype = 'single',
     initparsopt,
     idparsopt,
     parsfix,
     idparsfix,
     idparsnoshift = 6:10,
     idparsmat = NULL,
     res = 100,
     ddmodel = 0,
     cond = 0,
     island_ontogeny = NA,
     eqmodel = 0,
     x_E = 0.95,
     x_I = 0.98,
     tol = c(1e-04, 1e-05, 1e-07),
     maxiter = 1000 * round((1.25)^length(idparsopt)),
     methode = 'lsodes',
     optimmethod = 'subplex',
     CS_version = 1,
     verbose = 0,
     tolint = c(1E-16,1E-10)
     )  
{
  if(datatype == 'single')
  {
     if(is.na(island_ontogeny))
     {
       out = DAISIE_ML1(datalist,initparsopt,idparsopt,parsfix,idparsfix,idparsnoshift,res,ddmodel,cond,eqmodel,x_E,x_I,tol,maxiter,methode,optimmethod,CS_version,verbose,tolint)
     } else
     {
       out = DAISIE_ML3(datalist,initparsopt,idparsopt,parsfix,idparsfix,res,ddmodel,cond,island_ontogeny,tol,maxiter,methode,optimmethod,CS_version,verbose,tolint)
     }
  } else
  {
     out = DAISIE_ML2(datalist,initparsopt,idparsopt,parsfix,idparsfix,idparsmat,res,ddmodel,cond,tol,maxiter,methode,optimmethod,verbose,tolint)
  }
  return(out)
}