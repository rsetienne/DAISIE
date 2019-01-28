DAISIE_ML = function(
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
     eqmodel = 0,
     x_E = 0.95,
     x_I = 0.98,
     tol = c(1e-04, 1e-05, 1e-07),
     maxiter = 1000 * round((1.25)^length(idparsopt)),
     methode = 'lsodes',
     optimmethod = 'subplex'
     )  
{
  if(datatype == 'single')
  {
     out = DAISIE_ML1(datalist,initparsopt,idparsopt,parsfix,idparsfix,idparsnoshift,res,ddmodel,cond,eqmodel,x_E,x_I,tol,maxiter,methode,optimmethod)
  } else
  {
     out = DAISIE_ML2(datalist,initparsopt,idparsopt,parsfix,idparsfix,idparsmat,res,ddmodel,cond,tol,maxiter,methode,optimmethod)
  }
  return(out)
}