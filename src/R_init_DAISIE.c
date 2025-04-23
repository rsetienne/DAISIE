//
//  Copyright (c) 2023, Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

// [[Rcpp::plugins(cpp14)]]

#include"config.h"
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(daisie_fill1d)(double *vec, int *DIMP, double *parms, int *II);
extern void F77_NAME(daisie_initmod)(void (*steadyparms)(int *, double *));
extern void F77_NAME(daisie_initmod3)(void (*steadyparms)(int *, double *));
extern void F77_NAME(daisie_runmod)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);
extern void F77_NAME(daisie_runmod1)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);
extern void F77_NAME(daisie_runmod2)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);
extern void F77_NAME(daisie_runmod3)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);

static const R_FortranMethodDef FortranEntries[] = {
  {"daisie_fill1d", (DL_FUNC) &F77_NAME(daisie_fill1d),  4},
  {"daisie_initmod", (DL_FUNC) &F77_NAME(daisie_initmod),  1},
  {"daisie_initmod3", (DL_FUNC) &F77_NAME(daisie_initmod3),  1},
  {"daisie_runmod", (DL_FUNC) &F77_NAME(daisie_runmod),  6},
  {"daisie_runmod1", (DL_FUNC) &F77_NAME(daisie_runmod1),  6},
  {"daisie_runmod2", (DL_FUNC) &F77_NAME(daisie_runmod2),  6},
  {"daisie_runmod3", (DL_FUNC) &F77_NAME(daisie_runmod3),  6},
  {NULL, NULL, 0}
};


/* C bindings */
extern SEXP daisie_odeint_iw_num_threads(SEXP);
extern SEXP daisie_odeint_cs_max_steps(SEXP);
extern SEXP daisie_odeint_iw(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP daisie_odeint_cs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"daisie_odeint_iw_num_threads", (DL_FUNC) &daisie_odeint_iw_num_threads, 1},
  {"daisie_odeint_cs_max_steps", (DL_FUNC) &daisie_odeint_cs_max_steps, 1},
  {"daisie_odeint_iw", (DL_FUNC) &daisie_odeint_iw, 6},
  {"daisie_odeint_cs", (DL_FUNC) &daisie_odeint_cs, 9},
  {NULL, NULL, 0}
};


void R_init_DAISIE(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, FortranEntries, NULL);
  R_useDynamicSymbols(dll, TRUE);
}
