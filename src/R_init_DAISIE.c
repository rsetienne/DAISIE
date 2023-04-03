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
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
}
