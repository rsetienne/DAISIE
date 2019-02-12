#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(daisie_fill1d)(double *vec, int *DIMP, double *parms, int *II);
extern void F77_NAME(daisie_initmod)(void (*steadyparms)(int *, double *));
extern void F77_NAME(daisie_runmod)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);
extern void F77_NAME(daisie_runmod2)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);

static const R_FortranMethodDef FortranEntries[] = {
  {"daisie_fill1d", (DL_FUNC) &F77_NAME(daisie_fill1d),  4},
  {"daisie_initmod", (DL_FUNC) &F77_NAME(daisie_initmod),  1},
  {"daisie_runmod", (DL_FUNC) &F77_NAME(daisie_runmod),  6},
  {"daisie_runmod2", (DL_FUNC) &F77_NAME(daisie_runmod2),  6},
  {NULL, NULL, 0}
};

void R_init_DAISIE(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
