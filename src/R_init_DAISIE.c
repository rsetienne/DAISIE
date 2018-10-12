#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(fill1d)(double *vec, int *DIMP, double *parms, int *II);
extern void F77_NAME(initmod)(void (*steadyparms)(int *, double *));
extern void F77_NAME(runmod)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);
extern void F77_NAME(runmod2)(int *neq, double *t, double *Conc, double *dConc, double *yout, int *ip);

static const R_FortranMethodDef FortranEntries[] = {
  {"DAISIE_fill1d", (DL_FUNC) &F77_NAME(fill1d),  4},
  {"DAISIE_initmod", (DL_FUNC) &F77_NAME(initmod),  1},
  {"DAISIE_runmod", (DL_FUNC) &F77_NAME(runmod),  6},
  {"DAISIE_runmod2", (DL_FUNC) &F77_NAME(runmod2),  6},
  {NULL, NULL, 0}
};

void R_init_DAISIE(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}