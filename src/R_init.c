#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_RDS(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

