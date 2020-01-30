#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "spNNGP.h"

static const R_CallMethodDef CallEntries[] = {
    {"rNNGP", (DL_FUNC) &rNNGP, 28},
    {"sNNGP", (DL_FUNC) &sNNGP, 29},
    {"sNNGPLogit", (DL_FUNC) &sNNGPLogit, 27},
    {"cNNGP", (DL_FUNC) &cNNGP, 18},
    {"cSLGP", (DL_FUNC) &cSLGP, 21},
    {"rNNGPPredict", (DL_FUNC) &rNNGPPredict, 17},
    {"sNNGPPredict", (DL_FUNC) &sNNGPPredict, 19},
    {"PGLogit", (DL_FUNC) &PGLogit, 9},
    {"rNNGPReplicated", (DL_FUNC) &rNNGPReplicated, 14},
    {NULL, NULL, 0}
};

void R_init_spNNGP(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
