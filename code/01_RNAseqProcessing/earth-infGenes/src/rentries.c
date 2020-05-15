// rentries.c: Register native routines for R.
//             The core of this routine is the function R_init_earth.

#define USING_R 1
#include "R.h"
#include "Rinternals.h" // for REALSXP etc.
#include "R_ext/Rdynload.h"
#ifndef _MSC_VER // microsoft
#ifndef bool
    typedef int bool;
#endif
#endif
#include "earth.h"

static R_NativePrimitiveArgType FreeR_t[] = {INTSXP};
static R_NativePrimitiveArgType EvalSubsetsUsingXtxR_t[] = {
    REALSXP,    // 01 double        PruneTerms[]
    REALSXP,    // 02 double        RssVec[]
    INTSXP,     // 03 const int*    pnCases
    INTSXP,     // 04 const int*    pnResp
    INTSXP,     // 05 const int*    pnMaxTerms
    REALSXP,    // 06 const double  bx[]
    REALSXP,    // 07 const double  y[]
    REALSXP     // 08 const double* pTrace
};
static R_NativePrimitiveArgType RegressR_t[] = {
    REALSXP,    // 01 double       Betas[]
    REALSXP,    // 02 double       Residuals[]
    REALSXP,    // 03 double       Rss[]
    REALSXP,    // 04 double       Diags[]
    INTSXP,     // 05 int*         pnRank
    INTSXP,     // 06 int          iPivots[]
    REALSXP,    // 07 const double x[]
    REALSXP,    // 08 const double y[]
    INTSXP,     // 09 const int*   pnCases
    INTSXP,     // 10 const int*   pnResp
    INTSXP,     // 11 int*         pnCols
    LGLSXP      // 12 const bool   UsedCols[]
};
static R_NativePrimitiveArgType ForwardPassR_t[] = {
    INTSXP,     // 01 int    FullSet[]
    REALSXP,    // 02 double bx[]
    REALSXP,    // 03 double Dirs[]
    REALSXP,    // 04 double Cuts[]
    INTSXP,     // 05 int*   piTermCond
    REALSXP,    // 06 const double x[]
    REALSXP,    // 07 const double y[]
    REALSXP,    // 08 const double yw[]
    REALSXP,    // 09 const double WeightsArg[]
    INTSXP,     // 10 const int* pnCases
    INTSXP,     // 11 const int* pnResp
    INTSXP,     // 12 const int* pnPreds
    INTSXP,     // 13 const int* pnMaxDegree
    REALSXP,    // 14 const double* pPenalty
    INTSXP,     // 15 const int* pnMaxTerms
    REALSXP,    // 16 const double* pThresh
    INTSXP,     // 17 const int* pnMinSpan
    INTSXP,     // 18 const int* pnEndSpan
    INTSXP,     // 19 const int* pnFastK
    REALSXP,    // 20 const double* pFastBeta
    REALSXP,    // 21 const double* pNewVarPenalty
    INTSXP,     // 22 const int  LinPreds[]
    CLOSXP,     // 23 const SEXP Allowed
    INTSXP,     // 24 const int* pnAllowedFuncArgs
    ENVSXP,     // 25 const SEXP Env
    REALSXP,    // 26 const double* pAdjustEndSpan
    INTSXP,     // 27 const int* pnAutoLinPreds
    INTSXP,     // 28 const int* pnUseBetaCache
    REALSXP,    // 29 const double* pTrace
    STRSXP,     // 30 const char* sPredNames[]
    REALSXP,    // 31 const double* MyNullDouble
    CLOSXP      // 32 const double* MyNullFunc
};
static R_CMethodDef cEntries[] = {
  {"FreeR",                (DL_FUNC)&FreeR,                 0, FreeR_t},
  {"EvalSubsetsUsingXtxR", (DL_FUNC)&EvalSubsetsUsingXtxR,  8, EvalSubsetsUsingXtxR_t},
  {"RegressR",             (DL_FUNC)&RegressR,             12, RegressR_t},
  {"ForwardPassR",         (DL_FUNC)&ForwardPassR,         32, ForwardPassR_t},
  {NULL,                   NULL,                            0, NULL}
};
extern void F77_SUB(bakwrd)(
    int    *NP,
    int    *NRBAR,
    double *D,
    double *RBAR,
    double *THETAB,
    int    *FIRST,
    int    *LAST,
    int    *VORDER,
    double *TOL,
    double *RSS,
    double *BOUND,
    int    *NVMAX,
    double *RESS,
    int    *IR,
    int    *NBEST,
    int    *LOPT,
    int    *IL,
    double *WK,
    int    *IWK,
    int    *IER);
extern void F77_SUB(forwrd)(
    int    *NP,
    int    *NRBAR,
    double *D,
    double *RBAR,
    double *THETAB,
    int    *FIRST,
    int    *LAST,
    int    *VORDER,
    double *TOL,
    double *RSS,
    double *BOUND,
    int    *NVMAX,
    double *RESS,
    int    *IR,
    int    *NBEST,
    int    *LOPT,
    int    *IL,
    double *WK,
    int    *IWK,
    int    *IER);
extern void F77_SUB(seqrep)(
    int    *NP,
    int    *NRBAR,
    double *D,
    double *RBAR,
    double *THETAB,
    int    *FIRST,
    int    *LAST,
    int    *VORDER,
    double *TOL,
    double *RSS,
    double *BOUND,
    int    *NVMAX,
    double *RESS,
    int    *IR,
    int    *NBEST,
    int    *LOPT,
    int    *IL,
    double *WK,
    int    *IWK,
    int    *IER);
void F77_SUB(xhaust)(
    int    *NP,
    int    *NRBAR,
    double *D,
    double *RBAR,
    double *THETAB,
    int    *FIRST,
    int    *LAST,
    int    *VORDER,
    double *TOL,
    double *RSS,
    double *BOUND,
    int    *NVMAX,
    double *RESS,
    int    *IR,
    int    *NBEST,
    int    *LOPT,
    int    *IL,
    double *WK,
    int    *DIMWK,
    int    *IWK,
    int    *DIMIWK,
    int    *IER);
extern void F77_SUB(initr)(
    int    *NP,
    int    *NVMAX,
    int    *NBEST,
    double *BOUND,
    double *RESS,
    int    *IR,
    int    *LOPT,
    int    *IL,
    int    *VORDER,
    double *RSS,
    int    *IER);
extern void F77_SUB(sing)(
    int    *NP,
    int    *NRBAR,
    double *D,
    double *RBAR,
    double *THETAB,
    double *SSERR,
    double *TOL,
    bool   *LINDEP,
    double *WORK,
    int    *IER);
extern void F77_SUB(ssleaps)(
    int    *NP,
    double *D,
    double *THETAB,
    double *SSERR,
    double *RSS,
    int    *IER);
extern void F77_SUB(tolset)(
    int    *NP,
    int    *NRBAR,
    double *D,
    double *RBAR,
    double *TOL,
    double *WORK,
    int    *IER);
void F77_SUB(makeqr)(
    int    *NP,
    int    *NN,
    double *WEIGHTS,
    double *TXMAT,
    double *YVEC,
    double *D,
    double *RBAR,
    double *THETAB,
    double *SSERR,
    int    *IER);
static R_NativePrimitiveArgType bakwrd_t[] = {
    INTSXP,     // 01 INTEGER NP
    INTSXP,     // 02 INTEGER NRBAR
    REALSXP,    // 03 DOUBLE  D
    REALSXP,    // 04 DOUBLE  RBAR
    REALSXP,    // 05 DOUBLE  THETAB
    INTSXP,     // 06 INTEGER FIRST
    INTSXP,     // 07 INTEGER LAST
    INTSXP,     // 08 INTEGER VORDER
    REALSXP,    // 09 DOUBLE  TOL
    REALSXP,    // 10 DOUBLE  RSS
    REALSXP,    // 11 DOUBLE  BOUND
    INTSXP,     // 12 INTEGER NVMAX
    REALSXP,    // 13 DOUBLE  RESS
    INTSXP,     // 14 INTEGER IR
    INTSXP,     // 15 INTEGER NBEST
    INTSXP,     // 16 INTEGER LOPT
    INTSXP,     // 17 INTEGER IL
    REALSXP,    // 18 DOUBLE  WK
    INTSXP,     // 19 INTEGER IWK
    INTSXP,     // 20 INTEGER IER
};
static R_NativePrimitiveArgType forwrd_t[] = {
    INTSXP,     // 01 INTEGER NP
    INTSXP,     // 02 INTEGER NRBAR
    REALSXP,    // 03 DOUBLE  D
    REALSXP,    // 04 DOUBLE  RBAR
    REALSXP,    // 05 DOUBLE  THETAB
    INTSXP,     // 06 INTEGER FIRST
    INTSXP,     // 07 INTEGER LAST
    INTSXP,     // 08 INTEGER VORDER
    REALSXP,    // 09 DOUBLE  TOL
    REALSXP,    // 10 DOUBLE  RSS
    REALSXP,    // 11 DOUBLE  BOUND
    INTSXP,     // 12 INTEGER NVMAX
    REALSXP,    // 13 DOUBLE  RESS
    INTSXP,     // 14 INTEGER IR
    INTSXP,     // 15 INTEGER NBEST
    INTSXP,     // 16 INTEGER LOPT
    INTSXP,     // 17 INTEGER IL
    REALSXP,    // 18 DOUBLE  WK
    INTSXP,     // 19 INTEGER IWK
    INTSXP,     // 20 INTEGER IER
};
static R_NativePrimitiveArgType seqrep_t[] = {
    INTSXP,     // 01 INTEGER NP
    INTSXP,     // 02 INTEGER NRBAR
    REALSXP,    // 03 DOUBLE  D
    REALSXP,    // 04 DOUBLE  RBAR
    REALSXP,    // 05 DOUBLE  THETAB
    INTSXP,     // 06 INTEGER FIRST
    INTSXP,     // 07 INTEGER LAST
    INTSXP,     // 08 INTEGER VORDER
    REALSXP,    // 09 DOUBLE  TOL
    REALSXP,    // 10 DOUBLE  RSS
    REALSXP,    // 11 DOUBLE  BOUND
    INTSXP,     // 12 INTEGER NVMAX
    REALSXP,    // 13 DOUBLE  RESS
    INTSXP,     // 14 INTEGER IR
    INTSXP,     // 15 INTEGER NBEST
    INTSXP,     // 16 INTEGER LOPT
    INTSXP,     // 17 INTEGER IL
    REALSXP,    // 18 DOUBLE  WK
    INTSXP,     // 19 INTEGER IWK
    INTSXP,     // 20 INTEGER IER
};
static R_NativePrimitiveArgType xhaust_t[] = {
    INTSXP,     // 01 INTEGER NP
    INTSXP,     // 02 INTEGER NRBAR
    REALSXP,    // 03 DOUBLE  D
    REALSXP,    // 04 DOUBLE  RBAR
    REALSXP,    // 05 DOUBLE  THETAB
    INTSXP,     // 06 INTEGER FIRST
    INTSXP,     // 07 INTEGER LAST
    INTSXP,     // 08 INTEGER VORDER
    REALSXP,    // 09 DOUBLE  TOL
    REALSXP,    // 10 DOUBLE  RSS
    REALSXP,    // 11 DOUBLE  BOUND
    INTSXP,     // 12 INTEGER NVMAX
    REALSXP,    // 13 DOUBLE  RESS
    INTSXP,     // 14 INTEGER IR
    INTSXP,     // 15 INTEGER NBEST
    INTSXP,     // 16 INTEGER LOPT
    INTSXP,     // 17 INTEGER IL
    REALSXP,    // 18 DOUBLE  WK
    INTSXP,     // 19 INTEGER DIMWK
    INTSXP,     // 20 INTEGER IWK
    INTSXP,     // 21 INTEGER DIMIWK
    INTSXP,     // 22 INTEGER IER
};
static R_NativePrimitiveArgType initr_t[] = {
    INTSXP,     // 01 INTEGER NP
    INTSXP,     // 02 INTEGER NVMAX
    INTSXP,     // 03 INTEGER NBEST
    REALSXP,    // 09 DOUBLE  BOUND
    REALSXP,    // 10 DOUBLE  RESS
    INTSXP,     // 04 INTEGER IR
    INTSXP,     // 06 INTEGER LOPT
    INTSXP,     // 05 INTEGER IL
    INTSXP,     // 07 INTEGER VORDER
    REALSXP,    // 11 DOUBLE  RSS
    INTSXP,     // 08 INTEGER IER
};
static R_NativePrimitiveArgType sing_t[] = {
    INTSXP,     // 01 DOUBLE  NP
    INTSXP,     // 02 DOUBLE  NRBAR
    REALSXP,    // 04 DOUBLE  D
    REALSXP,    // 05 DOUBLE  RBAR
    REALSXP,    // 06 DOUBLE  THETAB
    REALSXP,    // 07 DOUBLE  SSERR
    REALSXP,    // 08 DOUBLE  TOL
    LGLSXP,     // 10 LOGICAL LINDEP
    REALSXP,    // 09 DOUBLE  WORK
    INTSXP,     // 03 INTEGER IER
};
static R_NativePrimitiveArgType ssleaps_t[] = {
    INTSXP,     // 01 INTEGER NP
    REALSXP,    // 03 DOUBLE  D
    REALSXP,    // 04 DOUBLE  THETAB
    REALSXP,    // 05 DOUBLE  SSERR
    REALSXP,    // 06 DOUBLE  RSS
    INTSXP,     // 02 INTEGER IER
};
static R_NativePrimitiveArgType tolset_t[] = {
    INTSXP,     // 01 INTEGER NP
    INTSXP,     // 02 INTEGER NRBAR
    REALSXP,    // 04 DOUBLE  D
    REALSXP,    // 05 DOUBLE  RBAR
    REALSXP,    // 06 DOUBLE  TOL
    REALSXP,    // 07 DOUBLE  WORK
    INTSXP,     // 03 INTEGER IER
};
static R_NativePrimitiveArgType makeqr_t[] = {
    INTSXP,     // 01 INTEGER NP
    INTSXP,     // 02 INTEGER NN
    REALSXP,    // 04 DOUBLE WEIGHTS
    REALSXP,    // 05 DOUBLE TXMAT
    REALSXP,    // 06 DOUBLE YVEC
    REALSXP,    // 07 DOUBLE D
    REALSXP,    // 08 DOUBLE RBAR
    REALSXP,    // 09 DOUBLE THETAB
    REALSXP,    // 10 DOUBLE SSERR
    INTSXP,     // 03 INTEGER IER
};
static R_FortranMethodDef fortranEntries[] = {
    { "bakwrd",  (DL_FUNC)&F77_SUB(bakwrd),  20, bakwrd_t},
    { "forwrd",  (DL_FUNC)&F77_SUB(forwrd),  20, forwrd_t},
    { "seqrep",  (DL_FUNC)&F77_SUB(seqrep),  20, seqrep_t},
    { "xhaust",  (DL_FUNC)&F77_SUB(xhaust),  22, xhaust_t},
    { "initr",   (DL_FUNC)&F77_SUB(initr),   11, initr_t},
    { "sing",    (DL_FUNC)&F77_SUB(sing),    10, sing_t},
    { "ssleaps", (DL_FUNC)&F77_SUB(ssleaps),  6, ssleaps_t},
    { "tolset",  (DL_FUNC)&F77_SUB(tolset),   7, tolset_t},
    { "makeqr",  (DL_FUNC)&F77_SUB(makeqr),  10, makeqr_t},
    {NULL,       NULL,                        0, NULL}
};
void R_init_earth(DllInfo *dll) // called by R after R loads the earth package
{
    R_registerRoutines(dll, cEntries, NULL, fortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
