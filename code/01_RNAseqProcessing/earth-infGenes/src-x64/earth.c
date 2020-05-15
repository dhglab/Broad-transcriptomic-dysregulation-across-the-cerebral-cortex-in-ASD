// earth.c
//
// This code is derived from code in the Rational Fortran file dmarss.r which is
// part of the R and S mda package by Hastie and Tibshirani.
// Comments containing "TODO" mark known issues
//
// See the R earth documentation for descriptions of the principal data structures.
// See also www.milbo.users.sonic.net.  This code uses a subset of C99.
//
// Stephen Milborrow Feb 2007 Petaluma
//
// References:
//
// HastieTibs: Trevor Hastie and Robert Tibshirani
//      S library mda version 0.3.2 dmarss.r Ratfor code
//      Modifications for R by Kurt Hornik, Friedrich Leisch, Brian Ripley
//
// FriedmanMars: Multivariate Adaptive Regression Splines (with discussion)
//      Annals of Statistics 19/1, 1--141, 1991
//
// FriedmanFastMars: Friedman "Fast MARS"
//      Dep. of Stats. Stanford, Tech Report 110, May 1993
//
// Miller: Alan Miller (2nd ed. 2002) Subset Selection in Regression
//
//-----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// A copy of the GNU General Public License is available at
// http://www.r-project.org/Licenses
//
//-----------------------------------------------------------------------------

#if !STANDALONE
#define USING_R 1
#endif // STANDALONE

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <math.h>
#if _MSC_VER            // microsoft
    #include <crtdbg.h> // microsoft malloc debugging library
    #define _C_ "C"
    // disable warning: 'vsprintf': This function or variable may be unsafe
    #pragma warning(disable: 4996)
#else
    #define _C_
    #ifndef bool
        typedef int bool;
        #define false 0
        #define true  1
    #endif
#endif

#if USING_R             // R with gcc
    #include "R.h"
    #include "Rinternals.h" // needed for Allowed function handling
    #include "allowed.h"
    #include "R_ext/Rdynload.h"
    #define printf Rprintf
    #define FINITE(x) R_FINITE(x)
#else
    #define warning printf
    void error(const char* args, ...);
#if _MSC_VER // microsoft
    #define ISNAN(x)  _isnan(x)
    #define FINITE(x) _finite(x)
#else
    #define ISNAN(x)  isnan(x)
    #define FINITE(x) finite(x)
#endif
#endif

#include "earth.h"
#ifdef MATLAB
    #include "mex.h" // for printf
#endif

// linpack functions (dqrdc2 is a modified form of dqrdc used by R)
extern _C_ int dqrdc2_(double* x, int* ldx, int* n, int* p,
                       double* tol, int* rank,
                       double* qraux, int* pivot, double* work);
extern _C_ int dtrsl_(double* t, int* ldt, int* n, double* b, int* job, int* info);
extern _C_ int dqrsl_(double* x, int* ldx, int* n, int* k,
                      double* qraux, double* y,
                      double* qy, double* qty, double* b,
                      double* rsd, double* xb, int* job, int* info);
extern _C_ int daxpy_(const int* n, const double* alpha,
                      const double* dx, const int* incx,
                      double* dy, const int* incy);
extern _C_ double ddot_(const int* n,
                        const double* dx, const int* incx,
                        const double* dy, const int* incy);

#define ASSERT(x) \
    if(!(x)) error("internal assertion failed in file %s line %d: %s\n", \
                   __FILE__, __LINE__, #x)

#define sq(x)       ((x) * (x))
#ifndef max
#define max(a,b)    (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)    (((a) < (b)) ? (a) : (b))
#endif

#define INLINE      inline

#define USE_BLAS    1 // 1 is faster (tested on Windows XP Pentium with R BLAS)
                      // also, need USE_BLAS to use bxOrthCenteredT

#define FAST_MARS   1 // 1 to use techniques in FriedmanFastMars (see refs)

#define WEIGHTS     1 // 1 if case weights are supported

#if STANDALONE
#define IOFFSET     0
#else
#define IOFFSET     1 // printfs only: 1 to convert 0-based indices to 1-based in printfs
                      // use 1 for R style indices in messages to the user
                      // use 0 for C style indices in messages to the user
#endif

static const char*  VERSION     = "version 4.6.0"; // change if you modify this file!
static const double MIN_GRSQ    = -10.0;
static const double QR_TOL      = 1e-8;  // same as R lm
static const double MIN_BX_SOS  = .01;
static const double ALMOST_ZERO = 1e-10;
static const int    ONE         = 1;     // parameter for BLAS routines
#ifdef __INTEL_COMPILER
static const double POS_INF     = __builtin_inff();
#elif _MSC_VER                           // microsoft compiler
static const double POS_INF     = _HUGE;
#else                                    // assume gcc
static const double POS_INF     = (1.0 / 0.0);
#endif
static const int    MAX_DEGREE  = 100;

// Poor man's array indexing -- not pretty, but pretty useful.
//
// Note that we use column major ordering. C programs usually use row major
// ordering but we don't here because the functions in this file are called
// by R and call Fortran routines which use column major ordering.
//
// Note that nCases is size_t (not int), allowing array indices bigger than 2GB.
// We don't expect nCases itself to be that big, but it can be used in
// expressions that evaluate to more than 2GB.

#define Dirs_(iTerm,iPred)      Dirs[(iTerm) + (iPred)*(nMaxTerms)]
#define Cuts_(iTerm,iPred)      Cuts[(iTerm) + (iPred)*(nMaxTerms)]

#define bx_(i,iTerm)                bx             [(i) + (iTerm)*(nCases)]
#define bxUsed_(i,iTerm)            bxUsed         [(i) + (iTerm)*(nCases)]
#define bxOrth_(i,iTerm)            bxOrth         [(i) + (iTerm)*(nCases)]
#define bxOrthCenteredT_(iTerm,i)   bxOrthCenteredT[(iTerm) + (i)*(nMaxTerms)]
#define x_(i,iPred)                 x              [(i) + (iPred)*(nCases)]
#define xOrder_(i,iPred)            xOrder         [(i) + (iPred)*(nCases)]
#define y_(i,iResp)                 y              [(i) + (iResp)*(nCases)]
#define Residuals_(i,iResp)         Residuals      [(i) + (iResp)*(nCases)]
#define ycboSum_(iTerm,iResp)       ycboSum        [(iTerm) + (iResp)*(nMaxTerms)]
#define Betas_(iTerm,iResp)         Betas          [(iTerm) + (iResp)*(nUsedCols)]

// Global copies of some input parameters.  These stay constant for the entire MARS fit.
static double TraceGlobal;         // copy of Trace parameter
static int nMinSpanGlobal;         // copy of nMinSpan parameter
static int nEndSpanGlobal;         // copy of nEndSpan parameter
static double AdjustEndSpanGlobal; // copy of AdjustEndSpan parameter

static void FreeBetaCache(void);
static char* sFormatMemSize(const size_t MemSize, const bool Align);

//-----------------------------------------------------------------------------
// DON'T USE free, malloc, and calloc in this file.
// Use free1, malloc1, and calloc1 instead.
//
// malloc and its friends are redefined (a) so under microsoft C using
// crtdbg.h we can easily track alloc errors and (b) so FreeR() doesn't
// re-free any freed blocks and (c) so out of memory conditions are
// immediately detected.

// free1 is a macro so we can zero p (so we know it is released if FreeR is invoked)
#define free1(p)    \
{                   \
    if(p)           \
        free(p);    \
    p = NULL;       \
}

static void* malloc1(size_t size, const char* args, ...)
{
    void* p = malloc(size);
    if(!p || TraceGlobal == 1.5) {
        if(args == NULL)
            printf("malloc %s\n", sFormatMemSize(size, true));
        else {
            char s[1000];
            va_list va;
            va_start(va, args);
            vsprintf(s, args, va);
            va_end(va);
            printf("malloc %s: %s\n", sFormatMemSize(size, true), s);
        }
    }
    if(!p)
        error("Out of memory (could not allocate %s)", sFormatMemSize(size, false));
    return p;
}

static void* calloc1(size_t num, size_t size, const char* args, ...)
{
    void* p = calloc(num, size);
    if(!p || TraceGlobal == 1.5) {
        if(args == NULL)
            printf("calloc %s\n", sFormatMemSize(size, true));
        else {
            char s[1000];
            va_list va;
            va_start(va, args);
            vsprintf(s, args, va);
            va_end(va);
            printf("calloc %s: %s\n", sFormatMemSize(size, true), s);
        }
    }
    if(!p)
        error("Out of memory (could not allocate %s)", sFormatMemSize(size, false));
    return p;
}

// After calling this, on program termination we will get a report if there are
// writes outside the borders of allocated blocks or if there are non-freed blocks.

#if _MSC_VER && _DEBUG          // microsoft C and debugging enabled?
static void InitMallocTracking(void)
{
    _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_WNDW);
    _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDOUT);
    int Flag = _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG);
    Flag |= (_CRTDBG_ALLOC_MEM_DF|
             // following commented out because it makes earth run very slowly
             // _CRTDBG_CHECK_ALWAYS_DF| // call _CrtCheckMemory at every alloc and free
             _CRTDBG_LEAK_CHECK_DF|
             _CRTDBG_DELAY_FREE_MEM_DF);
    _CrtSetDbgFlag(Flag);
}
#endif

//-----------------------------------------------------------------------------
// These are malloced blocks.  They unfortunately have to be declared globally so
// under R if the user interrupts we can free them using on.exit(.C("FreeR"))

static int*    xOrder;        // local to FindTerm
static bool*   WorkingSet;    // local to FindTerm and EvalSubsets
static double* xbx;           // local to FindTerm
static double* CovSx;         // local to FindTerm
static double* CovCol;        // local to FindTerm
static double* ycboSum;       // local to FindTerm (used to be called CovSy)
static double* bxOrth;        // local to ForwardPass
static double* yMean;         // local to ForwardPass

// Transposed and mean centered copy of bxOrth, for fast update in FindKnot.
// It's faster because there is better data locality as iTerm increases, so
// better L1 cache use.  This is used only if USE_BLAS is true.

static double* bxOrthCenteredT; // local to ForwardPass

static double* bxOrthMean;      // local to ForwardPass
static int*    nDegree;         // local to Earth or ForwardPassR
static int*    nUses;           // local to Earth or ForwardPassR
#if USING_R
static int*    iDirs;           // local to ForwardPassR
static bool*   BoolFullSet;     // local to ForwardPassR
#endif
#if FAST_MARS
static void FreeQ(void);
#endif

#if USING_R
void FreeR(void)                // for use by R
{
    free1(WorkingSet);
    free1(CovSx);
    free1(CovCol);
    free1(ycboSum);
    free1(xOrder);
    free1(bxOrthMean);
    free1(bxOrthCenteredT);
    free1(bxOrth);
    free1(yMean);
    free1(BoolFullSet);
    free1(iDirs);
    free1(nUses);
    free1(nDegree);
    FreeBetaCache();
#if FAST_MARS
    FreeQ();
#endif
}
#endif

//-----------------------------------------------------------------------------
static char* sFormatMemSize(const size_t MemSize, const bool Align)
{
    static char s[100];
    double Size = (double)MemSize;
    if(Size >= 1e9)
        sprintf(s, Align? "%6.3f GB": "%.3g GB", Size / ((size_t)1 << 30));
    else if(Size >= 1e6)
        sprintf(s, Align? "%6.0f MB": "%.3g MB", Size / ((size_t)1 << 20));
    else if(Size >= 1e3)
        sprintf(s, Align? "%6.0f kB": "%.3g kB", Size / ((size_t)1 << 10));
    else
        sprintf(s, Align? "%6.0f  B": "%g B", Size);
    return s;
}

//-----------------------------------------------------------------------------
static void tprintf(const int trace, const char *args, ...) // printf with trace check
{
    if(TraceGlobal >= trace) {
        char s[1000];
        va_list va;
        va_start(va, args);
        vsprintf(s, args, va);
        va_end(va);
        printf("%s", s);
    }
}

//-----------------------------------------------------------------------------
// Gets called periodically to service the R framework.
// Will not return if the user interrupts.

#if USING_R

static INLINE void ServiceR(void)
{
    R_FlushConsole();
    R_CheckUserInterrupt();     // may never return
}

#endif

//-----------------------------------------------------------------------------
#if FAST_MARS

typedef struct tQueue {
    int     iParent;            // parent term
    double  RssDelta;
    int     nTermsForRssDelta;  // number of terms when RssDelta was calculated
    double  AgedRank;
} tQueue;

static tQueue* Q;       // indexed on iTerm (this Q is used for queue updates)
static tQueue* SortedQ; // indexed on iParent rank (this Q is used to get next iParent)
static int     nQMax;   // number of elements in Q

static void InitQ(const int nMaxTerms)
{
    nQMax = 0;
    Q       = (tQueue*)malloc1(nMaxTerms * sizeof(tQueue),
                            "Q\t\t\tnMaxTerms %d sizeof(tQueue) %d",
                            nMaxTerms, sizeof(tQueue));
    SortedQ = (tQueue*)malloc1(nMaxTerms * sizeof(tQueue),
                            "SortedQ\t\tnMaxTerms %d sizeof(tQueue) %d",
                            nMaxTerms, sizeof(tQueue));
    for(int iTerm = 0; iTerm < nMaxTerms; iTerm++) {
        Q[iTerm].iParent = iTerm;
        Q[iTerm].nTermsForRssDelta = -99;   // not strictly needed, nice for debugging
        Q[iTerm].RssDelta = -1;
        Q[iTerm].AgedRank = -1;
    }
}

static void FreeQ(void)
{
    free1(Q);
    free1(SortedQ);
}

static void PrintSortedQ(int nFastK)     // for debugging
{
    printf("\n\nSortedQ  QIndex Parent nTermsForRssDelta AgedRank  RssDelta\n");
    for(int i = 0; i < nQMax; i++) {
        printf("            %3d    %3d   %15d    %5.1f  %g\n",
            i+IOFFSET,
            SortedQ[i].iParent+IOFFSET,
            SortedQ[i].nTermsForRssDelta+IOFFSET,
            SortedQ[i].AgedRank,
            SortedQ[i].RssDelta);
        if(i == nFastK-1)
            printf("FastK %d ----------------------------------------------------\n",
                nFastK);
    }
}

// Sort so highest RssDeltas are at low indices.
// Secondary sort key is iParent.  Not strictly needed, but removes
// possible differences in qsort implementations (which "sort"
// identical keys unpredictably).

static int CompareQ(const void* p1, const void* p2)     // for qsort
{
    double Diff = ((tQueue*)p2)->RssDelta - ((tQueue*)p1)->RssDelta;
    if(Diff < 0)
        return -1;
    else if(Diff > 0)
        return 1;

    // Diff is 0, so sort now on iParent

    int iDiff = ((tQueue*)p1)->iParent - ((tQueue*)p2)->iParent;
    if(iDiff < 0)
        return -1;
    else if(iDiff > 0)
        return 1;
    return 0;
}

// Sort so lowest AgedRanks are at low indices.
// If AgedRanks are the same then sort on RssDelta and iParent.

static int CompareAgedQ(const void* p1, const void* p2) // for qsort
{
    double Diff = ((tQueue*)p1)->AgedRank - ((tQueue*)p2)->AgedRank;
    if(Diff < 0)
        return -1;
    else if(Diff > 0)
        return 1;

    // Diff is 0, so sort now on RssDelta

    Diff = ((tQueue*)p2)->RssDelta - ((tQueue*)p1)->RssDelta;
    if(Diff < 0)
        return -1;
    else if(Diff > 0)
        return 1;

    // Diff is still 0, so sort now on iParent

    int iDiff = ((tQueue*)p1)->iParent - ((tQueue*)p2)->iParent;
    if(iDiff < 0)
        return -1;
    else if(iDiff > 0)
        return 1;
    return 0;
}

static void AddTermToQ(
    const int iTerm,        // in
    const int nTerms,       // in
    const double RssDelta,  // in
    const bool Sort,        // in
    const int nMaxTerms,    // in
    const double FastBeta)  // in: ageing Coef, 0 is no ageing, FastMARS recommends 1
{
    ASSERT(iTerm < nMaxTerms);
    ASSERT(nQMax < nMaxTerms);
    Q[nQMax].nTermsForRssDelta = nTerms;
    Q[nQMax].RssDelta = max(Q[iTerm].RssDelta, RssDelta);
    nQMax++;
    if(Sort) {
        memcpy(SortedQ, Q, nQMax * sizeof(tQueue));
        qsort(SortedQ, nQMax, sizeof(tQueue), CompareQ);         // sort on RssDelta
        if(FastBeta > 0) {
            for(int iRank = 0; iRank < nQMax; iRank++)
                SortedQ[iRank].AgedRank =
                    iRank + FastBeta * (nTerms - SortedQ[iRank].nTermsForRssDelta);
            qsort(SortedQ, nQMax, sizeof(tQueue), CompareAgedQ); // sort on aged rank
        }
    }
}

static void UpdateRssDeltaInQ(const int iParent, const int nTermsForRssDelta,
                              const double RssDelta)
{
    ASSERT(iParent == Q[iParent].iParent);
    ASSERT(iParent < nQMax);
    Q[iParent].nTermsForRssDelta = nTermsForRssDelta;
    Q[iParent].RssDelta = RssDelta;
}

static int GetNextParent(   // returns -1 if no more parents
    const bool InitFlag,    // use true to init, thereafter false
    const int  nFastK)
{
    static int iQ;          // index into sorted queue
    int iParent = -1;
    if(InitFlag) {
        if(TraceGlobal == 6)
            printf("\n|Considering parents ");
        iQ = 0;
    } else {
        if(iQ < min(nQMax, nFastK)) {
            iParent = SortedQ[iQ].iParent;
            iQ++;
        }
        if(TraceGlobal == 6 && iParent >= 0)
            printf("%d [%g] ", iParent+IOFFSET, SortedQ[iQ].RssDelta);
    }
    return iParent;
}

#endif // FAST_MARS

//-----------------------------------------------------------------------------
// This helps reduce the effects of numerical err, mostly when testing.

static INLINE double MaybeZero(double x)
{
    return x < ALMOST_ZERO? 0. : x;
}

//-----------------------------------------------------------------------------
// GetOrder() gets the sort indices of vector x, so
// x[sorted[i]] <= x[sorted[i+1]].  Ties may be reordered. The returned
// indices are 0 based (as in C not as in R).
//
// This function is similar to the R library function rsort_with_index(),
// but is defined here to minimize R dependencies.
// Informal tests show that this is faster than rsort_with_index().

static const double* pxGlobal;

static int Compare(const void* p1, const void* p2)  // for qsort
{
    const int i1 = *(int*)p1;
    const int i2 = *(int*)p2;
    double Diff = pxGlobal[i1] - pxGlobal[i2];
    if(Diff < 0)
        return -1;
    else if(Diff > 0)
        return 1;
    else
        return 0;
}

static void GetOrder(
    int sorted[],      // out: vec with nx elements
    const double x[],  // in: x is a vec with nx elems
    const int nx)      // in: number of elems in x
{
    for(int i = 0; i < nx; i++)
        sorted[i] = i;
    pxGlobal = x;
    qsort(sorted, nx, sizeof(int), Compare);
}


//-----------------------------------------------------------------------------
// Get order indices for an x array of dimensions nRows x nCols.
//
// Returns an nRows x nCols integer array of indices, where each column
// corresponds to a column of x.  See GetOrder() for ordering details.
//
// Caller must free the returned array.

static int* GetArrayOrder(
    const double x[],   // in
    const size_t nRows, // in
    const int    nCols) // in
{
    int* xOrder = (int*)malloc1(nRows * nCols * sizeof(int),
                            "xOrder\t\tnRows %d nCols %d sizeof(int) %d",
                            nRows, nCols, sizeof(int));

    for(int iCol = 0; iCol < nCols; iCol++) {
        GetOrder(xOrder + iCol*nRows, x + iCol*nRows, (int)nRows);
#if USING_R
        if(nRows > (int)1e4)
            ServiceR();
#endif
    }
    return xOrder;
}

//-----------------------------------------------------------------------------
// return the number of TRUEs in the boolean vector UsedCols

static int GetNbrUsedCols(const bool UsedCols[], const int nLen)
{
    int nTrue = 0;

    for(int iCol = 0; iCol < nLen; iCol++)
        if(UsedCols[iCol])
            nTrue++;

    return nTrue;
}

//-----------------------------------------------------------------------------
// Copy used columns in x to *pxUsed and return the number of used columns
// UsedCols[i] is true for each each used column index in x
// Caller must free *pxUsed

static int CopyUsedCols(
    double** pxUsed,        // out: caller must free
    const double x[],       // in: nCases x nCols
    const size_t nCases,    // in
    const int nCols,        // in
    const bool UsedCols[])  // in
{
    const int nUsedCols = GetNbrUsedCols(UsedCols, nCols);
    double* xUsed = (double*)malloc1(nCases * nUsedCols * sizeof(double),
                        "xUsed\t\t\tnCases %d nUsedCols %d sizeof(double) %d",
                        (int)nCases, nUsedCols, sizeof(double));

    int iUsed = 0;
    for(int iCol = 0; iCol < nCols; iCol++)
        if(UsedCols[iCol]) {
            memcpy(xUsed + iUsed * nCases,
                x + iCol * nCases, nCases * sizeof(double));
            iUsed++;
        }
    *pxUsed = xUsed;
    return nUsedCols;
}

//-----------------------------------------------------------------------------
// Print a summary of the model, for debug tracing

#if STANDALONE
static void PrintSummary(
    const int    nMaxTerms,   // in
    const int    nTerms,      // in: number of cols in bx, some may be unused
    const int    nPreds,      // in: number of predictors
    const int    nResp,       // in: number of cols in y
    const bool   UsedCols[],  // in: specifies used columns in bx
    const int    Dirs[],      // in
    const double Cuts[],      // in
    const double Betas[],     // in: if NULL will print zeros
    const int    nDegree[])   // in: degree of each term, degree of intercept is 0
{
    printf("   nFacs       Beta\n");

    const int nUsedCols = GetNbrUsedCols(UsedCols, nTerms);
    int iUsed = -1;
    for(int iTerm = 0; iTerm < nTerms; iTerm++) {
        if(UsedCols[iTerm]) {
            iUsed++;
            printf("%2.2d  %2d    ", iTerm, nDegree[iTerm]);
            for(int iResp = 0; iResp < nResp; iResp++)
                printf("%9.3g ", (Betas? Betas_(iUsed, iResp): 0));
            printf("| ");
            }
        else {
            printf("%2.2d  --    ", iTerm);
            for(int iResp = 0; iResp < nResp; iResp++)
                printf("%9s ", "--");
            printf("| ");
        }
        int iPred;
        for(iPred = 0; iPred < nPreds; iPred++)
            if(Dirs_(iTerm,iPred) == 0)
                printf(" . ");
            else
                printf("%2d ", Dirs_(iTerm,iPred));

        printf("|");

        for(iPred = 0; iPred < nPreds; iPred++)
            if(Dirs_(iTerm,iPred) == 0)
                printf("    .    ");
            else if(Dirs_(iTerm,iPred) == 2)
                printf("  linear ");
            else
                printf("%8.3g ", Cuts_(iTerm,iPred));

        printf("\n");
    }
    printf("\n");
}
#endif // STANDALONE

//-----------------------------------------------------------------------------
// Set Diags to the diagonal values of inverse(X'X),
// where X is referenced via the matrix R, from a previous call to dqrdc2
// with (in practice) bx.  The net result is that Diags is the diagonal
// values of inverse(bx'bx).  We assume that R is created from a full rank X.
//
// TODO This could be simplified

static void CalcDiags(
    double Diags[],      // out: nCols x 1
    const double R[],    // in: nCases x nCols, QR from prev call to dqrdc2
    const size_t nCases, // in
    const int nCols)     // in
{
    #define R_(i,j)     R [(i) + (j) * nCases]
    #define R1_(i,j)    R1[(i) + (const size_t)(j) * nCols]
    #define B_(i,j)     B [(i) + (const size_t)(j) * nCols]

    double* R1 = (double*)malloc1(nCols * nCols * sizeof(double),  // nCols rows of R
                            "R1\t\t\tnCols %d nCols %d sizeof(double) %d",
                            nCols, nCols, sizeof(double));

    double* B =  (double*)calloc1(nCols * nCols, sizeof(double),   // rhs of R1 * x = B
                            "B\t\t\tnCols %d nCols %d sizeof(double) %d",
                            nCols, nCols, sizeof(double));
    int i, j;
    for(i = 0; i < nCols; i++) {    // copy nCols rows of R into R1
        for(j =  0; j < nCols; j++)
            R1_(i,j) = R_(i,j);
        B_(i,i) = 1;                // set diag of B to 1
    }
    int job = 1;            // 1 means solve R1 * x = B where R1 is upper triangular
    int info = 0;
    for(i = 0; i < nCols; i++) {
        dtrsl_(             // LINPACK function
            R1,             // in: t, matrix of the system, untouched
            (int*)&nCols,   // in: ldt (typecast discards const)
            (int*)&nCols,   // in: n
            &B_(0,i),       // io: b, on return has solution x
            &job,           // in:
            &info);         // io:

        ASSERT(info == 0);
    }
    // B is now inverse(R1).  Calculate B x B.

    for(i = 0; i < nCols; i++)
        for(j =  0; j < nCols; j++) {
            double Sum = 0;
            for(int k = max(i,j); k < nCols; k++)
                Sum += B_(i,k) * B_(j,k);
            B_(i,j) = B_(j,i) = Sum;
        }
    for(i = 0; i < nCols; i++)
         Diags[i] = B_(i,i);
    free1(B);
    free1(R1);
}

//-----------------------------------------------------------------------------
// Regress y on the used columns of x, in the standard way (using QR).
// UsedCols[i] is true for each each used col i in x; unused cols are ignored.
//
// The returned Betas argument is computed from, and is indexed on,
// the compacted x vector, not on the original x.
//
// The returned iPivots should only be used if *pnRank != nUsedCols.
// The entries of iPivots refer to columns in the full x (and are 0 based).
// Entries in iPivots at *pnRank and above specify linearly dependent columns in x.
//
// To maximize compatibility we call the same routines as the R function lm.

static void Regress(
    double       Betas[],       // out: nUsedCols * nResp, can be NULL
    double       Residuals[],   // out: nCases * nResp, can be NULL
    double*      pRss,          // out: RSS, summed over all nResp, can be NULL
    double       Diags[],       // out: diags of inv(transpose(x) * x), can be NULL
    int*         pnRank,        // out: nbr of indep cols in x, can be NULL
    int          iPivots[],     // out: nCols, can be NULL
    const double x[],           // in: nCases x nCols, must include intercept
    const double y[],           // in: nCases x nResp
    const size_t nCases,        // in: number of rows in x and in y
    const int    nResp,         // in: number of cols in y
    const int    nCols,         // in: number of columns in x, some may not be used
    const bool   UsedCols[])    // in: specifies used columns in x
{
    double* xUsed;
    int nUsedCols = CopyUsedCols(&xUsed, x, nCases, nCols, UsedCols);
    bool MustFreeResiduals = false;
    if(Residuals == NULL) {
        Residuals = (double*)malloc1(nCases * nResp * sizeof(double),
                                "Residuals\t\tnCases %d nResp %d sizeof(double) %d",
                                (int)nCases, nResp, sizeof(double));
        MustFreeResiduals = true;
    }
    bool MustFreePivots = false;
    if(iPivots == NULL) {
        iPivots = (int*)malloc1(nUsedCols * sizeof(int),
                            "iPivots\t\tnUsedCols %d sizeof(int) %d",
                            nUsedCols, sizeof(int));
        MustFreePivots = true;
    }
    int iCol;
    for(iCol = 0; iCol < nUsedCols; iCol++)
        iPivots[iCol] = iCol+1;

    double* qraux = (double*)malloc1(nUsedCols * sizeof(double),
                                "qraux\t\t\tnUsedCols %d sizeof(double) %d",
                                nUsedCols, sizeof(double));

    // work size must be nUsedCols*2 for dqrdc2, and nCases*nUsedCols for
    // dqrsl where it is used for qy, qty, and rsd

    double* work = (double*)malloc1(
                    max((size_t)nUsedCols * 2, nCases * nUsedCols) * sizeof(double),
                    "work\t\t\tnCases %d nUsedCols %d sizeof(double) %d",
                    (int)nCases, nUsedCols, sizeof(double));

    int nCases1 = (int)nCases; // type convert from size_t

    int nRank;
    dqrdc2_(                // R function, QR decomp based on LINPACK dqrdc
        xUsed,              // io:  x, on return upper tri of x is R of QR
        &nCases1,           // in:  ldx
        &nCases1,           // in:  n
        &nUsedCols,         // in:  p
        (double*)&QR_TOL,   // in:  tol
        &nRank,             // out: k, num of indep cols of x
        qraux,              // out: qraux
        iPivots,            // out: jpvt
        work);              // work

    double Rss = 0;
    const bool NeedResiduals = !MustFreeResiduals || pRss;

    int job = (Betas? 100: 0) + (NeedResiduals? 10: 0);

    if(job) { // job will be zero if all we need are the iPivots from dqrdc2
        for(int iResp = 0; iResp < nResp; iResp++) {
            int info;
            dqrsl_(                 // LINPACK function
                xUsed,              // in:  x, generated by dqrdc2
                &nCases1,           // in:  ldx
                &nCases1,           // in:  n
                &nRank,             // in:  k
                qraux,              // in:  qraux
                (double*)(y + iResp * nCases),  // in:  y
                NULL,               // out: qy, unused here
                work,               // out: qty, required if rsd in job
                Betas?              // out: b, only needed if user asked for them
                    (double*)(&Betas_(0,iResp)): work,
                NeedResiduals?      // out: rsd
                    (double*)(&Residuals_(0,iResp)): NULL,
                NULL,               // out: xb = yHat, unused here
                &job,               // in:  job
                &info);             // in:  info

            ASSERT(info == 0);

            // compute Residuals and Rss (sum over all responses)

            if(NeedResiduals)
                for(int i = 0; i < (const int)nCases; i++)
                    Rss += sq(Residuals_(i, iResp));
        }
        if(pRss)
            *pRss = Rss;
    }
    if(nRank != nUsedCols &&
        !MustFreePivots ) { // only bother if caller wants iPivots back
        // adjust iPivots for missing cols in UsedCols and for 1 offset

        int* PivotOffset = (int*)malloc1(nCols * sizeof(int),
                                    "PivotOffset\t\t\tnCols %d sizeof(int) %d",
                                    nCols, sizeof(int));
        int nOffset = 0, iOld = 0;
        for(iCol = 0; iCol < nCols; iCol++) {
            if(!UsedCols[iCol])
                nOffset++;
            else {
                PivotOffset[iOld] = nOffset;
                if(++iOld > nUsedCols)
                    break;
            }
        }
        for(iCol = 0; iCol < nUsedCols; iCol++)
            iPivots[iCol] = iPivots[iCol] - 1 + PivotOffset[iPivots[iCol] - 1];
        free1(PivotOffset);
    }
    if(pnRank)
        *pnRank = nRank;
    if(Diags)
        CalcDiags(Diags, xUsed, nCases, nUsedCols);
    if(MustFreePivots)
        free1(iPivots);
    if(MustFreeResiduals)
        free1(Residuals);
    free1(xUsed);
    free1(qraux);
    free1(work);
}

//-----------------------------------------------------------------------------
// This routine is for testing Regress from R, to compare results to R's lm().

#if USING_R
void RegressR(                  // for testing earth routine Regress from R
    double       Betas[],       // out: (nUsedCols+1) * nResp, +1 is for intercept
    double       Residuals[],   // out: nCases * nResp
    double       Rss[],         // out: RSS, summed over all nResp
    double       Diags[],       // out: diags of inv(transpose(x) * x)
    int*         pnRank,        // out: nbr of indep cols in x
    int          iPivots[],     // out: nCols
    const double x[],           // in: nCases x nCols
    const double y[],           // in: nCases x nResp
    const int*   pnCases,       // in: number of rows in x and in y
    const int*   pnResp,        // in: number of cols in y
    int*         pnCols,        // in: number of columns in x, some may not be used
    const bool   UsedCols[])    // in: specifies used columns in x
{
    const size_t nCases1 = *pnCases; // type convert

    Regress(Betas, Residuals, Rss, Diags, pnRank, iPivots,
        x, y, nCases1, *pnResp, *pnCols, UsedCols);
}
#endif

//-----------------------------------------------------------------------------
// Regress y on bx to get Residuals and Betas.  If bx isn't of full rank,
// remove dependent cols, update UsedCols, and regress again on the bx with
// removed cols.

static void RegressAndFix(
    double       Betas[],      // out: nMaxTerms x nResp, can be NULL
    double       Residuals[],  // out: nCases x nResp, can be NULL
    double       Diags[],      // out: if !NULL set to diags of inv(transpose(bx) * bx)
    bool         UsedCols[],   // io:  will remove cols if necessary, nMaxTerms x 1
    const        double bx[],  // in:  nCases x nMaxTerms
    const double y[],          // in:  nCases x nResp
    const size_t nCases,       // in
    const int    nResp,        // in: number of cols in y
    const int    nTerms)       // in: number of cols in bx, some may not be used
{
    int nRank;
    int* iPivots = (int*)malloc1(nTerms * sizeof(int),
                            "iPivots\t\tnTerms %d sizeof(int) %d",
                            nTerms, sizeof(int));
    Regress(Betas, Residuals, NULL, Diags, &nRank, iPivots,
        bx, y, nCases, nResp, nTerms, UsedCols);
    int nUsedCols = GetNbrUsedCols(UsedCols, nTerms);
    const int nDeficient = nUsedCols - nRank;
    if(nDeficient) {    // rank deficient?
        // Remove linearly dependent columns.
        // The lin dep columns are at index nRank and higher in iPivots.

        for(int iCol = nRank; iCol < nUsedCols; iCol++)
            UsedCols[iPivots[iCol]] = false;

        Regress(Betas, Residuals, NULL, Diags, &nRank, NULL,
            bx, y, nCases, nResp, nTerms, UsedCols);
        nUsedCols = nUsedCols - nDeficient;
        if(nRank != nUsedCols)
            warning("Could not fix rank deficient bx: nUsedCols %d nRank %d",
                nUsedCols,  nRank);
        else {
            tprintf(1,
                "Fixed rank deficient bx by removing %d term%s, %d term%s remain%s\n",
                nDeficient, ((nDeficient==1)? "": "s"),
                nUsedCols,  ((nUsedCols==1)? "": "s"), ((nUsedCols==1)? "s": ""));
            tprintf(4, "\n");
        }
    }
    free1(iPivots);
}

//-----------------------------------------------------------------------------
static INLINE double Mean(const double x[], size_t n)
{
    double mean = 0;
    for(size_t i = 0; i < n; i++)
        mean += x[i] / n;
    return mean;
}

//-----------------------------------------------------------------------------
// get mean centered sum of squares

static INLINE double SumOfSquares(const double x[], const double mean, size_t n)
{
    double ss = 0;
    for(size_t i = 0; i < n; i++)
        ss += sq(x[i] - mean);
    return ss;
}

//-----------------------------------------------------------------------------
static INLINE double GetGcv(const int nTerms, // nbr basis terms including intercept
                const size_t nCases, double Rss, const double Penalty)
{
    double Cost;
    if(Penalty == -1)   // special case: terms and knots are free
        Cost = 0;
    else {
        const double nKnots = ((double)nTerms-1) / 2;
        Cost = (nTerms + Penalty * nKnots) / nCases;
    }
    // test against Cost ensures that GCVs are non-decreasing as nbr of terms increases
    return Cost >= 1? POS_INF : Rss / (nCases * sq(1 - Cost));
}

//-----------------------------------------------------------------------------
// Check if model term type is already in model, to avoid a linear dependence.
// TODO The code in this routine doesn't seem to make sense.

static bool GetNewFormFlag(const int iPred, const int iTerm,
                        const int Dirs[], const bool UsedCols[],
                        const int nTerms, const int nPreds, const int nMaxTerms)
{
    bool IsNewForm = true;
    for(int i = 1; i < nTerms; i++) // start at 1 to skip intercept
        if(UsedCols[i]) {
            IsNewForm = false;
            if(Dirs_(i,iPred) == 0) // unused in term i
                return true;
            // TODO if the following code is commented out, test suite passes (!)
            for(int j = 0; j < nPreds; j++)
                if(j != iPred && (Dirs_(i,j) != 0) != (Dirs_(iTerm,j) != 0))
                    return true;
        }
    return IsNewForm;
}

//-----------------------------------------------------------------------------
static double GetCut(
    const int i,
    const int iPred,
    const size_t nCases,
    const double x[],
    const int xOrder[])
{
    if(i < 0 || i >= (const int)nCases)
        error("GetCut i %d: i < 0 || i >= nCases", i);
    const int ix = xOrder_(i,iPred);
    if(ix < 0 || ix >= (const int)nCases)
        error("GetCut ix %d: ix < 0 || ix >= nCases", ix);
    return x_(ix,iPred);
}

//-----------------------------------------------------------------------------
// The BetaCache is used when searching for a new term pair, via FindTerm.
// Most of the calculation for the orthogonal regression betas is repeated
// with the same data, and thus we can save time by caching betas.
// (The "Betas" are the regression coefficients.)
//
// iParent    is the term that forms the base for the new term
// iPred      is the predictor for the new term
// iOrthCol   is the column index in the bxOrth matrix

static double* BetaCacheGlobal; // [iOrthCol,iParent,iPred]
                                // dim nPreds x nMaxTerms x nMaxTerms

static void InitBetaCache(const bool UseBetaCache,
                          const int nMaxTerms, const int nPreds)
{
    int nCache =  nMaxTerms * nMaxTerms * nPreds;
    if(!UseBetaCache) {
        BetaCacheGlobal = NULL;
    // 3e9 below is somewhat arbitrary but seems about right (in 2011)
    } else if(nCache * sizeof(double) > 3e9) {
            printf(
"\nNote: earth's beta cache would require %s, so forcing Use.beta.cache=FALSE.\n"
"      Invoke earth with Use.beta.cache=FALSE to make this message go away.\n\n",
                sFormatMemSize(nCache * sizeof(double), false));
            BetaCacheGlobal = NULL;
    } else {
        tprintf(5, "BetaCache %s\n", // print cache size
                sFormatMemSize(nCache * sizeof(double), false));

        BetaCacheGlobal = (double*)malloc1(nCache * sizeof(double),
            "BetaCacheGlobal\tnMaxTerms %d nMaxTerms %d nPreds %d sizeof(double) %d",
            nMaxTerms, nMaxTerms, nPreds, sizeof(double));

        for(int i = 0; i < nCache; i++) // mark all entries as uninitialized
            BetaCacheGlobal[i] = POS_INF;
    }
}

static void FreeBetaCache(void)
{
    if(BetaCacheGlobal)
        free1(BetaCacheGlobal);
}

//-----------------------------------------------------------------------------
// Init a new bxOrthCol to the residuals from regressing y on the used columns
// of the orthogonal matrix bxOrth.  The length (i.e. sum of squares divided
// by nCases) of each column of bxOrth must be 1 with mean 0 (except the
// first column which is the intercept).
//
// In practice this function is called with the params shown in {braces}
// and is called only by InitBxOrthCol.
//
// This function must be fast.
//
// In calculation of Beta, we used to have
//     xty += pbxOrth[i] * y[i];
// and now we have
//    xty += pbxOrth[i] * bxOrthCol[i];
// i.e. we use the "modified" instead of the "classic" Gram Schmidt.
// This is less susceptible to numerical error, although it is rare
// to see the effect in practice in earth models (but you can see it in the
// final model in the test suite function test.zigzag in test.weights.R).

static INLINE void OrthogResiduals(
    double bxOrthCol[],     // out: nCases x 1      { bxOrth[,nTerms] }
    const double y[],       // in:  nCases x nResp  { bx[,nTerms], xbx }
    const double bxOrth[],  // in:  nTerms x nPreds { bxOrth }
    const size_t nCases,    // in
    const int nTerms,       // in: nTerms in model, i.e. number of used cols in bxOrth
    const bool UsedTerms[], // in: UsedTerms[i] is true if col is used, unused cols ignored
                            //     Following parameters are only for the beta cache
    const int iParent,      // in: if >= 0, use BetaCacheGlobal {FindTerm iTerm, addTermP -1}
    const int iPred,        // in: predictor index i.e. col index in input matrix x
    const int nMaxTerms)    // in:
{
    double* pCache;
    if(iParent >= 0 && BetaCacheGlobal)
        pCache = BetaCacheGlobal + iParent*nMaxTerms + iPred*sq(nMaxTerms);
    else
        pCache = NULL;

    memcpy(bxOrthCol, y, nCases * sizeof(double));

    for(int iTerm = 0; iTerm < nTerms; iTerm++)
        if(UsedTerms[iTerm]) {
            const double* pbxOrth = &bxOrth_(0, iTerm);
            double Beta;
            if(pCache && pCache[iTerm] != POS_INF)
                Beta = pCache[iTerm];
            else {
                double xty = 0;
                for(int i = 0; i < (const int)nCases; i++)
                    xty += pbxOrth[i] * bxOrthCol[i]; // see header comment
                Beta = xty;  // no need to divide by xtx, it is 1
                ASSERT(FINITE(Beta));
                if(pCache)
                    pCache[iTerm] = Beta;
            }
            if(USE_BLAS) {
                const double NegBeta = -Beta;
                const int nCases1 = (int)nCases; // type convert from size_t
                daxpy_(&nCases1, &NegBeta, pbxOrth, &ONE, bxOrthCol, &ONE);
            } else for(int i = 0; i < (const int)nCases; i++)
                bxOrthCol[i] -= Beta * pbxOrth[i];
        }
}

//-----------------------------------------------------------------------------
// Init the rightmost column of bxOrth i.e. the column indexed by nTerms.
// The new col is the normalized residuals from regressing y on the
// lower (i.e. already existing) cols of bxOrth.
// Also updates bxOrthCenteredT and bxOrthMean.
//
// In practice this function is called only with the params shown in {braces}

static INLINE void InitBxOrthCol(
    double bxOrth[],          // io: col nTerms is changed, other cols not touched
    double bxOrthCenteredT[], // io: kept in sync with bxOrth
    double bxOrthMean[],      // io: element at nTerms is updated
    bool* pGoodCol,           // io: true if col sum-of-squares is greater than MIN_BX_SOS
    const double* y,          // in: { AddCandLinTerm xbx, AddTermPair bx[,nTerms] }
    const int nTerms,         // in: column goes in at index nTerms, 0 is the intercept
    const bool WorkingSet[],  // in
    const size_t nCases,      // in
    const int nMaxTerms,      // in
    const int iCacheTerm,     // in: if >= 0, use BetaCacheGlobal {FindTerm iTerm, AddTermP -1}
                              //     if < 0 then recalc Betas from scratch
    const int iPred)          // in: predictor index i.e. col index in input matrix x
{
    *pGoodCol = true;
    int i;

    if(nTerms == 0) {           // column 0, the intercept
        double len = 1 / sqrt((double)nCases);
        for(i = 0; i < (const int)nCases; i++)
            bxOrth_(i,0) = len;
        bxOrthMean[0] = len;
    } else if(nTerms == 1) {    // column 1, the first basis function, y = xbx = x[,1]
        double yMean = Mean(y, nCases);
        for(i = 0; i < (const int)nCases; i++)
            bxOrth_(i,1) = y[i] - yMean;
    } else
        OrthogResiduals(&bxOrth_(0,nTerms), // resids go in rightmost col of bxOrth at nTerms
            y, bxOrth, nCases, nTerms, WorkingSet, iCacheTerm, iPred, nMaxTerms);

    if(nTerms > 0) {
        // normalize the column to length 1 and init bxOrthMean[nTerms]

        double bxOrthSS = SumOfSquares(&bxOrth_(0,nTerms), 0, nCases);
        if(bxOrthSS <= MIN_BX_SOS)
            *pGoodCol = false;
        // iCacheTerm will be negative unless called by AddCandidateLinearTerm
        if(bxOrthSS > (iCacheTerm < 0? 0: MIN_BX_SOS)) {
            bxOrthMean[nTerms] = Mean(&bxOrth_(0,nTerms), nCases);
            const double len = sqrt(bxOrthSS);
            for(i = 0; i < (const int)nCases; i++)
                bxOrth_(i,nTerms) /= len;
        } else {
            bxOrthMean[nTerms] = 0;
            memset(&bxOrth_(0,nTerms), 0, nCases * sizeof(double));
        }
    }
    for(i = 0; i < (const int)nCases; i++) // keep bxOrthCenteredT in sync
        bxOrthCenteredT_(nTerms,i) = bxOrth_(i,nTerms) - bxOrthMean[nTerms];
}

//-----------------------------------------------------------------------------
// Add a new term pair to the arrays.
// Each term in the new term pair is a copy of an existing parent term but extended
// by multiplying it by a new hinge function at the selected knot.
// If the upper term in the term pair is invalid then we still add the upper
// term but mark it as false in FullSet.

static void AddTermPair(
    int    Dirs[],              // io
    double Cuts[],              // io
    double bx[],                // io: MARS basis matrix
    double bxOrth[],            // io
    double bxOrthCenteredT[],   // io
    double bxOrthMean[],        // io
    bool   FullSet[],           // io
    bool*  pIsNewForm,          // io
    int    nDegree[],           // io: degree of each term, degree of intercept is 0
    int    nUses[],             // io: nbr of times each predictor is used in the model
    const int nTerms,           // in: new term pair goes in at index nTerms and nTerms1
    const int iBestParent,      // in: parent term
    const int iBestCase,        // in
    const int iBestPred,        // in
    const int nPreds,           // in
    const size_t nCases,        // in
    const int nMaxTerms,        // in
    const bool LinPredIsBest,   // in: true if pred should enter linearly (no knot)
    const int LinPreds[],       // in: user specified preds which must enter linearly
    const double x[],           // in
    const int xOrder[],         // in
    const bool Weighted)        // in
{
    const int nTerms1 = nTerms+1;

    // copy the parent term to the new term pair

    int iPred;
    for(iPred = 0; iPred < nPreds; iPred++) {
        Dirs_(nTerms, iPred) = Dirs_(nTerms1,iPred) = Dirs_(iBestParent,iPred);
        Cuts_(nTerms, iPred) = Cuts_(nTerms1,iPred) = Cuts_(iBestParent,iPred);
    }
    // incorporate the new hinge function

    nDegree[nTerms]  = nDegree[nTerms1] = nDegree[iBestParent] + 1;

    int DirEntry = 1;
    if(LinPreds[iBestPred] || LinPredIsBest) { // changed in earth 4.0.0
        ASSERT(LinPredIsBest);
        DirEntry = 2;
    }
    Dirs_(nTerms, iBestPred) = DirEntry;
    Dirs_(nTerms1,iBestPred) = -1; // will be ignored if adding only one hinge

    const double BestCut = GetCut(iBestCase, iBestPred, nCases, x, xOrder);

    Cuts_(nTerms, iBestPred) = Cuts_(nTerms1,iBestPred) = BestCut;

    // Fill in new columns of bx, at nTerms and nTerms+1 (left and right hinges).
#if WEIGHTS
    // in FindWeightedPredGivenParent, we used the two columns
    // in bx as a scratch buffer, so zero them again
    ASSERT(nTerms+1 < nMaxTerms);
    memset(&bx_(0,nTerms),   0, nCases * sizeof(double));
    memset(&bx_(0,nTerms+1), 0, nCases * sizeof(double));
#endif
    int i;
    if(DirEntry == 2) { // linpred?
        for(i = 0; i < (const int)nCases; i++)      // add single term
            bx_(i,nTerms) = bx_(i,iBestParent) * x_(i,iBestPred);
    } else for(i = 0; i < (const int)nCases; i++) { // add term pair
        const int iOrdered = xOrder_(i, iBestPred);
        const double xi = x_(iOrdered, iBestPred);
        if(i > iBestCase)
            bx_(iOrdered, nTerms) = bx_(iOrdered, iBestParent) * (xi - BestCut);
        else
            bx_(iOrdered, nTerms1) = bx_(iOrdered, iBestParent) * (BestCut - xi);
    }
    nUses[iBestPred]++;

    // init the col in bxOrth at nTerms and init bxOrthMean[nTerms]

    FullSet[nTerms] = true;
    bool GoodCol;
    InitBxOrthCol(bxOrth, bxOrthCenteredT, bxOrthMean, &GoodCol,
        &bx_(0,nTerms), nTerms, FullSet, nCases, nMaxTerms, -1, nPreds);
                            // -1 means don't use BetaCacheGlobal, calc Betas afresh
#if 0 // TODO commented out because this sometimes happens
    if(!GoodCol)
        printf("GoodCol is false in AddTermPair\n");
#endif
    // init the col in bxOrth at nTerms1 and init bxOrthMean[nTerms1]

    if(!LinPredIsBest && *pIsNewForm)
        FullSet[nTerms1] = true;

    if(FullSet[nTerms1]) {
        InitBxOrthCol(bxOrth, bxOrthCenteredT, bxOrthMean, &GoodCol,
            &bx_(0,nTerms1), nTerms1, FullSet, nCases, nMaxTerms, -1, iPred);

        if(Weighted && !GoodCol) {
            // !GoodCol usuall happens only when weights are used.  For the
            // nonweighted situation we would have already cleared IsNewForm in
            // FindPredGivenParent (inside AddCandidateLinearTerm), although
            // we sometimes get a non GoodCol here without weights.
            tprintf(6, "clear IsNewForm\n");
            *pIsNewForm = false;
            FullSet[nTerms1] = false;
        }
    }
    if(!FullSet[nTerms1]) {
        // If the term is not valid, then we don't wan't to use it as the
        // base for a new term later (in FindTerm).  Enforce this by setting
        // nDegree to a value greater than any posssible nMaxDegree.
        nDegree[nTerms1] = MAX_DEGREE + 1;
        memset(&bxOrth_(0,nTerms1), 0, nCases * sizeof(double));
        bxOrthMean[nTerms1] = 0;
        for(i = 0; i < (const int)nCases; i++) // keep bxOrthCenteredT in sync
            bxOrthCenteredT_(nTerms1,i) = 0;
    }
}

//-----------------------------------------------------------------------------
static int GetNbrUsed(   // Nm in Friedman's notation
    const size_t nCases, // in
    const int iParent,   // in
    const double bx[])   // in: MARS basis matrix
{
    int nUsed = 0;
    if(bx == NULL)
        nUsed = (const int)nCases;
    else for(int i = 0; i < (const int)nCases; i++)
        if(bx_(i,iParent) > 0)
            nUsed++;
    return nUsed;
}

//-----------------------------------------------------------------------------
static int GetEndSpan(
    const int    nPreds,
    const int    nDegree, // in: degree of current term (for adjusting endspan)
    const size_t nCases)
{
    int nEndSpan = 1;
    if(nEndSpanGlobal > 0)                      // user specified endspan?
        nEndSpan = nEndSpanGlobal;
    else if(nEndSpanGlobal == 0) {              // auto?
        // eqn 45 FriedmanMars (see refs)
        static const double log_2 = 0.69315;    // log(2)
        static const double temp1 = 7.32193;    // 3 + log(20)/log(2);
        nEndSpan = (int)(temp1 + log((double)nPreds) / log_2);
    } else                                      // negative endspan illegal
        error("endspan %d < 0", nEndSpanGlobal);
    if(nDegree >= 2) // .5 below makes (int) cast act as round
        nEndSpan += (int)(AdjustEndSpanGlobal * nEndSpan + .5);
    if(nEndSpan > (const int)nCases / 2 - 1)  // always at least one knot, so above adjustment
        nEndSpan = (const int)nCases / 2 - 1; // doesn't completely inhibit degree2 terms
    nEndSpan = max(1, nEndSpan);
    return nEndSpan;
}

//-----------------------------------------------------------------------------
static void GetSpanParams(
    int* pnMinSpan,       // out: number cases between knots
    int* pnEndSpan,       // out: number of cases ignored on each end
    int* pnStartSpan,     // out: number of cases from end until first knot
    const size_t nCases,  // in
    const int    nPreds,  // in
    const int    nDegree, // in: degree of current term (for adjusting endspan)
    const int    iParent, // in
    const double bx[])    // in: MARS basis matrix, can be NULL
{
    const int nEndSpan = GetEndSpan(nPreds, nDegree, nCases);

    int nStartSpan = 0, nMinSpan = 0;
    if(nMinSpanGlobal < 0) {           // treat negative minspan as number of knots
        // get nMinSpan
        nMinSpan = (int)(ceil(nCases / (1.-nMinSpanGlobal))); // convert nknots to minspan
        // get nStartSpan
        nStartSpan = nMinSpan;
        while(nStartSpan < nEndSpan)
            nStartSpan += nMinSpan;
        nStartSpan--;
        nStartSpan = max(1, nStartSpan);
    } else {
        // get nMinSpan
        if(nMinSpanGlobal > 0)          // user specified minspan?
            nMinSpan = nMinSpanGlobal;
        else if(nMinSpanGlobal == 0) {  // auto?
            // eqn 43 in FriedmanMars paper (see refs)
            const int nUsed = GetNbrUsed(nCases, iParent, bx); // Nm in Friedmans notation
            static const double temp1 = 2.9702;     // -log(-log(0.95)
            static const double temp2 = 1.7329;     // 2.5 * log(2)
            nMinSpan = (int)((temp1 + log((double)(nPreds * nUsed))) / temp2);
        }
        nMinSpan = max(1, nMinSpan);
        // get nStartSpan
        const int nAvail = max(0, (const int)nCases - 2 * nEndSpan);
        nStartSpan = nAvail / 2;        // if space for only one knot, put it in center
        if(nAvail > nMinSpan) {         // space for more than one knot?
            const int nDiv = nAvail / nMinSpan;
            if(nAvail == nDiv * nMinSpan)
                nStartSpan = nMinSpan / 2;
            else
                nStartSpan = (nAvail -  nDiv * nMinSpan) / 2;
        }
        // TODO consider moving the following line of code into the "}" above
        nStartSpan = max(1, nEndSpan + nStartSpan);
    }
    *pnStartSpan = nStartSpan;
    *pnMinSpan = nMinSpan;
    *pnEndSpan = nEndSpan;
}

//-----------------------------------------------------------------------------
// The caller has selected a candidate predictor iPred and a candidate iParent.
// This function now selects a knot.  If it finds a knot it will
// update *piBestCase and pRssDeltaForParPredPair.
//
// The general idea: scan backwards through all (ordered) values (i.e. potential
// knots) for the given predictor iPred, calculating RssDelta.
// If RssDelta > *pRssDeltaForParPredPair (and all else is ok), then
// select the knot (by updating *piBestCase and *pRssDeltaForParPredPair).
//
// We want to add a term pair at index iNewCol and iNewCol+1.
// There are currently nTerms in the model.
//
// This function must be fast.

static INLINE void FindKnot(
    int*    piBestCase,              // out: possibly updated, row index in x
    double* pRssDeltaForParPredPair, // io: updated if knot is better
    double  CovCol[],                // scratch buffer, overwritten, nTerms x 1
    double  ycboSum[],               // scratch buffer, overwritten, nMaxTerms x nResp
    double  CovSx[],                 // scratch buffer, overwritten, nTerms x 1
    double* ybxSum,                  // scratch buffer, overwritten, nResp x 1
    const int iNewCol,               // in: tentative knot goes into bx[iNewCol]
    const int iParent,               // in: parent term
    const int iPred,                 // in: predictor index
    const size_t nCases,             // in
    const int nResp,                 // in: number of cols in y
    const int nMaxTerms,             // in
    const double RssDeltaLin,        // in: change in RSS if predictor iPred enters linearly
    const double MaxLegalRssDelta,   // in: FindKnot rejects any changes in Rss greater than this
    const double bx[],               // in: MARS basis matrix
    const double bxOrth[],           // in
    const double bxOrthCenteredT[],  // in
    const double bxOrthMean[],       // in
    const double x[],                // in: nCases x nPreds
    const double y[],                // in: nCases x nResp
    const int xOrder[],              // in
    const double yMean[],            // in: vector nResp x 1
    const int nStartSpan,            // in: number of cases from end until first knot
    const int nMinSpan,              // in: number cases between knots
    const int nEndSpan,              // in: number of cases ignored on each end
    const double NewVarAdjust,       // in: 1 if not a new var, 1/(1+NewVarPenalty) if new var
    const double RssBeforeNewTerm)   // in: used only when trace >= 8
{
    tprintf(8, "--FindKnotBegin-- iPred %d iNewCol %d RssBeforeAddingHinge %g "
        "nMinSpan %d nEndSpan %d nStartSpan %d\n",
        iPred+IOFFSET, iNewCol+IOFFSET, RssBeforeNewTerm,
        nMinSpan, nEndSpan, nStartSpan);

    ASSERT(MaxLegalRssDelta > 0);

    // Tol was .01 prior to earth 4.4.0 but that caused zigzag runout in test.weights.R.
    // The comparison against iNewCol provides some measure of back compatibility and
    // helps prevent overfitting in small models.  The 15 is fairly arb but was chosen
    // to keep small models from having nearby knots, for the zizgag function in
    // test.weights.R, and for the spkmap neural data.
    // TODO This isn't a clean solution.
    const double Tol = iNewCol < 15? .01 : 1e-5;

    int iResp;
    for(iResp = 0; iResp < nResp; iResp++)
        ycboSum_(iNewCol, iResp) = 0;
    memset(CovCol, 0, (iNewCol+1) * sizeof(double));
    memset(CovSx,  0, (iNewCol+1) * sizeof(double));
    memset(ybxSum, 0, nResp * sizeof(double));
    double bxSum = 0, bxSqSum = 0, bxSqxSum = 0, bxxSum = 0, st = 0;
    int iSpan = nStartSpan;
    for(int i = (const int)nCases-2; i >= nEndSpan; i--) { // -2 allows for ix1
        // may Mars have mercy on the poor soul who enters here
        const int    ix0 = xOrder_(i, iPred);   // get the x's in descending order
        const double x0  = x_(ix0,iPred);       // the knot (printed as Cut in trace prints)
        const int    ix1 = xOrder_(i+1, iPred);
        const double x1  = x_(ix1, iPred);      // case next to the cut
        const double bx1 = bx_(ix1, iParent);
        const double bxSq = sq(bx1);
        const double xDelta = x1 - x0;          // will always be non negative
        if(USE_BLAS) {
            daxpy_(&iNewCol, &bx1, &bxOrthCenteredT_(0,ix1), &ONE, CovSx,  &ONE);
            daxpy_(&iNewCol, &xDelta, CovSx, &ONE, CovCol, &ONE);
        } else for(int it = 0; it < iNewCol; it++) {
            CovSx[it]  += (bxOrth_(ix1,it) - bxOrthMean[it]) * bx1;
            CovCol[it] += xDelta * CovSx[it];
        }
        bxSum    += bx1;
        bxSqSum  += bxSq;
        bxxSum   += bx1 * x1;
        bxSqxSum += bxSq * x1;
        const double su = st;
        st = bxxSum - bxSum * x0;

        CovCol[iNewCol] += xDelta * (2 * bxSqxSum - bxSqSum * (x0 + x1)) +
                           (sq(su) - sq(st)) / nCases;

        if(nResp == 1) {    // treat nResp==1 as a special case, for speed
            ybxSum[0] += (y_(ix1, 0) - yMean[0]) * bx1;
            ycboSum_(iNewCol, 0) += xDelta * ybxSum[0];
        } else for(iResp = 0; iResp < nResp; iResp++) {
            ybxSum[iResp] += (y_(ix1, iResp) - yMean[iResp]) * bx1;
            ycboSum_(iNewCol, iResp) += xDelta * ybxSum[iResp];
        }
        if(bx1 > 0 && CovCol[iNewCol] > 0 && --iSpan == 0) {
            iSpan = nMinSpan;
            double RssDelta = 0;
            bool Best = false, TolGood = true;
            // calculate RssDelta and see if this knot beats the previous best
            RssDelta = 0;
            double temp1, temp2;
            for(iResp = 0; iResp < nResp; iResp++) {
                if(USE_BLAS) {
                    temp1 =
                        ycboSum_(iNewCol,iResp) -
                        ddot_(&iNewCol, &ycboSum_(0,iResp), &ONE, CovCol, &ONE);
                    temp2 =
                        CovCol[iNewCol] - ddot_(&iNewCol, CovCol, &ONE, CovCol, &ONE);
                } else {
                    temp1 = ycboSum_(iNewCol,iResp);
                    temp2 = CovCol[iNewCol];
                    for(int it = 0; it < iNewCol; it++) {
                        temp1 -= ycboSum_(it,iResp) * CovCol[it];
                        temp2 -= sq(CovCol[it]);
                    }
                }
                // TODO HastieTibs code has a comment saying the following has to be fixed
                if(temp2 / CovCol[iNewCol] > Tol)
                    RssDelta += sq(temp1) / temp2;
                else
                    TolGood = false;
            }
            RssDelta = NewVarAdjust * (RssDeltaLin + RssDelta);
            if(RssDelta > *pRssDeltaForParPredPair &&
                    RssDelta < MaxLegalRssDelta) {
                *piBestCase = i;
                *pRssDeltaForParPredPair = RssDelta;
                Best = true;
            }
            if(TraceGlobal >= 8) {
                const double RssWithKnot = RssBeforeNewTerm - RssDelta;
                printf("--FindKnot--Case %4d RssWithKnot %12.5g RssDelta %12.5g "
                    "Cut % 12.5g ",
                    i+IOFFSET, MaybeZero(RssWithKnot), MaybeZero(RssDelta), x0);
                // separate trace to make it easier to check compat
                // with FindWeightedKnot when trace==8
                tprintf(9, "bx1G %d CovColG %d TolG %d MaxG %d ",
                    bx1 > 0, CovCol[iNewCol] > 0, TolGood,
                    RssDelta < MaxLegalRssDelta);
                printf("%s\n", Best? " best": "");
            }
        } else if(TraceGlobal >= 9)
            printf("--FindKnot--Case %4d iSpan %d bx1 % 8.4f\n",
                i+IOFFSET, iSpan, bx1);
    } // for
    tprintf(8, "--FindKnotEnd--\n");
}

//-----------------------------------------------------------------------------
// Add a candidate term at bx[,nTerms], with the parent term multiplied by
// the predictor iPred entering linearly.  Do this by setting the knot at
// the lowest value xMin of x, since max(0,x-xMin)==x-xMin for all x.  The
// change in RSS caused by adding this term forms the base RSS delta which
// we will try to beat in the search in FindKnot.
//
// This also initializes CovCol, bxOrth[,nTerms], and ycboSum[nTerms,]

static INLINE void AddCandidateLinearTerm(
    double* pRssDeltaLin,     // out: change to RSS caused by adding new term
    bool*   pIsNewForm,       // out: true on entry, may be cleared by InitBxOrthCol
    double  xbx[],            // out: nCases x 1
    double  CovCol[],         // out: nMaxTerms x 1
    double ycboSum[],         // io: nMaxTerms x nResp
    double bxOrth[],          // io
    double bxOrthCenteredT[], // io
    double bxOrthMean[],      // io
    const int iPred,          // in
    const int iParent,        // in
    const double x[],         // in: nCases x nPreds
    const double y[],         // in: nCases x nResp, scaled y
    const size_t nCases,      // in
    const int nResp,          // in: number of cols in y
    const int nTerms,         // in
    const int nMaxTerms,      // in
    const double bx[],        // in: MARS basis matrix
    const bool FullSet[])     // in
{
    // set xbx to x[,iPred] * bx[,iParent]
    // note: when iParent==1, bx_[,iParent] is all ones, therefore xbx is x

    int i;
    for(i = 0; i < (const int)nCases; i++)
        xbx[i] = x_(i,iPred) * bx_(i,iParent);

    // Init bxOrth[,nTerms] and bxOrthMean[nTerms] for the candidate term.
    // Clears both those columns, bxOrthMean, and *pIsNewForm if
    // column sum-of-squares is less than MIN_BX_SOS.

    InitBxOrthCol(bxOrth, bxOrthCenteredT, bxOrthMean, pIsNewForm,
        xbx, nTerms, FullSet, nCases, nMaxTerms, iParent, iPred);

    // init CovCol and ycboSum[nTerms], for use by FindKnot later

    memset(CovCol, 0, (nTerms-1) * sizeof(double));
    CovCol[nTerms] = 1;
    int iResp;
    for(iResp = 0; iResp < nResp; iResp++) {
        ycboSum_(nTerms, iResp) = 0;
        for(i = 0; i < (const int)nCases; i++)
            ycboSum_(nTerms, iResp) += (y_(i, iResp) - yMean[iResp]) *
                                       bxOrth_(i,nTerms);
    }
    // calculate change to RSS caused by adding candidate new term

    *pRssDeltaLin = 0;
    for(iResp = 0; iResp < nResp; iResp++) {
        double yboSum = 0;
        for(i = 0; i < (const int)nCases; i++)
            yboSum += y_(i,iResp) * bxOrth_(i,nTerms);
        *pRssDeltaLin += sq(yboSum);
    }
}

//-----------------------------------------------------------------------------
// The caller has selected a candidate parent term iParent.
// This function now selects a predictor, and a knot for that predictor.
//
// TODO These functions have a ridiculous number of parameters, I know.
//
// TODO A note on the comparison against ALMOST_ZERO below:
// It's not a clean solution but seems to work ok.
// It was added after we saw different results on different
// machines for certain datasets e.g. (tested on earth 1.4.0)
// ldose  <- rep(0:5, 2) - 2
// ldose1 <- c(0.1, 1.2, 2.3, 3.4, 4.5, 5.6, 0.3, 1.4, 2.5, 3.6, 4.7, 5.8)
// sex3 <- factor(rep(c("male", "female", "andro"), times=c(6,4,2)))
// fac3 <- factor(c("lev2", "lev2", "lev1", "lev1", "lev3", "lev3",
//                  "lev2", "lev2", "lev1", "lev1", "lev3", "lev3"))
// numdead <- c(1,4,9,13,18,20,0,2,6,10,12,16)
// numdead2 <- c(2,3,10,13,19,20,0,3,7,11,13,17)
// pair <- cbind(numdead, numdead2)
// df <- data.frame(sex3, ldose, ldose1, fac3)
// am <-  earth(df, pair, trace=6, pmethod="none", degree=2)

static INLINE void FindPredGivenParent(
    int*    piBestCase,             // out: untouched unless an improving term is found
    int*    piBestPred,             // out: ditto
    int*    piBestParent,           // out: existing term on which we are basing the new term
    double* pBestRssDeltaForTerm,   // io: untouched unless an improving term is found
    double* pBestRssDeltaForParent, // io: used only by FAST_MARS
    bool*   pIsNewForm,             // out
    bool*   pLinPredIsBest,         // out: true if pred should enter linearly (no knot)

    double bxOrth[],                // io
    double bxOrthCenteredT[],       // io
    double bxOrthMean[],            // io
    double xbx[],                   // io: nCases x 1
    double CovSx[],                 // io: nMaxTerms x 1
    double CovCol[],                // io: nMaxTerms x 1
    double ycboSum[],               // io: nMaxTerms x nResp

    const double bx[],              // in: MARS basis matrix
    const double yMean[],           // in: vector nResp x 1
    const double RssBeforeNewTerm,  // in
    const double MaxLegalRssDelta,  // in: FindKnot rejects any changes in Rss greater than this

    const int iParent,              // in
    const double x[],               // in: nCases x nPreds, unweighted x
    const double y[],               // in: nCases x nResp, unweighted but scaled y
    const size_t nCases,            // in
    const int nResp,                // in: number of cols in y
    const int nPreds,               // in
    const int nTerms,               // in
    const int nMaxTerms,            // in

    const bool FullSet[],           // in
    const int xOrder[],             // in: order of each column of x array
    const int nUses[],              // in: nbr of times each pred is used in the model
    const int Dirs[],               // in
    const double NewVarPenalty,     // in: penalty for adding a new variable (default is 0)
    const int LinPreds[],           // in: nPreds x 1, 1 if predictor must enter linearly

    const int nMinSpan,             // in
    const int nEndSpan,             // in
    const int nStartSpan)           // in
{
    double* ybxSum = (double*)malloc1(nResp * sizeof(double),  // working var for FindKnot
                        "ybxSum\t\tnResp %d sizeof(double) %d",
                        nResp, sizeof(double));

    bool UpdatedBestRssDelta = false;
    for(int iPred = 0; iPred < nPreds; iPred++) {
        tprintf(9, "\n");
        if(Dirs_(iParent,iPred) != 0) {     // predictor is in parent term?
            tprintf(7,
"|Parent %-2d Pred %-2d %44.44s skip (pred is in parent)\n",
                    iParent+IOFFSET, iPred+IOFFSET, " ");
#if USING_R
        } else if(!IsAllowed(iPred, iParent, Dirs, nPreds, nMaxTerms)) {
            tprintf(7,
"|Parent %-2d Pred %-2d %44.44s skip (not allowed by \"allowed\" func)\n",
                    iParent+IOFFSET, iPred+IOFFSET, " ");
#endif
        } else {
            // we apply the penalty if the variable is entering for the first time
            const double NewVarAdjust = 1 / (1 + (nUses[iPred] == 0? NewVarPenalty: 0));
            double RssDeltaLin = 0; // change in RSS for iPred entering linearly
            double UnadjustedRssDeltaLin = 0;
            double RssBeforeKnot = RssBeforeNewTerm;
            bool IsNewForm = GetNewFormFlag(iPred, iParent, Dirs,
                                            FullSet, nTerms, nPreds, nMaxTerms);
            if(IsNewForm) {
                // Create a candidate term at bx[,nTerms],
                // with iParent and iPred entering linearly
                // This may clear IsNewForm.
                AddCandidateLinearTerm(&UnadjustedRssDeltaLin, &IsNewForm,
                    xbx, CovCol, ycboSum, bxOrth, bxOrthCenteredT, bxOrthMean,
                    iPred, iParent, x, y,
                    nCases, nResp, nTerms, nMaxTerms, bx, FullSet);
                RssDeltaLin = NewVarAdjust * UnadjustedRssDeltaLin;
                tprintf(8,
"\n|Parent %-2d Pred %-2d Case   -1 Cut % 12.4g< Rss %-12.5g",
                        iParent+IOFFSET, iPred+IOFFSET,
                        GetCut(0, iPred, nCases, x, xOrder),
                        MaybeZero(RssBeforeNewTerm - RssDeltaLin));
                tprintf(9, " RssDeltaLin %-12.5g ", RssDeltaLin);
                if(fabs(RssDeltaLin - *pBestRssDeltaForTerm) < ALMOST_ZERO) {
                    RssDeltaLin = *pBestRssDeltaForTerm; // see header note
                    tprintf(7, "RssDelta %g is ALMOST_ZERO\n",
                            RssDeltaLin - *pBestRssDeltaForTerm);
                }
                if(RssDeltaLin > *pBestRssDeltaForParent)
                    *pBestRssDeltaForParent = RssDeltaLin;
                RssBeforeKnot -= RssDeltaLin;
                if(RssDeltaLin > *pBestRssDeltaForTerm) {
                    // The new term (with predictor entering linearly) beats other
                    // candidate terms so far.
                    tprintf(9, "best for term (lin pred) ");
                    UpdatedBestRssDelta = true;
                    *pBestRssDeltaForTerm = RssDeltaLin;
                    *pLinPredIsBest = true;
                    *piBestCase     = 0;        // knot is at the lowest value of x
                    *piBestPred     = iPred;
                    *piBestParent   = iParent;
                }
                tprintf(8, "\n");
            } else
                tprintf(8, "|Parent %-2d Pred %-2d no new form\n",
                        iParent+IOFFSET, iPred+IOFFSET);
            double RssDeltaForParPredPair = RssDeltaLin;
            if(!LinPreds[iPred]) {
                int iBestCase = -1;
                FindKnot(&iBestCase, &RssDeltaForParPredPair,
                        CovCol, ycboSum, CovSx, ybxSum,
                        (IsNewForm? nTerms + 1: nTerms),
                        iParent, iPred, nCases, nResp, nMaxTerms,
                        UnadjustedRssDeltaLin, MaxLegalRssDelta,
                        bx, bxOrth, bxOrthCenteredT, bxOrthMean,
                        x, y, xOrder, yMean,
                        nStartSpan, nMinSpan, nEndSpan, NewVarAdjust,
                        RssBeforeNewTerm); // RssBeforeNewTerm is for tracing
                if(RssDeltaForParPredPair > *pBestRssDeltaForParent)
                    *pBestRssDeltaForParent = RssDeltaForParPredPair;
                if(RssDeltaForParPredPair > *pBestRssDeltaForTerm) {
                    UpdatedBestRssDelta = true;
                    *pBestRssDeltaForTerm = RssDeltaForParPredPair;
                    *pLinPredIsBest = false;
                    *piBestCase     = iBestCase;
                    *piBestPred     = iPred;
                    *piBestParent   = iParent;
                    *pIsNewForm     = IsNewForm;
                    tprintf(7,
"|Parent %-2d Pred %-2d Case %4d Cut % 12.4g  Rss %-12.5g RssDelta %-12.5g%s\n",
                        iParent+IOFFSET, iPred+IOFFSET,
                        iBestCase+IOFFSET,
                        GetCut(iBestCase, iPred, nCases, x, xOrder),
                        MaybeZero(RssBeforeNewTerm - *pBestRssDeltaForTerm),
                        MaybeZero(*pBestRssDeltaForTerm),
                        RssDeltaForParPredPair > *pBestRssDeltaForTerm?
                            " best for term": "");
                } else
                    tprintf(7,
"|Parent %-2d Pred %-2d Case %4d Cut % 12.4g< Rss %-12.5g RssDelta %-12.5g\n",
                        iParent+IOFFSET, iPred+IOFFSET, -1,
                        GetCut(0, iPred, nCases, x, xOrder),
                        MaybeZero(RssBeforeNewTerm - RssDeltaLin),
                        0.);
            }
        } // else
    } // for iPred
    free1(ybxSum);
    if(UpdatedBestRssDelta && NewVarPenalty != 0. && nUses[*piBestPred] == 0) {
        // we applied NewVarPenalty earlier, now un-apply it
        *pBestRssDeltaForTerm *= 1 + NewVarPenalty;
    }
}

//-----------------------------------------------------------------------------
#if WEIGHTS
// TODO This could be made faster by caching previous results?

static INLINE void InitHinge(   // return TRUE if column is not all zeros
    double bxCol[],             // out: this column will be initialized
    const int iHinge,           // in: hinge index (row in x)
    const double bxParentCol[], // in: column in bx for iParent
    const double xCol[],        // in: column of x for iPred
    const int xOrderCol[],      // in: column in xOrder for iPred
    const size_t nCases)        // in
{
    const double Cut = xCol[xOrderCol[iHinge]];
    for(int i = (const int)nCases-1; i > iHinge; i--) {
        const int ix = xOrderCol[i];
        bxCol[ix] = bxParentCol[ix] * (xCol[ix] - Cut);
    }
}
#endif

//-----------------------------------------------------------------------------
#if WEIGHTS
static double GetRegressionRss(
    double x[],          // io: nCases x nCols, gets overwritten
    const double y[],    // in: nCases x nResp
    const size_t nCases, // in: number of rows in x and in y
    const int nResp,     // in: number of cols in y
    int nCols,           // in: number of columns in x
    double Residuals[],  // in: nCases, working storage
    int    iPivots[],    // in: nCols, working storage
    double qraux[],      // in: nCols, working storage
    double work[])       // in: max(nCols * 2, nCases), working storage
{
    for(int iCol = 0; iCol < nCols; iCol++)
        iPivots[iCol] = iCol+1;

    int nCases1 = (const int)nCases; // type convert from size_t

    int nRank;
    dqrdc2_(                // R function, QR decomp based on LINPACK dqrdc
        x,                  // io:  x, on return upper tri of x is R of QR
        &nCases1,           // in:  ldx
        &nCases1,           // in:  n
        &nCols,             // in:  p
        (double*)&QR_TOL,   // in:  tol
        &nRank,             // out: k, num of indep cols of x
        qraux,              // out: qraux
        iPivots,            // out: jpvt
        work);              // work

    double Rss = 0;

    int job = 10; // specify 10 because all we need are the residuals
    int info;
    for(int iResp = 0; iResp < nResp; iResp++) {
        dqrsl_(                 // LINPACK function
            x,                  // in:  x, generated by dqrdc2
            &nCases1,           // in:  ldx
            &nCases1,           // in:  n
            &nRank,             // in:  k
            qraux,              // in:  qraux
            (double*)(y + iResp * nCases),  // in:  y
            NULL,               // out: qy, unreferenced here
            work,               // out: qty, unused here but needed for dqrsl
            NULL,               // out: b, unreferenced here
            (double*)Residuals, // out: rsd
            NULL,               // out: xb, unreferenced here
            &job,               // in:  job
            &info);             // in:  info

        ASSERT(info == 0);

        for(int i = 0; i < (const int)nCases; i++)
            Rss += sq(Residuals[i]);
    }
    return Rss;
}
#endif // WEIGHTS

//-----------------------------------------------------------------------------
#if WEIGHTS
static INLINE void FindWeightedKnot(
    int*    piBestCase,     // out: updated, row index in x
    double* pRssBestKnot,   // io: updated, on entry is LinRss if did lin pred
    const bool UsedCols[],  // in
    int iNewCol,            // in: tentative knot goes into bx[iNewCol]
    const int iParent,      // in: parent term
    const int iPred,        // in: predictor index
    const size_t nCases,    // in
    const int nResp,        // in: number of cols in yw
    double bx[],            // in: MARS basis matrix, columns nTerm and nTerm+1 filled in
    const double x[],       // in: nCases x nPreds, unweighted x matrix
    const double yw[],      // in: nCases x nResp, weighted and scaled y matrix
    const int xOrder[],     // in: order of each column of _unweighted_ x array
    const int nStartSpan,   // in: number of cases from end until first knot
    const int nMinSpan,     // in: number cases between knots
    const int nEndSpan,     // in: number of cases ignored on each end
    const double RssBeforeNewTerm) // in: used only for tracing
{
    tprintf(8, "--FindKnotBegin-- iPred %d iNewCol %d RssBeforeAddingHinge %g "
        "nMinSpan %d nEndSpan %d nStartSpan %d\n",
        iPred+IOFFSET, iNewCol+IOFFSET, RssBeforeNewTerm,
        nMinSpan, nEndSpan, nStartSpan);

    // we malloc everything for GetRegressionRss once here instead of in the loop

    double* bxUsed;
    const int nUsedCols = CopyUsedCols(&bxUsed, bx, nCases, iNewCol+1, UsedCols);

    double* bxTemp = (double*)malloc1(nCases * nUsedCols * sizeof(double),
                            "bxTemp\t\tnCases %d nUsedCols %d sizeof(double) %d",
                            (const int)nCases, nUsedCols, sizeof(double));

    double* Residuals = (double*)malloc1(nCases * sizeof(double),
                            "Residuals\t\tnCases %d sizeof(double) %d",
                            (const int)nCases, sizeof(double));

    int* iPivots = (int*)malloc1(nUsedCols * sizeof(int),
                            "iPivots\t\tnUsedCols %d sizeof(int) %d",
                            nUsedCols, sizeof(int));

    double* qraux = (double*)malloc1(nUsedCols * sizeof(double),
                            "qraux\t\t\tnUsedCols %d sizeof(double) %d",
                            nUsedCols, sizeof(double));

    // in GetRegressionRss, work must be p*2 for dqrdc2, and
    // nCases in dqrsl where it is used temp storage for qty

    double* work = (double*)malloc1(
                            max(nUsedCols * 2, (const int)nCases) * sizeof(double),
                            "work\t\t\tnCases %d sizeof(double) %d",
                            (const int)nCases, sizeof(double));

    // zero the current column of bxUsed, we will fill it in with the hinge functions
    memset(&bxUsed_(0,iNewCol), 0, nCases * sizeof(double));
    // for-loop indices and trace prints are compatible with FindKnot (when trace<=8)
    int iSpan = nStartSpan;
    for(int i = (const int)nCases-2; i >= nEndSpan; i--) {
        // TODO This should be bx0 (not bx1)?  But for compat with FindKnot we use bx1.
        //      In test.mods.R, using bx1 or bx0 here give almost identical results.
        const int ix1 = xOrder_(i+1, iPred);
        const double bx1 = bx_(ix1,iParent);
        if(bx1 > 0 && --iSpan == 0) {
            iSpan = nMinSpan;
            bool Best = false;
            InitHinge(&bxUsed_(0, iNewCol), // init this column of bxUsed
                      i, &bx_(0, iParent), &x_(0, iPred), &xOrder_(0, iPred), nCases);
            // bxTemp is needed because GetRegressionRss destroys its first arg
            memcpy(bxTemp, bxUsed, nCases * nUsedCols * sizeof(double));
            const double KnotRss =
                GetRegressionRss(bxTemp, yw, nCases, nResp, nUsedCols,
                                 Residuals, iPivots, qraux, work);
            // using ALMOST_ZERO here gives results closer to FindKnot
            if(KnotRss < *pRssBestKnot - ALMOST_ZERO) {
                *piBestCase = i;
                *pRssBestKnot = KnotRss;
                Best = true;
            }
            if(TraceGlobal >= 8) {
                const double RssDelta = RssBeforeNewTerm - KnotRss;
                printf(
"--FindKnot--Case %4d RssWithKnot %12.5g RssDelta %12.5g Cut % 12.5g%s\n",
                    i+IOFFSET, MaybeZero(KnotRss), MaybeZero(RssDelta),
                    GetCut(i, iPred, nCases, x, xOrder), Best? " best": "");
            }
        } else if(TraceGlobal >= 9)
            printf("--FindKnot--Case %4d iSpan %d bx1 % 8.4f\n", i+IOFFSET, iSpan, bx1);
    }
    free1(work);
    free1(qraux);
    free1(iPivots);
    free1(Residuals);
    free1(bxTemp);
    free1(bxUsed);

    tprintf(8, "--FindKnotEnd--\n");
}
#endif // WEIGHTS

//-----------------------------------------------------------------------------
#if WEIGHTS
static INLINE void FindWeightedPredGivenParent(
    int*    piBestCase,             // out: untouched unless an improving term is found
    int*    piBestPred,             // out: ditto
    int*    piBestParent,           // out: existing term on which we are basing the new term
    double* pBestRssDeltaForTerm,   // io: untouched unless an improving term is found
    double* pBestRssDeltaForParent, // io: used only by FAST_MARS
    bool*   pIsNewForm,             // out
    bool*   pLinPredIsBest,         // out: true if pred should enter linearly (no knot)

    double bx[],                    // in: MARS basis matrix
    const double RssBeforeNewTerm,  // in

    const int iParent,              // in
    const double x[],               // in: nCases x nPreds, unweighted x
    const double yw[],              // in: nCases x nResp, weighted and scaled y
    const size_t nCases,            // in
    const int nResp,                // in: number of cols in y
    const int nPreds,               // in
    const int nTerms,               // in
    const int nMaxTerms,            // in

    const bool FullSet[],           // in
    const int xOrder[],             // in: order of each column of _unweighted_ x array
    const int nUses[],              // in: nbr of times each pred is used in the model
    const int Dirs[],               // in
    const double NewVarPenalty,     // in: penalty for adding a new variable (default is 0)
    const int LinPreds[],           // in: nPreds x 1, 1 if predictor must enter linearly

    const int nMinSpan,             // in
    const int nEndSpan,             // in
    const int nStartSpan)           // in
{
    if(NewVarPenalty != 0) // NewVarPenalty is not yet fully tested when weights are used
        error("newvar.penalty is not yet implemented with weights");
    bool UpdatedBestRssDelta = false;

    bool* UsedCols = (bool*)calloc1(nMaxTerms, sizeof(bool),
                        "UsedCols\t\tnMaxTerms %d sizeof(bool) %d",
                        nMaxTerms, sizeof(bool));
    for(int iTerm = 0; iTerm < nTerms; iTerm++)
        UsedCols[iTerm] = true;
    for(int iPred = 0; iPred < nPreds; iPred++) {
        tprintf(8, "\n");
        if(Dirs_(iParent,iPred) != 0) {     // predictor is in parent term?
            tprintf(7, "|Parent %-2d Pred %-2d %44.44s skip (pred is in parent)\n",
                iParent+IOFFSET, iPred+IOFFSET, " ");
#if USING_R
        } else if(!IsAllowed(iPred, iParent, Dirs, nPreds, nMaxTerms)) {
            tprintf(7,
"|Parent %-2d Pred %-2d %44.44s skip (not allowed by \"allowed\" func)\n",
                iParent+IOFFSET, iPred+IOFFSET, " ");
#endif
        } else {
#if USING_R
            // TODO we don't release UsedCols here if user interrupts
            ServiceR();
#endif
            // const double NewVarAdjust = 1 + (nUses[iPred] == 0? NewVarPenalty: 0);
            bool IsNewForm = GetNewFormFlag(iPred, iParent, Dirs,
                                FullSet, nTerms, nPreds, nMaxTerms);
            ASSERT(nTerms+1 < nMaxTerms);
            UsedCols[nTerms] = UsedCols[nTerms+1] = false;
            double RssBeforeKnot = RssBeforeNewTerm;
            const int iNewCol = IsNewForm? nTerms+1: nTerms;
            if(IsNewForm) {
                // Add the new predictor as a linear term in bx[,nTerms].
                // This updates RssBeforeKnot, which we will
                // try to beat in FindWeightedKnot.
                // Note that unlike FindPredGivenParent we never clear IsNewForm
                // here (that will be done later if necessary in AddTermPair).
                for(int i = 0; i < (const int)nCases; i++)
                    bx_(i, nTerms) = bx_(i, iParent) * x_(i, iPred);
                UsedCols[nTerms] = true;
                Regress(NULL, NULL, &RssBeforeKnot, NULL, NULL, NULL,
                    bx, yw, nCases, nResp, nMaxTerms, UsedCols);
                tprintf(8,
"|Parent %-2d Pred %-2d Case   -1 Cut % 12.4g< Rss %-12.5g\n",
                        iParent+IOFFSET, iPred+IOFFSET,
                        GetCut(0, iPred, nCases, x, xOrder),
                        MaybeZero(RssBeforeKnot));
            } else {
                tprintf(8, "\n|Parent %-2d Pred %-2d no new form\n",
                    iParent+IOFFSET, iPred+IOFFSET);
#if 1 // TODO remove this slow check when weights code has been fully tested
                Regress(NULL, NULL, &RssBeforeKnot, NULL, NULL, NULL,
                    bx, yw, nCases, nResp, nMaxTerms, UsedCols);
                if(fabs(RssBeforeKnot - RssBeforeNewTerm) > RssBeforeKnot * 1e-6)
                    error(
"fabs(RssBeforeKnot %g - RssBeforeNewTerm %g) %g > %g",
                           RssBeforeKnot, RssBeforeNewTerm,
                           fabs(RssBeforeKnot - RssBeforeNewTerm),
                           RssBeforeKnot * 1e-6);
#endif
            }
            double RssBestKnot = RssBeforeKnot;
            int iBestCase = 0;
            if(!LinPreds[iPred]) {
#if 1 // TODO remove this slow check when weights code has been fully tested
                double RssTemp;
                Regress(NULL, NULL, &RssTemp, NULL, NULL, NULL,
                    bx, yw, nCases, nResp, nMaxTerms, UsedCols);
                if(fabs(RssTemp - RssBeforeKnot) > RssBeforeKnot * 1e-6)
                    error("fabs(RssTemp - RssBeforeKnot) > %g",
                          RssTemp, RssBeforeKnot,
                          fabs(RssTemp - RssBeforeKnot),
                          RssBeforeKnot * 1e-6);
#endif
                ASSERT(iNewCol > 0 && iNewCol < nMaxTerms);
                UsedCols[iNewCol] = true;

                FindWeightedKnot(&iBestCase, &RssBestKnot,
                        UsedCols, iNewCol,
                        iParent, iPred, nCases, nResp,
                        bx, x, yw, xOrder, nStartSpan, nMinSpan, nEndSpan,
                        RssBeforeNewTerm);
            }
            // TODO must use NewVarAdjust here?
            const bool LinPredIsBest = RssBeforeKnot <= RssBestKnot;
            const double Rss = (LinPredIsBest? RssBeforeKnot: RssBestKnot);
            double RssDeltaForTerm = RssBeforeNewTerm - Rss;
            if(RssDeltaForTerm > *pBestRssDeltaForParent)
                *pBestRssDeltaForParent = RssDeltaForTerm;
            if(RssDeltaForTerm > *pBestRssDeltaForTerm) {
                UpdatedBestRssDelta   = true;
                *pBestRssDeltaForTerm = RssDeltaForTerm;
                *pLinPredIsBest       = LinPredIsBest;
                *piBestCase           = iBestCase;
                *piBestPred           = iPred;
                *piBestParent         = iParent;
                *pIsNewForm           = IsNewForm;
                tprintf(7,
"|Parent %-2d Pred %-2d Case %4d Cut % 12.4g  Rss %-12.5g RssDelta %-12.5g%s\n",
                    iParent+IOFFSET,
                    iPred+IOFFSET,
                    iBestCase+IOFFSET,
                    GetCut(iBestCase, iPred, nCases, x, xOrder),
                    MaybeZero(Rss),
                    MaybeZero(*pBestRssDeltaForTerm),
                    RssDeltaForTerm > *pBestRssDeltaForTerm?
                        " best for term": "");
            } else
                tprintf(7,
"|Parent %-2d Pred %-2d Case %4d Cut % 12.4g  Rss %-12.5g RssDelta %-12.5g\n",
                    iParent+IOFFSET,
                    iPred+IOFFSET,
                    -1,
                    GetCut(0, iPred, nCases, x, xOrder),
                    MaybeZero(Rss),
                    MaybeZero(RssBeforeKnot));
        }
    } // for iPred
    if(UpdatedBestRssDelta && nUses[*piBestPred] == 0) {
        // de-adjust for NewVarPenalty (only makes a difference if NewVarPenalty != 0)
        const double NewVarAdjust = 1 + NewVarPenalty;
        *pBestRssDeltaForTerm *= NewVarAdjust;
    }
    free1(UsedCols);
}
#endif // WEIGHTS

//-----------------------------------------------------------------------------
// Find a new term to add to the model, if possible, and return the
// selected case (i.e. knot), predictor, and parent term indices.
//
// The new term is a copy of an existing parent term but extended
// by multiplying the parent by a new hinge function at the selected knot.
//
// Actually, this usually finds a term _pair_, with left and right hinges.
//
// There are currently nTerms in the model.  We want to add a term at index nTerms.

static void FindTerm(
    int*    piBestCase,           // out: return -1 if no new term available, else row index
    int*    piBestPred,           // out:
    int*    piBestParent,         // out

    double* pBestRssDeltaForTerm, // out: adding new term reduces RSS this much
                                  //      will be set to 0 if no possible new term
    bool*   pIsNewForm,           // out
    bool*   pLinPredIsBest,       // out: true if pred should enter linearly (no knot)

    double bxOrth[],              // io: column nTerms overwritten
    double bxOrthCenteredT[],     // io: kept in sync with bxOrth
    double bxOrthMean[],          // io: element nTerms overwritten

#if WEIGHTS
    double bx[],                  // io: cols at nTerms and nTerms+1 used as scratch, will be set to 0
#else
    const double bx[],            // in: MARS basis matrix
#endif
    const double x[],             // in: nCases x nPreds, unweighted x
    const double y[],             // in: nCases x nResp, unweighted but scaled y
    const double yw[],            // in: nCases x nResp, weighted and scaled y, can be NULL
    const size_t nCases,          // in:
    const int nResp,              // in: number of cols in y
    const int nPreds,             // in:
    const int nMaxDegree,         // in:
    const int nTerms,             // in:
    const int nMaxTerms,          // in:

    const double yMean[],         // in: vector nResp x 1
    const double RssBeforeNewTerm, // in
    const double MaxLegalRssDelta, // in: FindKnot rejects any changes in Rss greater than this
    const bool FullSet[],         // in:
    const int xOrder[],           // in:
    const int nDegree[],          // in: degree of each term, degree of intercept is 0

    const int nUses[],            // in: nbr of times each predictor is used in the model
    const int Dirs[],             // in:
    const int nFastK,             // in: Fast MARS K
    const double NewVarPenalty,   // in: penalty for adding a new variable (default is 0)
    const int LinPreds[])         // in: nPreds x 1, 1 if predictor must enter linearly
{
#if !FAST_MARS // prevent compiler warning: unused parameter
    int Dummy = nFastK;
    ASSERT(Dummy != -999);
#endif
#if !WEIGHTS // prevent compiler warning: unused parameter
    double Dummy1 = RssBeforeNewTerm;
    ASSERT(Dummy1 != -999);
#endif
    tprintf(7, "\n----------------------------------------------------------"
       "-------------------");
    tprintf(7, "\n|FindTerm: Searching for new term %-3d  "
            "RssDelta 0 MaxLegalRssDelta %g\n",
            nTerms+IOFFSET, MaxLegalRssDelta);

    *piBestCase = -1;
    *pBestRssDeltaForTerm = 0;
    *pLinPredIsBest = false;
    *pIsNewForm = false;
    int i;

    xbx = (double*)malloc1(nCases * sizeof(double),
                "xbx\t\t\tnCases %d sizeof(double) %d",
                (const int)nCases, sizeof(double));

    CovSx  = (double*)malloc1(nMaxTerms * sizeof(double),
                "CovSx\t\t\tnMaxTerms %d sizeof(double) %d",
                nMaxTerms, sizeof(double));

    CovCol = (double*)calloc1(nMaxTerms, sizeof(double),
                "CovCol\t\tnMaxTerms %d sizeof(double) %d",
                nMaxTerms, sizeof(double));

    ycboSum  = (double*)calloc1(nMaxTerms * nResp, sizeof(double),
                "ycboSum\t\tnMaxTerms %d nResp %d sizeof(double) %d",
                nMaxTerms, nResp, sizeof(double));

    for(int iResp = 0; iResp < nResp; iResp++)
        for(int iTerm = 0; iTerm < nTerms; iTerm++)
            for(i = 0; i < (const int)nCases; i++)
                ycboSum_(iTerm,iResp) +=
                    (y_(i,iResp) - yMean[iResp]) * bxOrth_(i,iTerm);

#if USING_R
    const int nServiceR = (int)1e6 / nCases;
#endif
    int iParent;
#if FAST_MARS
    GetNextParent(true, nFastK); // init queue iterator
    while((iParent = GetNextParent(false, nFastK)) > -1) {
#else
    for(iParent = 0; iParent < nTerms; iParent++) {
#endif
#if USING_R
        static int iServiceR = 0;
        if(++iServiceR > nServiceR) {
            ServiceR();
            iServiceR = 0;
        }
#endif
        // Assume a bad RssDelta for iParent.  This pushes parent terms that
        // can't be used to the bottom of the FastMARS queue.  (A parent can't
        // be used if its degree is too big or all predictors are in the parent.)

        double BestRssDeltaForParent = -1;    // used only by FAST_MARS

        if(nDegree[iParent] >= nMaxDegree)
            tprintf(7, "|Parent %-2d %52.52s skip (degree of term would be %d)\n",
                iParent+IOFFSET, " ", nDegree[iParent]+1);
        else {
            int nMinSpan, nEndSpan, nStartSpan;
            GetSpanParams(&nMinSpan, &nEndSpan, &nStartSpan,
                          nCases, nPreds, nDegree[iParent]+1, iParent, bx);
#if WEIGHTS
            if(yw)
                FindWeightedPredGivenParent(
                    piBestCase, piBestPred, piBestParent, pBestRssDeltaForTerm,
                    &BestRssDeltaForParent, pIsNewForm, pLinPredIsBest,
                    bx, RssBeforeNewTerm,
                    iParent, x, yw, nCases, nResp, nPreds, nTerms, nMaxTerms,
                    FullSet, xOrder, nUses, Dirs, NewVarPenalty, LinPreds,
                    nMinSpan, nEndSpan, nStartSpan);
            else
#endif
                FindPredGivenParent(
                    piBestCase, piBestPred, piBestParent, pBestRssDeltaForTerm,
                    &BestRssDeltaForParent, pIsNewForm, pLinPredIsBest,
                    bxOrth, bxOrthCenteredT, bxOrthMean, xbx, CovSx, CovCol, ycboSum,
                    bx, yMean, RssBeforeNewTerm, MaxLegalRssDelta,
                    iParent, x, y, nCases, nResp, nPreds, nTerms, nMaxTerms,
                    FullSet, xOrder, nUses, Dirs, NewVarPenalty, LinPreds,
                    nMinSpan, nEndSpan, nStartSpan);
#if FAST_MARS
            UpdateRssDeltaInQ(iParent, nTerms, BestRssDeltaForParent);
#endif
        }
    } // iParent
    tprintf(7, "\n");
    // free in opposite order to alloc to help operating system memory manager
    free1(ycboSum);
    free1(CovCol);
    free1(CovSx);
    free1(xbx);
}

//-----------------------------------------------------------------------------
static void PrintForwardProlog(
    const size_t nCases,      // in
    const int nPreds,         // in
    const int nMaxTerms,      // in
    const char* sPredNames[], // in: predictor names, can be NULL
    const bool HasWeights)    // in
{
    if(TraceGlobal == 1)
        printf("Forward pass term %d", IOFFSET);
    else if(TraceGlobal == 1.5)
        printf("Forward pass term %d\n", IOFFSET);
    else if(TraceGlobal >= 2) {
        int nMinSpan, nEndSpan, nStartSpan;
        GetSpanParams(&nMinSpan, &nEndSpan, &nStartSpan,
                      nCases, nPreds, 1 /*nDegree*/, 0 /*iParent*/, NULL /*bx*/);
        char sx[100];
        strcpy(sx,  sFormatMemSize(nCases * nPreds    * sizeof(double), false));
        char sbx[100];
        strcpy(sbx, sFormatMemSize(nCases * nMaxTerms * sizeof(double), false));
        printf("Forward pass: minspan %d endspan %d    x[%d,%d] %s    bx[%d,%d] %s%s\n\n",
            nMinSpan, nEndSpan,
            (int)nCases, nPreds, sx,
            (int)nCases, nMaxTerms,  sbx,
            HasWeights? " weighted": "");
        printf("         GRSq    RSq     DeltaRSq Pred ");
        if(sPredNames)
            printf("    PredName  ");
        printf("       Cut  Terms   Par Deg\n");

        // following matches printfs in PrintForwardStep
        if(sPredNames) // in: predictor names, can be NULL
            printf("%-4d%9.4f %6.4f                   %12.12s\n",
                IOFFSET, 0., 0., "(Intercept)");
        else
            printf("%-4d%9.4f %6.4f                               %d\n",
                IOFFSET, 0., 0., IOFFSET);
    }
}

//-----------------------------------------------------------------------------
static void PrintForwardStep(
    const int nTerms,
    const int nUsedTerms,
    const int iBestCase,
    const int iBestPred,
    const int iBestParent,
    const int nDegree,

    const double RSq,
    const double RSqDelta,
    const double Gcv,
    const double GcvNull,

    const size_t nCases,
    const int xOrder[],
    const double x[],
    const bool IsTermPair,
    const char* sPredNames[])   // in: predictor names, can be NULL
{
    if(TraceGlobal == 6)
        printf("\n\n");
    if(TraceGlobal == 1) {
        printf(", ");
        if(nTerms % 30 == 29)
            printf("\n     ");
        printf("%d", nTerms+IOFFSET);
    } else if(TraceGlobal == 1.5)
        printf("Forward pass term %d\n", nTerms+IOFFSET);
    else if(TraceGlobal >= 2) {
        tprintf(7,
"         GRSq    RSq     DeltaRSq Pred     PredName         Cut  Terms   Par Deg\n");
        printf("%-4d%9.4f %6.4f %12.4g ",
            nTerms+IOFFSET, 1 - Gcv / GcvNull, RSq, RSqDelta);
        if(iBestPred < 0) // *piBestCase not updated in FindKnot (no DeltaRSq)
            printf("   -                                       ");
        else {
            printf("%4d", iBestPred+IOFFSET);
            if(sPredNames) {
                if(sPredNames[iBestPred] && sPredNames[iBestPred][0])
                    printf(" %12.12s ", sPredNames[iBestPred]);
                else
                    printf(" %12.12s ", " ");
            }
            if(iBestCase == -1)
                printf("       none  ");
            else
                printf("% 11.5g%c ",
                    GetCut(iBestCase, iBestPred, nCases, x, xOrder),
                        (iBestCase==0? '<': ' ')); // print '<' if knot is at min
            if(IsTermPair) // two new used terms?
                printf("%-3d %-3d ", nUsedTerms-2+IOFFSET, nUsedTerms-1+IOFFSET);
            else
                printf("%-3d     ", nUsedTerms-1+IOFFSET);
            if(iBestParent != 0) // print parent if it isn't the intercept
                tprintf(2, "%3d ", iBestParent+IOFFSET);
            else
                tprintf(2, "    ", iBestParent+IOFFSET);
            printf("%3d ", nDegree);
        }
    }
#if !USING_R // no flush needed when using R_printf
    if(TraceGlobal != 0)
        fflush(stdout);
#endif
}

//-----------------------------------------------------------------------------
static int ForwardEpilog( // returns reason we stopped adding terms
    const int nTerms, const int nMaxTerms,
    const double Thresh,
    const double RSq, const double RSqDelta,
    const double Gcv, const double GcvNull,
    const int iBestCase,
    const bool FullSet[])
{
    tprintf(7,
"\n-----------------------------------------------------------------------------\n");

    const double GRSq = 1 - Gcv / GcvNull;
    int iTermCond = 0;

    char sUsed[100] = "";
    const int nUsed = GetNbrUsedCols(FullSet, nMaxTerms);
    if(nUsed != nTerms)
        sprintf(sUsed, ", %d term%s used", nUsed, nUsed == 1? "": "s");
    char sTerms[100];
    sprintf(sTerms, "%d term%s%s", nTerms, nTerms == 1? "": "s", sUsed);

    // NOTE 1: this code must match the loop termination conditions in ForwardPass
    // NOTE 2: if you update this, also update print.termcond in the R code

    // treat very low nMaxTerms as a special case
    // because RSDelta etc. not yet completely initialized
    if(nMaxTerms < 3) {
        iTermCond = 1;
        tprintf(1, "\nReached maximum number of terms %d\n", nMaxTerms);
    } else if(Thresh != 0 && GRSq < MIN_GRSQ) {
        if(GRSq < -1000) {
            iTermCond = 2;
            tprintf(1, "\nGRSq -Inf at %s\n", sTerms);
        } else {
            iTermCond = 3;
            if(TraceGlobal >= 1)
                printf("\nReached minimum GRSq %g at %s (GRSq %.2g)\n",
                       MIN_GRSQ, sTerms, GRSq);
        }
    } else if(Thresh != 0 && RSqDelta < Thresh) {
        iTermCond = 4;
        if(TraceGlobal >= 1)
            printf("\nRSq changed by less than %g at %s (DeltaRSq %.2g)\n",
                   Thresh, sTerms, RSqDelta);
    } else if(RSq >= 1-Thresh) {
        iTermCond = 5;
        if(TraceGlobal >= 1)
            printf("\nReached maximum RSq %.4f at %s (RSq %.4f)\n",
                   1-Thresh, sTerms, RSq);
    } else if(iBestCase < 0) { // TODO seems fishy, happens with linpreds so should give appropriate msg?
        iTermCond = 6;
        tprintf(1,
                "\nNo new term increases RSq (perhaps reached numerical limits) at %s\n",
                sTerms);
    } else {
        iTermCond = 7;
#if USING_R
        tprintf(1, "\nReached nk %d\n", nMaxTerms);
#else
        tprintf(1, "\nReached maximum number of terms %d\n", nMaxTerms);
#endif
    }
    if(TraceGlobal >= 1)
        printf("After forward pass GRSq %.3f RSq %.3f\n", GRSq, RSq);
    tprintf(2, "Forward pass complete: %s\n", sTerms);
    tprintf(3, "\n");

    return iTermCond;
}

//-----------------------------------------------------------------------------
static void CheckVec(
    const double x[],
    const size_t nCases,
    const int nCols,
    const char sVecName[])
{
    int iCol, i;

    for(iCol = 0; iCol < nCols; iCol++)
        for(i = 0; i < (const int)nCases; i++) {
#if USING_R
             if(ISNA(x[i + iCol * nCases])) {
                 if(nCols > 1)
                     error("%s[%d,%d] is NA",
                         sVecName, i+IOFFSET, iCol+IOFFSET);
                 else
                     error("%s[%d] is NA", sVecName, i+IOFFSET);
             }
#endif
             if(ISNAN(x[i + iCol * nCases])) {
                 if(nCols > 1)
                     error("%s[%d,%d] is NaN",
                         sVecName, i+IOFFSET, iCol+IOFFSET);
                 else
                     error("%s[%d] is NaN", sVecName, i+IOFFSET);
             }
             if(!FINITE(x[i + iCol * nCases])) {
                 if(nCols > 1)
                     error("%s[%d,%d] is not finite",
                         sVecName, i+IOFFSET, iCol+IOFFSET);
                 else
                     error("%s[%d] is not finite", sVecName, i+IOFFSET);
             }
    }
}

//-----------------------------------------------------------------------------
static double CheckRssNull(
    double RssNull,
    const double y[],
    const int iResp,
    const int nResp,
    const size_t nCases)
{
    if(RssNull < 1e-8 * nCases) {   // 1e-8 is arbitrary
        if(nResp)
            tprintf(1, "Variance of y[,%d] is zero (values are all equal to %g)\n",
                    iResp+IOFFSET, y_(0,iResp));
        else
            tprintf(1, "Variance of y is zero (values are all equal to %g)\n",
                    y_(0,iResp));
        RssNull = 1e-8 * nCases; // prevent later divide by zero
    }
    return RssNull;
}

//-----------------------------------------------------------------------------
static double GetRssNull(
    const double y[],          // in: nCases x nResp, unweighted but scaled y
    const double WeightsArg[], // in: nCases x 1, can be NULL
    const size_t nCases,       // in: number of rows in x and elements in y
    const int nResp)           // in: number of cols in y
{
    double RssNull = 0;
    if(WeightsArg)
        for(int iResp = 0; iResp < nResp; iResp++) {
            double SumY = 0, SumWeights = 0;
            int i;
            for(i = 0; i < (const int)nCases; i++) {
                SumY += WeightsArg[i] * y_(i,iResp);
                SumWeights += WeightsArg[i];
            }
            const double WeightedMean = SumY / SumWeights;
            for(i = 0; i < (const int)nCases; i++) {
                RssNull += WeightsArg[i] * sq(y_(i,iResp) - WeightedMean);
            }
            RssNull = CheckRssNull(RssNull, y, iResp, nResp, nCases);
        }
    else // no weights
        for(int iResp = 0; iResp < nResp; iResp++) {
            const double yMean = Mean(&y_(0,iResp), nCases);
            RssNull += SumOfSquares(&y_(0,iResp), yMean, nCases);
            RssNull = CheckRssNull(RssNull, y, iResp, nResp, nCases);
        }

    return RssNull;
}

//-----------------------------------------------------------------------------
// The limits below are somewhat arbitrary and generous.  They are intended to
// catch gross errors on the part of the caller, and to prevent crashes because
// of 0 sizes etc.  We use error rather than ASSERT because these are user
// settable params and we want to be informative from the user's perspective.
// The errors are reported using the variable names in the R code.

static void CheckForwardPassArgs(
    const double x[],
    const double y[],
    const double yw[],
    const double WeightsArg[],
    const size_t nCases,
    const int nResp,
    const int nPreds,
    const int nMaxDegree,
    const int nMaxTerms,
    const double Penalty,
    const double Thresh,
    const double FastBeta,
    const double NewVarPenalty,
    const int  LinPreds[],
    const double AdjustEndSpan,
    const bool AutoLinPreds,
    const bool UseBetaCache)
{
    const int nCases1 = (const int)nCases; // type convert from size_t
    if(nCases1 < 2)
        error("the x matrix must have at least two rows");
    if(nCases1 > 1e9)
        error("too many rows %d in the input matrix, max allowed is 1e9", nCases1);
    if(nResp < 1)
        error("the number of responses %d is less than 1", nResp);
        if(nPreds < 1)
        error("the number of predictors %d is less than 1", nPreds);
    if(nPreds > 1e5)
        error("the number of predictors %d is greater than 1e5", nPreds);
    if(nMaxDegree <= 0)
        error("degree %d is not greater than 0", nMaxDegree);
    if(nMaxDegree > MAX_DEGREE)
        error("degree %d is greater than %d", nMaxDegree, MAX_DEGREE);
    if(nMaxTerms < 1)       // prevent internal misbehavior
        error("nk %d is less than 1", nMaxTerms);
    if(nMaxTerms > 1000)
        error("nk %d is greater than 1000", nMaxTerms);
    if(Penalty < 0 && Penalty != -1)
        error("penalty %g is less than 0 and the only legal value less than 0 is -1 "
            "(meaning terms and knots are free)", Penalty);
    if(Penalty > 1000)
        error("penalty %g is greater than 1000", Penalty);
    if(Thresh < 0)
        error("thresh %g is less than 0", Thresh);
    if(Thresh >= 1)
        error("thresh %g >= 1", Thresh);
    if(nMinSpanGlobal > nCases1)
        error("minspan %d is greater than the number of cases %d", nMinSpanGlobal, nCases1);
    if(nEndSpanGlobal > nCases1)
        error("endspan %d is greater than the number of cases %d", nEndSpanGlobal, nCases1);
    else if(nEndSpanGlobal < 0)
        error("endspan %d is less than 0", nEndSpanGlobal);
    if(FastBeta < 0)
        error("fast.beta %g is less than 0", FastBeta);
    if(FastBeta > 1000)
        error("fast.beta %g is greater than 1000", FastBeta);
    if(TraceGlobal < 0)
        warning("trace %g is less than 0", TraceGlobal);
    if(TraceGlobal > 10)
        warning("trace %g is greater than 10", TraceGlobal);
    if(NewVarPenalty < 0)
        warning("newvar.penalty %g is less than 0", NewVarPenalty);
    if(NewVarPenalty > 100)
        warning("newvar.penalty %g is greater than 100", NewVarPenalty);
    if(AdjustEndSpan < 0 || AdjustEndSpan > 10)
        error("Endspan.penalty is %g but should be between 0 and 10", AdjustEndSpan);
    if(AutoLinPreds != 0 && AutoLinPreds != 1)
        error("Auto.linpreds  is neither TRUE nor FALSE");
    if(UseBetaCache != 0 && UseBetaCache != 1)
        warning("Use.Beta.Cache is neither TRUE nor FALSE");

    CheckVec(x, nCases, nPreds, "x");
    CheckVec(y, nCases, nResp,  "y");
#if WEIGHTS
    if(yw) {
        ASSERT(WeightsArg);
        CheckVec(yw, nCases, nResp,  "yw");
        for(int i = 0; i < (const int)nCases; i++) {
            CheckVec(WeightsArg, nCases, 1, "weights");
            if(WeightsArg[i] < ALMOST_ZERO)
                error("weights[%d] is not greater than zero",
                      i+IOFFSET);
        }
    }
#else
    ASSERT(!yw);
    ASSERT(!WeightsArg || abs(WeightsArg[0] - WeightsArg[1]) < 1e-8);
#endif
    for(int iPred = 0; iPred < nPreds; iPred++)
        if(LinPreds[iPred] != 0 && LinPreds[iPred] != 1)
            error("linpreds[%d] is not 0 or 1", iPred+IOFFSET);
}

//-----------------------------------------------------------------------------
// Forward pass
//
// After initializing the intercept term, the main for loop adds terms in pairs.
// In the for loop, nTerms is the index of the potential new term; nTerms+1
// the index of its partner.
// The upper term in the term pair may not be useable.  If so we still
// increment nTerms by 2 but don't set the flag in FullSet.

static void ForwardPass(
    int*   pnTerms,              // out: highest used term number in full model
    int*   piTermCond,           // out: reason we terminated the forward pass
    bool   FullSet[],            // out: 1 * nMaxTerms, indices of lin indep cols of bx
    double bx[],                 // out: MARS basis matrix, nCases * nMaxTerms
    int    Dirs[],               // out: nMaxTerms * nPreds, -1,0,1,2 for iTerm, iPred
    double Cuts[],               // out: nMaxTerms * nPreds, cut for iTerm, iPred
    int    nDegree[],            // out: degree of each term, degree of intercept is 0
    int    nUses[],              // out: nbr of times each predictor is used in the model
    const double x[],            // in: nCases x nPreds, unweighted x
    const double y[],            // in: nCases x nResp, unweighted but scaled y
    const double yw[],           // in: nCases x nResp, weighted and scaled y, can be NULL
    const double WeightsArg[],   // in: nCases x 1, can be NULL
    const size_t nCases,         // in: number of rows in x and elements in y
    const int nResp,             // in: number of cols in y
    const int nPreds,            // in:
    const int nMaxDegree,        // in:
    const int nMaxTerms,         // in:
    const double Penalty,        // in: GCV penalty per knot
    const double Thresh,         // in: forward step threshold
    int nFastK,                  // in: Fast MARS K
    const double FastBeta,       // in: Fast MARS ageing coef
    const double NewVarPenalty,  // in: penalty for adding a new variable (default is 0)
    const int  LinPreds[],       // in: nPreds x 1, 1 if predictor must enter linearly
    const double AdjustEndSpan,  // in:
    const bool AutoLinPreds,     // in: assume predictor linear if knot is min predictor value
    const bool UseBetaCache,     // in: true to use the beta cache, for speed
    const char* sPredNames[])    // in: predictor names, can be NULL
{
    tprintf(5, "earth.c %s\n", VERSION);
    CheckForwardPassArgs(x, y, yw, WeightsArg, nCases, nResp, nPreds,
        nMaxDegree, nMaxTerms, Penalty, Thresh, FastBeta, NewVarPenalty,
        LinPreds, AdjustEndSpan, AutoLinPreds, UseBetaCache);

    if(nFastK <= 0)
        nFastK = 10000+1;   // bigger than any nMaxTerms
    if(nFastK < 3)          // avoid possible queue boundary conditions
        nFastK = 3;

    xOrder = GetArrayOrder(x, nCases, nPreds);
    InitBetaCache(UseBetaCache, nMaxTerms, nPreds);

    bxOrth          = (double*)malloc1(nCases * nMaxTerms * sizeof(double),
                        "bxOrth\t\tnCases %d nMaxTerms %d sizeof(double) %d",
                        (const int)nCases, nMaxTerms, sizeof(double));

    bxOrthCenteredT = (double*)malloc1(nMaxTerms * nCases * sizeof(double),
                        "bxOrthCenteredT\tnMaxTerms %d nCases %d sizeof(double) %d",
                        nMaxTerms, (const int)nCases, sizeof(double));

    bxOrthMean      = (double*)malloc1(nMaxTerms * nResp * sizeof(double),
                        "bxOrthMean\t\tnMaxTerms %d nResp %d sizeof(double) %d",
                        nMaxTerms, nResp, sizeof(double));

    yMean           = (double*)malloc1(nResp * sizeof(double),
                        "yMean\t\t\tnResp %d sizeof(double) %d",
                        nResp, sizeof(double));

    memset(FullSet, 0, nMaxTerms * sizeof(bool));
    memset(Dirs,    0, nMaxTerms * nPreds * sizeof(int));
    memset(Cuts,    0, nMaxTerms * nPreds * sizeof(double));
    memset(nDegree, 0, nMaxTerms * sizeof(int));
    memset(nUses,   0, nPreds * sizeof(int));
    memset(bx,      0, nCases * nMaxTerms * sizeof(double));

    // Intercept columns of bx and bxOrth.
    // Note that we use the weights here, and since every term is a multiple
    // of this intercept term or a term derived from it, we don't need to use
    // the weighted x when forming bx.
    if(WeightsArg)
        for(int i = 0; i < (const int)nCases; i++)
            bx_(i,0) = sqrt(WeightsArg[i]);
    else
        for(int i = 0; i < (const int)nCases; i++)
            bx_(i,0) = 1;

    bool GoodCol;
    InitBxOrthCol(bxOrth, bxOrthCenteredT, bxOrthMean, &GoodCol,
        &bx_(0,0), 0 /*nTerms*/, FullSet, nCases, nMaxTerms, -1, -1);
    if(!GoodCol) // should never happen
        tprintf(1, "GoodCol is false in ForwardPass\n");
    GoodCol = true;
    FullSet[0] = true;  // intercept

    for(int iResp = 0; iResp < nResp; iResp++)
        yMean[iResp]  = Mean(&y_(0,iResp), nCases);
    const double RssNull = GetRssNull(y, WeightsArg, nCases, nResp);
    double Rss = RssNull, RssDelta = RssNull, RSq = 0, RSqDelta = 0;
    int nUsedTerms = 1;     // number of used basis terms including intercept, for GCV calc
    double Gcv = 0, GcvNull = GetGcv(nUsedTerms, nCases, RssNull, Penalty);
    PrintForwardProlog(nCases, nPreds, nMaxTerms, sPredNames, yw != NULL);
#if FAST_MARS
    InitQ(nMaxTerms);
    AddTermToQ(0, 1, RssNull, true, nMaxTerms, FastBeta); // intercept term into Q
#endif
    int nTerms = 1, iBestCase = -1;
    if(nMaxTerms >= 3) while(1) { // start after intercept, add terms in pairs
        int iBestPred = -1, iBestParent = -1;
        bool IsNewForm, LinPredIsBest;
#if USING_R
        ServiceR();
#endif
        if(Rss <= 0)
            error("assertion failed: Rss <= 0 (y is all const?)");
        ASSERT(RssDelta > 0);
        // Changed factor from 2 to 10 in version 4.2.0 (2 was too conservative).
        // Note that only the code without weights uses this.
        const double MaxLegalRssDelta = min(1.01 * Rss, 10 * RssDelta);

        FindTerm(&iBestCase, &iBestPred, &iBestParent,
            &RssDelta, &IsNewForm, &LinPredIsBest,
            bxOrth, bxOrthCenteredT, bxOrthMean,
            bx, x, y, yw, nCases, nResp, nPreds, nMaxDegree, nTerms, nMaxTerms,
            yMean, Rss, MaxLegalRssDelta, FullSet, xOrder, nDegree,
            nUses, Dirs, nFastK, NewVarPenalty, LinPreds);

        // following code added for Auto.linpreds (earth version 4.6.0, Dec 2017)
        if((LinPredIsBest && iBestCase != 0) || // paranoia, should never happen
           (!LinPredIsBest && iBestCase == 0))
            printf("\nLinPredIsBest %d yet iBestCase%-5d\n", LinPredIsBest, iBestCase);
        if(!AutoLinPreds && !LinPreds[iBestPred])
            LinPredIsBest = false;

        if(iBestCase >= 0)
            AddTermPair(Dirs, Cuts, bx, bxOrth, bxOrthCenteredT, bxOrthMean,
                FullSet, &IsNewForm, nDegree, nUses,
                nTerms, iBestParent, iBestCase, iBestPred, nPreds, nCases,
                nMaxTerms, LinPredIsBest, LinPreds, x, xOrder, yw != NULL);

        const bool IsTermPair = iBestCase > 0 && IsNewForm;
        nUsedTerms++;
        if(IsTermPair)
            nUsedTerms++;     // add paired term
        Rss -= RssDelta;
        Rss = MaybeZero(Rss); // RSS can go slightly neg due to numerical error
        Gcv = GetGcv(nUsedTerms, nCases, Rss, Penalty);
        const double OldRSq = RSq;
        RSq = 1 - Rss / RssNull;
        RSqDelta = MaybeZero(RSq - OldRSq);

        PrintForwardStep(nTerms, nUsedTerms, iBestCase, iBestPred,
            iBestParent, iBestParent < 0? 0: nDegree[iBestParent]+1,
            RSq, RSqDelta, Gcv, GcvNull,
            nCases, xOrder, x, IsTermPair, sPredNames);

        // note the possible breaks in the code below

        if(iBestCase < 0) { // *piBestCase was not updated in FindKnot
            FullSet[nTerms] = FullSet[nTerms+1] = false;
            tprintf(2, "reject (no DeltaRsq)\n");
            break;
        }
        if(Thresh != 0 && RSqDelta < Thresh) {
            FullSet[nTerms] = FullSet[nTerms+1] = false;
            tprintf(2, "reject (small DeltaRSq)\n");
            break;
        }
        if(Thresh != 0 && 1 - Gcv / GcvNull < MIN_GRSQ) {
            FullSet[nTerms] = FullSet[nTerms+1] = false;
            tprintf(2, "reject (negative GRSq)\n");
            break;
        }
#if FAST_MARS
        if(!LinPredIsBest && IsNewForm) {   // good upper term?
            AddTermToQ(nTerms,   nTerms, POS_INF, false, nMaxTerms, FastBeta);
            AddTermToQ(nTerms+1, nTerms, POS_INF, true,  nMaxTerms, FastBeta);
        } else
            AddTermToQ(nTerms,   nTerms, POS_INF, true,  nMaxTerms, FastBeta);
        if(TraceGlobal == 6)
            PrintSortedQ(nFastK);
#endif
        nTerms += 2;
        if(RSq >= 1 - Thresh) {
            tprintf(2, "final (max RSq)\n");
            break;
        }
        if(nTerms >= nMaxTerms - 1) { // -1 allows for upper term in pair
            tprintf(2, "final (reached nk %d)\n", nMaxTerms);
            break;
        }
        tprintf(2, "\n");
    } // while(1)
    *piTermCond = ForwardEpilog(nTerms, nMaxTerms, Thresh, RSq, RSqDelta,
                                Gcv, GcvNull, iBestCase, FullSet);
    *pnTerms = nTerms;

    // free in opposite order to alloc to help operating system memory manager
#if FAST_MARS
    FreeQ();
#endif
    free1(yMean);
    free1(bxOrthMean);
    free1(bxOrthCenteredT);
    free1(bxOrth);
    FreeBetaCache();
    free1(xOrder);
}

//-----------------------------------------------------------------------------
// This is an interface from R to the C routine ForwardPass

#if USING_R
void ForwardPassR(              // for use by R
    int    FullSet[],           // out: nMaxTerms x 1, bool vec of lin indep cols of bx
    double bx[],                // out: MARS basis matrix, nCases x nMaxTerms
    double Dirs[],              // out: nMaxTerms x nPreds, elements are -1,0,1,2
    double Cuts[],              // out: nMaxTerms x nPreds, cut for iTerm,iPred
    int*   piTermCond,          // out: reason we terminated the forward pass
    const double x[],           // in: nCases x nPreds, unweighted x
    const double y[],           // in: nCases x nResp, unweighted but scaled y
    const double yw[],          // in: nCases x nResp, weighted and scaled y
    const double WeightsArg[],  // in: nCases x 1, never MyNullDouble
    const int* pnCases,         // in: number of rows in x and elements in y
    const int* pnResp,          // in: number of cols in y
    const int* pnPreds,         // in: number of cols in x
    const int* pnMaxDegree,     // in:
    const double* pPenalty,     // in:
    const int* pnMaxTerms,      // in:
    const double* pThresh,      // in: forward step threshold
    const int* pnMinSpan,       // in:
    const int* pnEndSpan,       // in:
    const int* pnFastK,         // in: Fast MARS K
    const double* pFastBeta,    // in: Fast MARS ageing coef
    const double* pNewVarPenalty,  // in: penalty for adding a new variable (default is 0)
    const int  LinPreds[],         // in: nPreds x 1, 1 if predictor must enter linearly
    const SEXP Allowed,            // in: constraints function, can be MyNullFunc
    const int* pnAllowedFuncArgs,  // in: number of arguments to Allowed function, 3 or 4
    const SEXP Env,                // in: environment for Allowed function
    const double* pAdjustEndSpan,  // in:
    const int* pnAutoLinPreds,     // in: assume predictor linear if knot is min predictor value
    const int* pnUseBetaCache,     // in: 1 to use the beta cache, for speed
    const double* pTrace,          // in: 0 none 1 overview 2 forward 3 pruning 4 more pruning
    const char* sPredNames[],      // in: predictor names in trace printfs
    const double* MyNullDouble,    // in: trick to avoid R check warnings on passing R_NilValue
    const SEXP MyNullFunc)         // in: trick to avoid R check warnings on passing R_NilValue
{
    TraceGlobal = *pTrace;
    nMinSpanGlobal = *pnMinSpan;
    nEndSpanGlobal = *pnEndSpan;
    AdjustEndSpanGlobal = *pAdjustEndSpan;

    size_t nCases = *pnCases; // type convert

    const int nResp = *pnResp;
    const int nPreds = *pnPreds;
    const int nMaxTerms = *pnMaxTerms;

    // nUses is the number of time each predictor is used in the model
    nUses = (int*)malloc1(*pnPreds * sizeof(int),
                    "nUses\t\t\t*pnPreds %d sizeof(int)",
                    *pnPreds, sizeof(int));

    //  nDegree is degree of each term, degree of intercept is considered to be 0
    nDegree = (int*)malloc1(nMaxTerms * sizeof(int),
                        "nDegree\t\tnMaxTerms %d sizeof(int) %d",
                        nMaxTerms, sizeof(int));

    iDirs = (int*)calloc1(nMaxTerms * nPreds, sizeof(int),
                        "iDirs\t\t\tnMaxTerms %d nPreds %d sizeof(int) %d",
                        nMaxTerms, nPreds, sizeof(int));

    // convert int to bool (may be redundant, depending on compiler)
    BoolFullSet = (bool*)malloc1(nMaxTerms * sizeof(bool),
                        "BoolFullSet\t\tnMaxTerms %d sizeof(bool) %d",
                        nMaxTerms, sizeof(bool));

    int iTerm;
    for(iTerm = 0; iTerm < nMaxTerms; iTerm++)
        BoolFullSet[iTerm] = FullSet[iTerm] != 0;

    // convert R my.null to C NULL

#if !WEIGHTS
    ASSERT(*(double*)yw == *(double*)MyNullDouble);
#endif
    if(*(double*)yw == *(double*)MyNullDouble)
        yw = NULL;
    ASSERT(*(double*)WeightsArg != *(double*)MyNullDouble);

    InitAllowedFunc(*(int*)Allowed == *(int*)MyNullFunc? NULL: Allowed,
                    *pnAllowedFuncArgs, Env, sPredNames, nPreds);

    int nTerms;
    ForwardPass(&nTerms, piTermCond,
            BoolFullSet, bx, iDirs, Cuts, nDegree, nUses,
            x, y, yw, WeightsArg, nCases, nResp, nPreds, *pnMaxDegree, nMaxTerms,
            *pPenalty, *pThresh, *pnFastK, *pFastBeta, *pNewVarPenalty,
            LinPreds, *pAdjustEndSpan, *pnAutoLinPreds,
            (bool)(*pnUseBetaCache != 0), sPredNames);

    FreeAllowedFunc();

    // remove linearly independent columns if necessary -- this updates BoolFullSet

    RegressAndFix(NULL, NULL, NULL, BoolFullSet,
        bx, yw? yw: y, nCases, nResp, nMaxTerms);

    for(iTerm = 0; iTerm < nMaxTerms; iTerm++)      // convert int to double
        for(int iPred = 0; iPred < nPreds; iPred++)
            Dirs[iTerm + iPred * nMaxTerms] =
                iDirs[iTerm + iPred * nMaxTerms];

    for(iTerm = 0; iTerm < nMaxTerms; iTerm++)      // convert bool to int
        FullSet[iTerm] = BoolFullSet[iTerm];

    free1(BoolFullSet);
    free1(iDirs);
    free1(nDegree);
    free1(nUses);
}
#endif // USING_R

//-----------------------------------------------------------------------------
// Step backwards through the terms, at each step deleting the term that
// causes the least RSS increase.  The subset of terms and RSS of each subset are
// saved in PruneTerms and RssVec (which are indexed on subset size).
//
// The crux of the method used here is that the change in RSS (for nResp=1)
// caused by removing predictor iPred is DeltaRss = sq(Betas[iPred]) / Diags[iPred]
// where Diags is the diagonal elements of the inverse of X'X.
// See for example Miller (see refs in file header) section 3.4 p44.
//
// For multiple responses we sum the above DeltaRss over all responses.
//
// This method is fast and simple but accuracy can be poor if inv(X'X) is
// ill conditioned.  The Alan Miller code in the R package "leaps" uses a more
// stable method, but does not support multiple responses.
//
// The "Xtx" in the name refers to the X'X matrix.

static void EvalSubsetsUsingXtx(
    bool   PruneTerms[],       // out: nMaxTerms x nMaxTerms
    double RssVec[],           // out: nMaxTerms x 1, RSS of each subset
    const size_t nCases,       // in
    const int    nResp,        // in: number of cols in y
    const int    nMaxTerms,    // in: number of MARS terms in full model
    const double bx[],         // in: nCases x nMaxTerms, all cols must be indep
    const double y[])          // in: nCases * nResp
{
    double* Betas = (double*)malloc1(nMaxTerms * nResp * sizeof(double),
                        "Betas\t\t\tnMaxTerms %d nResp %d sizeof(double) %d",
                        nMaxTerms, nResp, sizeof(double));

    double* Diags = (double*)malloc1(nMaxTerms * sizeof(double),
                        "Diags\t\t\tnMaxTerms %d sizeof(double) %d",
                        nMaxTerms, sizeof(double));

    WorkingSet = (bool*)malloc1(nMaxTerms * sizeof(bool),
                        "WorkingSet\t\tnMaxTerms %d sizeof(bool) %d",
                        nMaxTerms, sizeof(bool));

    for(int iTerm = 0; iTerm < nMaxTerms; iTerm++)
        WorkingSet[iTerm] = true;
    const double RssNull = GetRssNull(y, NULL, nCases, nResp);
    tprintf(4, "EvalSubsetsUsingXtx:\nnTerms iTerm DeltaRss    RSq\n");

    for(int nUsedCols = nMaxTerms; nUsedCols > 0; nUsedCols--) {
        int nRank;
        double Rss;
        Regress(Betas, NULL, &Rss, Diags, &nRank, NULL,
            bx, y, nCases, nResp, nMaxTerms, WorkingSet);

        if(nRank != nUsedCols)
            error("nRank %d != nUsedCols %d "
                "(probably because of lin dep terms in bx)\n",
                nRank, nUsedCols);

        RssVec[nUsedCols-1] = Rss;
        memcpy(PruneTerms + (nUsedCols-1) * nMaxTerms, WorkingSet,
               nMaxTerms * sizeof(bool));

        if(nUsedCols == 1)
            break;

        // set iDelete to the best term for deletion

        int iDelete = -1;   // term to be deleted
        int iTerm1 = 0;     // index taking into account false vals in WorkingSet
        double MinDeltaRss = POS_INF;
        for(int iTerm = 0; iTerm < nMaxTerms; iTerm++) {
            if(WorkingSet[iTerm]) {
                double DeltaRss = 0;
                for(int iResp = 0; iResp < nResp; iResp++)
                    DeltaRss += sq(Betas_(iTerm1, iResp)) / Diags[iTerm1];
                bool NewMin = false;
                if(iTerm > 0 && DeltaRss < MinDeltaRss) {   // new minimum?
                    MinDeltaRss = DeltaRss;
                    iDelete = iTerm;
                    NewMin = true;
                }
                if(iTerm != 0)
                    tprintf(4, "%6d %5d %8.5g %6.4f%s\n",
                            nUsedCols, iTerm+IOFFSET, DeltaRss,
                            1 - (Rss + DeltaRss) / RssNull, NewMin? " min" : "");
                iTerm1++;
            }
        }
        ASSERT(iDelete > 0);
        WorkingSet[iDelete] = false;
        tprintf(4, "\n");
    }
    free1(WorkingSet);
    free1(Diags);
    free1(Betas);
}

//-----------------------------------------------------------------------------
// This is invoked from R if y has multiple columns i.e. a multiple response model.
// It is needed because the alternative (Alan Miller's Fortran code) supports
// only one response.

#if USING_R
void EvalSubsetsUsingXtxR(      // for use by R
    double       PruneTerms[],  // out: specifies which cols in bx are in best set
    double       RssVec[],      // out: nTerms x 1
    const int*   pnCases,       // in
    const int*   pnResp,        // in: number of cols in y
    const int*   pnMaxTerms,    // in
    const double bx[],          // in: MARS basis matrix, all cols must be indep
    const double y[],           // in: nCases * nResp (possibly weighted)
    const double* pTrace)       // in
{
    TraceGlobal = *pTrace;

    const int nMaxTerms = *pnMaxTerms;
    bool* BoolPruneTerms = (bool*)malloc1(nMaxTerms * nMaxTerms * sizeof(bool),
                                "BoolPruneTerms\tMaxTerms %d nMaxTerms %d sizeof(bool) %d",
                                nMaxTerms, nMaxTerms, sizeof(bool));

    size_t nCases = *pnCases; // type convert

    EvalSubsetsUsingXtx(BoolPruneTerms, RssVec, nCases, *pnResp,
                        nMaxTerms, bx, y);

    // convert BoolPruneTerms to upper triangular matrix PruneTerms

    for(int iModel = 0; iModel < nMaxTerms; iModel++) {
        int iPrune = 0;
        for(int iTerm = 0; iTerm < nMaxTerms; iTerm++)
            if(BoolPruneTerms[iTerm + iModel * nMaxTerms])
                PruneTerms[iModel + iPrune++ * nMaxTerms] = iTerm + 1;
    }
    free1(BoolPruneTerms);
}
#endif // USING_R

//-----------------------------------------------------------------------------
#if STANDALONE && WEIGHTS
static void UnweightBx(
    double bx[],               // in: nCases x nMaxTerms
    const double WeightsArg[], // in
    const size_t nCases,       // in: number of rows in bx
    const int nMaxTerms)       // in: number of cols in bx
{
    if(WeightsArg) {
        for(int iTerm = 0; iTerm < nMaxTerms; iTerm++)
            for(int i = 0; i < (const int)nCases; i++)
                bx_(i, iTerm) /= sqrt(WeightsArg[i]);
    }
}
#endif // STANDALONE

//-----------------------------------------------------------------------------
#if STANDALONE
static void BackwardPass(
    double* pBestGcv,          // out: GCV of the best model i.e. BestSet columns of bx
    bool   BestSet[],          // out: nMaxTerms x 1, indices of best set of cols of bx
    double Residuals[],        // out: nCases x nResp
    double Betas[],            // out: nMaxTerms x nResp
    double bx[],               // in: nCases x nMaxTerms, will be unweighted if weights
    const double y[],          // in: nCases x nResp
    const double* WeightsArg,  // in; NULL or nCases
    const size_t nCases,       // in: number of rows in bx and elements in y
    const int nResp,           // in: number of cols in y
    const int nMaxTerms,       // in: number of cols in bx
    const double Penalty)      // in: GCV penalty per knot
{
    double* RssVec = (double*)malloc1(nMaxTerms * sizeof(double),
                        "RssVec\t\tnMaxTerms %d sizeof(double) %d",
                        nMaxTerms, sizeof(double));

    bool* PruneTerms = (bool*)malloc1(nMaxTerms * nMaxTerms * sizeof(bool),
                        "PruneTerms\t\tnMaxTerms %d nMaxTerms %d sizeof(bool) %d",
                        nMaxTerms, nMaxTerms, sizeof(bool));

    EvalSubsetsUsingXtx(PruneTerms, RssVec, nCases, nResp,
                        nMaxTerms, bx, y);

    // now we have the RSS for each model, so find the iModel which has the best GCV

    tprintf(3, "Backward pass:\nSubsetSize         GRSq          RSq\n");
    int iBestModel = -1;
    double GcvNull = GetGcv(1, nCases, RssVec[0], Penalty);
    double BestGcv = POS_INF;
    for(int iModel = 0; iModel < nMaxTerms; iModel++) {
        const double Gcv = GetGcv(iModel+1, nCases, RssVec[iModel], Penalty);
        if(Gcv < BestGcv) {
            iBestModel = iModel;
            BestGcv = Gcv;
        }
        double GRSq = 1 - BestGcv/GcvNull;
        // Prevent negative almost-zero issued by some compilers (earth version 4.4.5)
        // This prints as -0.0000 below and messes up the test scripts.
        if(GRSq > -1e-12 && GRSq < 0) // -1e-12 is fairly arb
            GRSq = 0;
        tprintf(3, "%10d %12.4f %12.4f\n", iModel+IOFFSET,
                GRSq, 1 - RssVec[iModel]/RssVec[0]);
    }
    tprintf(3, "\nBackward pass complete: selected %d terms of %d, "
            "GRSq %.3f RSq %.3f\n\n",
            iBestModel+IOFFSET, nMaxTerms,
            1 - BestGcv/GcvNull, 1 - RssVec[iBestModel]/RssVec[0]);

    // set BestSet to the model which has the best GCV

    ASSERT(iBestModel >= 0);
    memcpy(BestSet, PruneTerms + iBestModel * nMaxTerms, nMaxTerms * sizeof(bool));
    free1(PruneTerms);
    free1(RssVec);
    *pBestGcv = BestGcv;

#if WEIGHTS
    UnweightBx(bx, WeightsArg, nCases, nMaxTerms);
    // TODO should use weighted regression in RegressAndFix below
#endif

    // get final model Betas, Residuals, Rss
    RegressAndFix(Betas, Residuals, NULL, BestSet,
        bx, y, nCases, nResp, nMaxTerms);
}
#endif // STANDALONE

//-----------------------------------------------------------------------------
#if STANDALONE
static int DiscardUnusedTerms(
    double bx[],             // io: nCases x nMaxTerms
    int    Dirs[],           // io: nMaxTerms x nPreds
    double Cuts[],           // io: nMaxTerms x nPreds
    bool   WhichSet[],       // io: tells us which terms to discard
    int    nDegree[],        // io: degree of each term, degree of intercept is 0
    const int nMaxTerms,
    const int nPreds,
    const size_t nCases)
{
    int nUsed = 0, iTerm;
    for(iTerm = 0; iTerm < nMaxTerms; iTerm++)
        if(WhichSet[iTerm]) {
            memcpy(bx + nUsed * nCases, bx + iTerm * nCases, nCases * sizeof(double));
            for(int iPred = 0; iPred < nPreds; iPred++) {
                Dirs_(nUsed, iPred) = Dirs_(iTerm, iPred);
                Cuts_(nUsed, iPred) = Cuts_(iTerm, iPred);
            }
            nDegree[nUsed] = nDegree[iTerm];
            nUsed++;
        }
    memset(WhichSet, 0, nMaxTerms * sizeof(bool));
    for(iTerm = 0; iTerm < nUsed; iTerm++)
        WhichSet[iTerm] = true;
    return nUsed;
}
#endif // STANDALONE

//-----------------------------------------------------------------------------
#if STANDALONE
void Earth(
    double* pBestGcv,            // out: GCV of the best model i.e. BestSet columns of bx
    int*    pnTerms,             // out: max term nbr in final model, after removing lin dep terms
    int*    piTermCond,          // out: reason we terminated the foward pass
    bool    BestSet[],           // out: nMaxTerms x 1, indices of best set of cols of bx
    double  bx[],                // out: nCases x nMaxTerms
    int     Dirs[],              // out: nMaxTerms x nPreds, -1,0,1,2 for iTerm, iPred
    double  Cuts[],              // out: nMaxTerms x nPreds, cut for iTerm, iPred
    double  Residuals[],         // out: nCases x nResp
    double  Betas[],             // out: nMaxTerms x nResp
    const double x[],            // in: nCases x nPreds
    const double y[],            // in: nCases x nResp
    const double WeightsArg[],   // in: nCases x 1, can be NULL, not yet supported
    const size_t nCases,         // in: number of rows in x and elements in y
    const int nResp,             // in: number of cols in y
    const int nPreds,            // in: number of cols in x
    const int nMaxDegree,        // in: Friedman's mi
    const int nMaxTerms,         // in: includes the intercept term
    const double Penalty,        // in: GCV penalty per knot
    const double Thresh,         // in: forward step threshold
    const int nMinSpan,          // in: set to non zero to override internal calculation
    const int nEndSpan,          // in: set to non zero to override internal calculation
    const bool Prune,            // in: do backward pass
    const int nFastK,            // in: Fast MARS K
    const double FastBeta,       // in: Fast MARS ageing coef
    const double NewVarPenalty,  // in: penalty for adding a new variable
    const int LinPreds[],        // in: nPreds x 1, 1 if predictor must enter linearly
    const double AdjustEndSpan,  // in:
    const bool AutoLinPreds,     // in: assume predictor linear if knot is min predictor value
    const bool UseBetaCache,     // in: 1 to use the beta cache, for speed
    const double Trace,          // in: 0 none 1 overview 2 forward 3 pruning 4 more pruning
    const char* sPredNames[])    // in: predictor names in trace printfs, can be NULL
{
#if _MSC_VER && _DEBUG
    InitMallocTracking();
#endif
    TraceGlobal = Trace;
    nMinSpanGlobal = nMinSpan;
    nEndSpanGlobal = nEndSpan;
    AdjustEndSpanGlobal = AdjustEndSpan;

    // nUses is the number of time each predictor is used in the model
    nUses = (int*)malloc1(nPreds * sizeof(int),
                        "nUses\t\t\tnPreds %d sizeof(int) %d",
                        nPreds, sizeof(int));

    // nDegree is degree of each term, degree of intercept is considered to be 0
    nDegree = (int*)malloc1(nMaxTerms * sizeof(int),
                            "nDegree\t\tnMaxTerms %d sizeof(int) %d",
                            nMaxTerms, sizeof(int));

    double* yw = NULL;
#if WEIGHTS
    if(WeightsArg) {
        error("weights are not yet supported in STANDALONE earth"); // TODO
        yw = (double*)malloc1(nCases * nResp * sizeof(double),
                        "yw\t\t\tnCases %d nResp %d sizeof(double) %d",
                        (const int)nCases, nResp, sizeof(double));
        for(int iResp = 0; iResp < nResp; iResp++)
            for(int i = 0; i < (const int)nCases; i++) {
                const int j = iResp * (const int)nCases + i;
                yw[j] = sqrt(WeightsArg[i]) * y[j];
            }
    }
#else
    ASSERT(WeightsArg == NULL); // weights are not currently supported
#endif

    int nTerms = 0, iTermCond = 0;
    ForwardPass(&nTerms, &iTermCond,
        BestSet, bx, Dirs, Cuts, nDegree, nUses,
        x, y, yw, WeightsArg,
        nCases, nResp, nPreds, nMaxDegree, nMaxTerms,
        Penalty, Thresh, nFastK, FastBeta, NewVarPenalty,
        LinPreds, AdjustEndSpan, AutoLinPreds, UseBetaCache, sPredNames);

    // ensure bx is full rank by updating BestSet, and get Residuals and Betas

    RegressAndFix(Betas, Residuals, NULL, BestSet,
        bx, yw? yw: y, nCases, nResp, nMaxTerms);

    if(TraceGlobal >= 6)
        PrintSummary(nMaxTerms, nTerms, nPreds, nResp,
            BestSet, Dirs, Cuts, Betas, nDegree);

    int nMaxTerms1 = DiscardUnusedTerms(bx, Dirs, Cuts, BestSet, nDegree,
                        nMaxTerms, nPreds, nCases);
    if(Prune)
        BackwardPass(pBestGcv, BestSet, Residuals, Betas,
            bx, yw? yw: y, WeightsArg, nCases, nResp, nMaxTerms1, Penalty);
    else if(WeightsArg) {
        // TODO should use weighted regression in RegressAndFix
        // UnweightBx(bx, WeightsArg, nCases, nMaxTerms);
        RegressAndFix(Betas, Residuals, NULL, BestSet,
            bx, y, nCases, nResp, nMaxTerms);
    }
    if(TraceGlobal >= 6)
        PrintSummary(nMaxTerms, nMaxTerms1, nPreds, nResp,
            BestSet, Dirs, Cuts, Betas, nDegree);

    *pnTerms = nMaxTerms1;
    *piTermCond = iTermCond;
    if(yw)
        free1(yw);
    free1(nDegree);
    free1(nUses);
}
#endif // STANDALONE

//-----------------------------------------------------------------------------
// Return the max number of knots in any term.
// Lin dep factors are considered as having one knot (at the min value of the predictor)

#if STANDALONE
static int GetMaxKnotsPerTerm(
    const bool   UsedCols[],    // in
    const int    Dirs[],        // in
    const int    nPreds,        // in
    const int    nTerms,        // in
    const int    nMaxTerms)     // in
{
    int nKnotsMax = 0;
    for(int iTerm = 1; iTerm < nTerms; iTerm++)
        if(UsedCols[iTerm]) {
            int nKnots = 0; // number of knots in this term
            for(int iPred = 0; iPred < nPreds; iPred++)
                if(Dirs_(iTerm, iPred) != 0)
                    nKnots++;
            if(nKnots > nKnotsMax)
                nKnotsMax = nKnots;
        }
    return nKnotsMax;
}
#endif // STANDALONE

//-----------------------------------------------------------------------------
// print a string representing the earth expression, one term per line
// TODO spacing is not quite right and is overly complicated

#if STANDALONE
static void FormatOneResponse(
    const bool   UsedCols[],// in: nMaxTerms x 1, indices of best set of cols of bx
    const int    Dirs[],    // in: nMaxTerms x nPreds, -1,0,1,2 for iTerm, iPred
    const double Cuts[],    // in: nMaxTerms x nPreds, cut for iTerm, iPred
    const double Betas[],   // in: nMaxTerms x nResp
    const int    nPreds,
    const int    iResp,
    const int    nTerms,
    const int    nMaxTerms,
    const int    nDigits,   // number of significant digits to print
    const double MinBeta)   // terms with fabs(beta) less than this are not printed, 0 for all
{
    int iBestTerm = 0;
    int nKnotsMax = GetMaxKnotsPerTerm(UsedCols, Dirs, nPreds, nTerms, nMaxTerms);
    int nKnots = 0;
    char s[1000];
    ASSERT(nDigits >= 0);
    char sFormat[50];  sprintf(sFormat,  "%%-%d.%dg", nDigits+6, nDigits);
    char sFormat1[50]; sprintf(sFormat1, "%%%d.%dg",  nDigits+6, nDigits);
    int nPredWidth;
    if(nPreds > 100)
        nPredWidth = 3;
    else if(nPreds > 10)
        nPredWidth = 2;
    else
        nPredWidth = 1;
    char sPredFormat[20]; sprintf(sPredFormat, "%%%dd", nPredWidth);
    char sPad[500]; sprintf(sPad, "%*s", 28+nDigits+nPredWidth, " ");    // comment pad
    const int nUsedCols = nTerms;       // nUsedCols is needed for the Betas_ macro
    printf(sFormat, Betas_(0, iResp));  // intercept
    while(nKnots++ < nKnotsMax)
        printf("%s", sPad);
    printf(" // 0\n");

    for(int iTerm = 1; iTerm < nTerms; iTerm++)
        if(UsedCols[iTerm]) {
            iBestTerm++;
            if(fabs(Betas_(iBestTerm, iResp)) >= MinBeta) {
                printf("%+-9.3g", Betas_(iBestTerm, iResp));
                nKnots = 0;
                for(int iPred = 0; iPred < nPreds; iPred++) {
                    switch(Dirs_(iTerm, iPred)) {
                        case  0:
                            break;
                        case -1:
                            sprintf(s, " * max(0, %s - %*sx[%s])",
                                sFormat, nDigits+2, " ", sPredFormat);
                            printf(s, Cuts_(iTerm, iPred), iPred);
                            nKnots++;
                            break;
                        case  1:
                            sprintf(s, " * max(0, x[%s]%*s-  %s)",
                                sPredFormat,  nDigits+2, " ", sFormat1);
                            printf(s, iPred, Cuts_(iTerm, iPred));
                            nKnots++;
                            break;
                        case  2:
                            sprintf(s, " * x[%s]%*s                    ",
                                sPredFormat,  nDigits+2, " ");
                            printf(s, iPred);
                            nKnots++;
                            break;
                        default:
                            ASSERT(false);
                            break;
                    }
                }
                while(nKnots++ < nKnotsMax)
                    printf("%s", sPad);
                printf(" // %d\n", iBestTerm);
            }
        }
}

void FormatEarth(
    const bool   UsedCols[],// in: nMaxTerms x 1, indices of best set of cols of bx
    const int    Dirs[],    // in: nMaxTerms x nPreds, -1,0,1,2 for iTerm, iPred
    const double Cuts[],    // in: nMaxTerms x nPreds, cut for iTerm, iPred
    const double Betas[],   // in: nMaxTerms x nResp
    const int    nPreds,
    const int    nResp,     // in: number of cols in y
    const int    nTerms,
    const int    nMaxTerms,
    const int    nDigits,   // number of significant digits to print
    const double MinBeta)   // terms with fabs(betas) less than this are not printed, 0 for all
{
    for(int iResp = 0; iResp < nResp; iResp++) {
        if(nResp > 1)
            printf("Response %d:\n", iResp+IOFFSET);
        FormatOneResponse(UsedCols, Dirs, Cuts, Betas, nPreds, iResp,
            nTerms, nMaxTerms, nDigits, MinBeta);
    }
}
#endif // STANDALONE

//-----------------------------------------------------------------------------
// return the value predicted by an earth model, given  a vector of inputs x

#if STANDALONE
static double PredictOneResponse(
    const double x[],        // in: vector nPreds x 1 of input values
    const bool   UsedCols[], // in: nMaxTerms x 1, indices of best set of cols of bx
    const int    Dirs[],     // in: nMaxTerms x nPreds, -1,0,1,2 for iTerm, iPred
    const double Cuts[],     // in: nMaxTerms x nPreds, cut for iTerm, iPred
    const double Betas[],    // in: nMaxTerms x 1
    const int    nPreds,     // in: number of cols in x
    const int    nTerms,
    const int    nMaxTerms)
{
    double yHat = Betas[0];
    int iTerm1 = 0;
    for(int iTerm = 1; iTerm < nTerms; iTerm++)
        if(UsedCols[iTerm]) {
            iTerm1++;
            double Term = Betas[iTerm1];
            for(int iPred = 0; iPred < nPreds; iPred++) {
                switch(Dirs_(iTerm, iPred)) {
                    case  0: break;
                    case -1: Term *= max(0, Cuts_(iTerm, iPred) - x[iPred]); break;
                    case  1: Term *= max(0, x[iPred] - Cuts_(iTerm, iPred)); break;
                    case  2: Term *= x[iPred]; break;
                    default: ASSERT("bad direction" == NULL); break;
                }
            }
            yHat += Term;
        }
    return yHat;
}

void PredictEarth(
    double       y[],        // out: vector nResp
    const double x[],        // in: vector nPreds x 1 of input values
    const bool   UsedCols[], // in: nMaxTerms x 1, indices of best set of cols of bx
    const int    Dirs[],     // in: nMaxTerms x nPreds, -1,0,1,2 for iTerm, iPred
    const double Cuts[],     // in: nMaxTerms x nPreds, cut for iTerm, iPred
    const double Betas[],    // in: nMaxTerms x nResp
    const int    nPreds,     // in: number of cols in x
    const int    nResp,      // in: number of cols in y
    const int    nTerms,
    const int    nMaxTerms)
{
    for(int iResp = 0; iResp < nResp; iResp++)
        y[iResp] = PredictOneResponse(x, UsedCols, Dirs, Cuts,
                       Betas + iResp * nTerms, nPreds, nTerms, nMaxTerms);
}
#endif // STANDALONE

//-----------------------------------------------------------------------------
// Example main routine
// See earth/inst/slowtests/test.earthc.c for more complex examples

#if STANDALONE && MAIN
void error(const char *args, ...)       // params like printf
{
    char s[1000];
    va_list va;
    va_start(va, args);
    vsprintf(s, args, va);
    va_end(va);
    printf("\nError: %s\n", s);
    exit(-1);
}

// extern here prevents clang -Wmissing-prototypes warning
extern void xerbla_(char *srname, int* info);

void xerbla_(char *srname, int* info)   // needed by BLAS and LAPACK routines
{
    char buf[7];
    strncpy(buf, srname, 6);
    buf[6] = 0;
    error("BLAS/LAPACK routine %6s gave error code %d", buf, -(*info));
}

int main(void)
{
    const int    nMaxTerms = 21;  // called "nk" in the R code
    const size_t nCases = 100;    // note that nCases is size_t, not int
                                  // this allows e.g. mallocs below to be bigger than 2GB
    const int    nResp = 1;       // number of responses i.e. number of y columns
    const int    nPreds = 1;      // number of predictors i.e. number of x columns
    const int    nMaxDegree = 1;  // called "degree" in the R code
    const double Penalty = (nMaxDegree > 1)? 3: 2;
    const double Thresh = .001;
    const int    nMinSpan = 0;    // 0 means auto
    const int    nEndSpan = 0;    // 0 means auto
    const bool   Prune = true;
    const int    nFastK = 20;
    const double FastBeta = 1;
    const double NewVarPenalty = 0;
    int*         LinPreds = (int*)calloc1(nPreds, sizeof(int), NULL); // "linpreds" in R code
    const double AdjustEndSpan = 2.0;
    const bool   AutoLinPreds = true;
    const bool   UseBetaCache = true;
    const double Trace = 3;
    const char** sPredNames = NULL;

    double  BestGcv;
    int     nTerms;
    int     iTermCond;
    bool*   BestSet   = (bool*)  malloc1(nMaxTerms * sizeof(bool), NULL);
    double* bx        = (double*)malloc1(nCases    * nMaxTerms * sizeof(double), NULL);
    int*    Dirs      = (int*)   malloc1(nMaxTerms * nPreds    * sizeof(int),    NULL);
    double* Cuts      = (double*)malloc1(nMaxTerms * nPreds    * sizeof(double), NULL);
    double* Residuals = (double*)malloc1(nCases    * nResp     * sizeof(double), NULL);
    double* Betas     = (double*)malloc1(nMaxTerms * nResp     * sizeof(double), NULL);

    double* x = (double*)malloc1(nCases * nPreds   * sizeof(double), NULL);
    double* y = (double*)malloc1(nCases * nResp * sizeof(double),    NULL);

    ASSERT(nResp == 1);         // code below only works for nResp == 1

    for(int i = 0; i < (const int)nCases; i++) {
        const double xi = i / (double)nCases;
        x[i] = xi;
        y[i] = sin(4 * xi);     // target function, change this to whatever you want
    }
    Earth(&BestGcv, &nTerms, &iTermCond, BestSet, bx, Dirs, Cuts, Residuals, Betas,
        x, y, NULL /*WeightsArg*/, nCases, nResp, nPreds,
        nMaxDegree, nMaxTerms, Penalty, Thresh, nMinSpan, nEndSpan, Prune,
        nFastK, FastBeta, NewVarPenalty, LinPreds,
        AdjustEndSpan, AutoLinPreds, UseBetaCache, Trace, sPredNames);

    printf("Expression:\n");
    FormatEarth(BestSet, Dirs, Cuts, Betas, nPreds, nResp, nTerms, nMaxTerms, 3, 0);

    double x1 = 0.1234, y1;
    PredictEarth(&y1,
        &x1, BestSet, Dirs, Cuts, Betas, nPreds, nResp, nTerms, nMaxTerms);
    printf("\nf(%g) = %g\n", x1, y1);

    free1(y);
    free1(x);
    free1(Betas);
    free1(Residuals);
    free1(Cuts);
    free1(Dirs);
    free1(bx);
    free1(BestSet);
    free1(LinPreds);

    return 0;
}
#endif // STANDALONE && MAIN
