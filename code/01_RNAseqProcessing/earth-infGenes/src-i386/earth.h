// earth.h: externs for earth.c
//
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

#if !defined(EARTH_H)
#define EARTH_H

#if USING_R

void FreeR(void);

void ForwardPassR(              // for use by R
    int    FullSet[],           // out: nMaxTerms x 1, bool vec of lin indep cols of bx
    double bx[],                // out: MARS basis matrix, nCases x nMaxTerms
    double Dirs[],              // out: nMaxTerms x nPreds, elements are -1,0,1,2
    double Cuts[],              // out: nMaxTerms x nPreds, cut for iTerm,iPred
    int*   piTermCond,          // out: reason we terminated the forward pass
    const double x[],           // in: nCases x nPreds, unweighted x
    const double y[],           // in: nCases x nResp, unweighted but scaled y
    const double yw[],          // in: nCases x nResp, weighted and scaled y
    const double WeightsArg[],  // in: nCases x 1, never MyNull
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
    const SEXP Allowed,            // in: constraints function, can be MyNull
    const int* pnAllowedFuncArgs,  // in: number of arguments to Allowed function, 3 or 4
    const SEXP Env,                // in: environment for Allowed function
    const double* pAdjustEndSpan,  // in:
    const int* pnAutoLinPreds,      // in: assume predictor linear if knot is min predictor value
    const int* pnUseBetaCache,     // in: 1 to use the beta cache, for speed
    const double* pTrace,          // in: 0 none 1 overview 2 forward 3 pruning 4 more pruning
    const char* sPredNames[],      // in: predictor names in trace printfs, can be MyNull
    const double* MyNullDouble,    // in: trick to avoid R check warnings on passing R_NilValue
    const SEXP MyNullFunc);        // in: ditto

void EvalSubsetsUsingXtxR(      // for use by R
    double        PruneTerms[], // out: specifies which cols in bx are in best set
    double        RssVec[],     // out: nTerms x 1
    const int*    pnCases,      // in
    const int*    pnResp,       // in: number of cols in y
    const int*    pnMaxTerms,   // in
    const double  bx[],         // in: MARS basis matrix, all cols must be indep
    const double  y[],          // in: nCases * nResp (possibly weighted)
    const double* pTrace);      // in

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
    const bool   UsedCols[]);   // in: specifies used columns in x

#endif // USING_R

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
    const double AdjustEndSpan,  // in: for adjusting endspan for interaction terms
    const bool AutoLinPreds,      // in: assume predictor linear if knot is max predictor value
    const bool UseBetaCache,     // in: 1 to use the beta cache, for speed
    const double Trace,          // in: 0 none 1 overview 2 forward 3 pruning 4 more pruning
    const char* sPredNames[]);   // in: predictor names in trace printfs, can be NULL

void FormatEarth(
    const bool   UsedCols[], // in: nMaxTerms x 1, indices of best set of cols of bx
    const int    Dirs[],     // in: nMaxTerms x nPreds, -1,0,1,2 for iTerm, iPred
    const double Cuts[],     // in: nMaxTerms x nPreds, cut for iTerm, iPred
    const double Betas[],    // in: nMaxTerms x nResp
    const int    nPreds,
    const int    nResp,      // in: number of cols in y
    const int    nTerms,
    const int    nMaxTerms,
    const int    nDigits,    // number of significant digits to print
    const double MinBeta);   // terms with fabs(betas) less than this are not printed, 0 for all

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
    const int    nMaxTerms);

#endif // STANDALONE

#endif // EARTH_H
