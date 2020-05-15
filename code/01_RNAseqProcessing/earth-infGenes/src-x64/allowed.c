// allowed.c: routines for the "allowed" parameter of the R function earth().

#include "R.h"
#include "Rinternals.h"
#define printf Rprintf

#ifndef _MSC_VER // microsoft
#ifndef bool
    typedef int bool;
    #define false 0
    #define true  1
#endif
#endif

#define Dirs_(iTerm,iPred) Dirs[(iTerm) + (iPred)*(nMaxTerms)]

static SEXP AllowedFunc;
static SEXP AllowedEnv;
static int  nArgs;
static bool First;

// Initialize the R function AllowedFunc from the Allowed function
// argument which was passed into ForwardPassR.
// For efficiency, we initialize once here rather than in IsAllowed.
//
// The caller of ForwardPassR has already checked that Allowed is
// a function and has three arguments: degree, pred, parents.
//
// The "allowed" function has the following prototype, where
// namesx and first are optional.
//
//     allowed <- function(degree, pred, parents, namesx, first)
//     {
//         ...
//         TRUE   # return TRUE if allowed
//     }
//
// where "degree" is the MARS term degree, with pred in the term.
//       "pred" is column index in the input matrix x
//       "parents" is an integer vector of parent predictors
//                 (it's a copy of Dirs[iParent,]
//       "namesx" is optional and is the colnames of the x arg
//                to earth, after factor expansion
//       "first" is optional and is 1 the first time "allowed"
//               is invoked for the current model

void InitAllowedFunc(
        SEXP Allowed, // can be NULL
        int nAllowedFuncArgs, SEXP Env,
        const char* sPredNames[], int nPreds)
{
    if(Allowed == NULL)
        AllowedFunc = NULL;
    else {
        if(nAllowedFuncArgs < 3 || nAllowedFuncArgs > 5)
            error("Bad nAllowedFuncArgs %d", nAllowedFuncArgs);

        AllowedEnv = Env;
        nArgs = nAllowedFuncArgs;

        // the UNPROTECT for the PROTECT below is in FreeAllowedFunc()
        PROTECT(AllowedFunc = allocList(1 + nAllowedFuncArgs));

        SEXP s = AllowedFunc;   // 1st element is the function
        SETCAR(s, Allowed);
        SET_TYPEOF(s, LANGSXP);

        s = CDR(s);             // 2nd element is "degree"
        SETCAR(s, allocVector(INTSXP, 1));

        s = CDR(s);             // 3rd element is "pred"
        SETCAR(s, allocVector(INTSXP, 1));

        s = CDR(s);             // 4th element is "parents"
        SETCAR(s, allocVector(INTSXP, nPreds));

        if(nAllowedFuncArgs >= 4) {
            SEXP namesx;
            s = CDR(s);        // 5th element is "namesx"
            SETCAR(s, namesx = allocVector(STRSXP, nPreds));
            PROTECT(namesx);
            if(sPredNames == NULL)
                error("Bad sPredNames");
            for(int i = 0; i < nPreds; i++)
                SET_STRING_ELT(namesx, i, mkChar(sPredNames[i]));
            UNPROTECT(1);
        }
        if(nAllowedFuncArgs >= 5) {
            s = CDR(s);        // 6th element is "first"
            SETCAR(s, allocVector(LGLSXP, 1));
        }
    }
    First = true;
}

void FreeAllowedFunc(void)
{
    if(AllowedFunc != NULL) {
        UNPROTECT(1);          // matches PROTECT in InitAllowedFunc
        AllowedFunc = NULL;
    }
}

static bool EvalAllowedFunc(void)
{
    if(AllowedFunc == NULL)
        error("EvalAllowedFunc: AllowedFunc == NULL");

    SEXP s = eval(AllowedFunc, AllowedEnv);

    bool allowed;
    switch(TYPEOF(s)) {         // be fairly permissive with return type
        case LGLSXP:
            allowed = (bool)(LOGICAL(s)[0] != 0);
            break;
        case INTSXP:
            allowed = INTEGER(s)[0] != 0;
            break;
        case REALSXP:
            allowed = (bool)(REAL(s)[0] != 0.);
            break;
        default:
            error("the \"allowed\" function returned a %s instead of a logical",
                  Rf_type2char(TYPEOF(s)));
            allowed = FALSE; // -Wall
            break;
    }
    if(LENGTH(s) != 1)
        error("the \"allowed\" function did not return a logical of length 1");

    return allowed;
}

// Return TRUE if the current iPred can be used in a term with iParent
// i.e. TRUE means no constraint.
//
// This calls the R function Allowed which was passed in as a parameter to
// ForwardPassR.  The fields of Allowed have been preallocated into
// AllowedFunc and so all we do here is fill in the values and call eval.

bool IsAllowed(
    const int iPred,        // in: candidate predictor
    const int iParent,      // in: candidate parent term
    const int Dirs[],       // in:
    const int nPreds,       // in:
    const int nMaxTerms)    // in:
{
    if(AllowedFunc == NULL)
       return TRUE;

    SEXP s = AllowedFunc;           // 1st element is the function
    s = CDR(s);                     // 2nd element is "degree"
    INTEGER(CADR(s))[0] = iPred+1;  // 3rd element is "pred"
    int* p = INTEGER(CADDR(s));     // 4th element is "parents"
    int i, nDegree = 1;
    for(i = 0; i < nPreds; i++) {
        p[i] = Dirs_(iParent, i);
        if(p[i])
            nDegree++;
    }
    INTEGER(CAR(s))[0] = nDegree;

    // optional 5th element already initialized to predictor names

    if(nArgs >= 5)                  // optional 6th element is "first"
        *(LOGICAL(CAD4R(s))) = First;
    First = false;

    return EvalAllowedFunc();
}
