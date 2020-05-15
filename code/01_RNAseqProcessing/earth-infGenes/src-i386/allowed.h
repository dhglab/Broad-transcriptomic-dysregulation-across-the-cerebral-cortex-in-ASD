// allowed.h: externs for allowed.c

void InitAllowedFunc(
        SEXP Allowed, // can be NULL
        int nAllowedFuncArgs, SEXP Env,
        const char* sPredNames[], int nPreds);

void FreeAllowedFunc(void);

bool IsAllowed(
    const int iPred,        // in: candidate predictor
    const int iParent,      // in: candidate parent term
    const int Dirs[],       // in:
    const int nPreds,       // in:
    const int nMaxTerms);   // in:
