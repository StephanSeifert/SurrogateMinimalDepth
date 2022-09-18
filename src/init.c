// This routine export the c function for R
#include "get_surg.h"
#include "R_ext/Rdynload.h"
#include "node.h"
#include "rpartproto.h"

SEXP getSurrogates(SEXP ncat2, SEXP wt, SEXP xmat, SEXP opt, SEXP var, SEXP split);

static const R_CallMethodDef CallEntries[] = {
	// registering native routines http://www.hep.by/gnu/r-patched/r-exts/R-exts_95.html
	{ "getSurrogates", (DL_FUNC) & getSurrogates, 6 }, { NULL, NULL, 0 }
};

#include <Rversion.h>
void
// registering native routines http://www.hep.by/gnu/r-patched/r-exts/R-exts_95.html
R_init_SurrogateMinimalDepth(DllInfo * dll) {
	R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
#if defined(R_VERSION) && R_VERSION >= R_Version(2, 16, 0)
	R_forceSymbols(dll, TRUE);
#endif
}
