// This routine is originally from rpart.
//
// common variables for the rpart routine
//
// Start with things that depend on R.h


#include <R.h>
#include <Rinternals.h>

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("SurrogateMinimalDepth", String)
#else
#define _(String) (String)
#endif

// Header for OpenMP
#include <omp.h>

// Most compilers with openMP support supply a
//  pre-defined compiler macro _OPENMP. Following
//  facilitates selective turning off (by testing
//  value or defining multiple versions OPENMP_ON1,
//  OPENMP_ON2...)
//  note also that there is no actual *need* to
//  protect #pragmas with #ifdef OPENMP_ON, since C
//  ignores undefined pragmas, but failing to do so
//  may produce warnings if openMP is not supported.
//  In contrast functions from omp.h must be protected.


// TODO DEBUG
//#define DEBUG
////#ifdef DEBUG
////Rprintf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
////#endif


#ifdef DEBUG
//for using the function sleep
#include <unistd.h>
// Header for measuring execution time
//#include "stopwatch.h"
#endif

#if defined _OPENMP
#define OPENMP_ON 1
#endif


// Memory defined with R_alloc is removed automatically
// That with "CALLOC" I have to remove myself. Use the
// latter for objects that need to persist between the
// s_to_rp1 and s_to_rp2 calls

#define ALLOC(a,b) R_alloc(a,b)
#define CALLOC(a,b) R_chk_calloc((size_t)(a), b)

// done with the R internals
// used for the variable "extra" in nodes
#define LEFT  (-1)
#define RIGHT  1

#ifdef MAINRP
#define EXTERN
#else
#define EXTERN extern
#endif


// As a sop to S, I need to keep the total number of external symbols
// somewhat smaller. So, pack most of them all into a structure.

EXTERN struct {
	double **xdata;
	double *wt;
	// matrix of sort indices
	int **sorts;
	// total number of subjects
	int n;
	// number of predictors
	int nvar;
	// max # of primary or surrogate splits to use
	int maxsur;
	// 0 = my style, 1=CART style
	int sur_agree;
	// to be allocated by the mainline, of length n
	int csplit;
	int *csplit_for_cat;
	double *lwt;
	double *rwt;
	int *left;
	int *right;
	int nthreads;
	int *numcat;
} rp;

EXTERN int nodesize;


// Categorical variables must be coded as 1,2,3, ...., and there may be
// missing categories. The upper limit is determined on the fly.
