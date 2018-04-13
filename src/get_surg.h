/*
 * This routine is orignaly from rpart.
 *
 * commom variables for the rpart routine
 *
 * Start with things that depend on R.h
 */

#include <R.h>
#include <Rinternals.h>

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("SurrogateMinimalDepth", String)
#else
#define _(String) (String)
#endif

/*
 * Memory defined with R_alloc is removed automatically
 * That with "CALLOC" I have to remove myself. Use the
 * latter for objects that need to persist between the
 * s_to_rp1 and s_to_rp2 calls
 */
#define ALLOC(a,b) R_alloc(a,b)
#define CALLOC(a,b) R_chk_calloc((size_t)(a), b)

/* done with the R internals */
/* used for the variable "extra" in nodes */
#define LEFT  (-1)              
#define RIGHT  1

#ifdef MAINRP
#define EXTERN
#else
#define EXTERN extern
#endif

/* 
 * As a sop to S, I need to keep the total number of external symbols
 * somewhat smaller. So, pack most of them all into a structure.
 */
EXTERN struct {
		double **xdata;
		double *xtemp;
		double *wt;
		double *lwt;
		/* scratch double vectors, of length ncat */
		double *rwt;
		/* variable type: 0=cont, 1+ =#categories */
		int *numcat;
		/* matrix of sort indices */
		int **sorts;
		/* total number of subjects */
		int n;
		/* number of predictors */
		int nvar;
		/* max # of primary or surrogate splits to use */
		int maxsur;
		/* 0 = my style, 1=CART style */
		int sur_agree;
		/* to be allocated by the mainline, of length n */
		int *tempvec;
		int *csplit;
		int *left;
		int *right;
	} rp;

	EXTERN int nodesize;

	/*
	 * Categorical variables must be coded as 1,2,3, ..., and there may be
	 * missing categories. The upper limit is determined on the fly.
	 */
