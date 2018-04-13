/*
 * This routine is orignaly from rpart.
 *
 * The main entry point for R to find surrogate variables for a node.
 *
 * Input variables:
 *      ncat    = # categories for each var, 0 for continuous variables.
 *      wt      = vector of case weights
 *      xmat    = matrix of continuous variables
 *      opt     = vector of options.  Same order as get_surg.control, as a vector
 *                of doubles.
 *      node    = primary node
 *
 * Returned: a list with elements
 *      dsplit = for each split, numeric variables (doubles)
 *      isplit = for each split, integer variables
 *      isplitcount =  for each node, integer variables
 *
 * Naming convention: ncat = pointer to an integer vector, ncatR = the
 * input R object (SEXP) containing that vector
 */
#define MAINRP
#include <math.h>
#include "get_surg.h"
#include "node.h"
#include "rpartproto.h"

SEXP getSurrogates(SEXP ncatR, SEXP wt, SEXP xmat, SEXP opt, SEXP node) {
	int *ncat;
	double *dptr;
	int *iptr;
	int i, j, k, n;
	int maxcat;
	pNode tree;

	/* return objects for R 
	 * end in "3" to avoid overlap with internal names 
	 */
	SEXP dsplit3, isplit3;

	/* output variables */
	int splitcount;
	double *ddsplit[3];
	int *iisplit[3];

	/* hand over arguments */
	splitcount = 0;

	ncat = INTEGER(ncatR);
	rp.numcat = ncat;

	rp.n = nrows(xmat);
	n = rp.n;
	rp.nvar = ncols(xmat);

	rp.wt = REAL(wt);

	/*
	 * create pointers to the matrix
	 * x and missmat are in column major order
	 * y is in row major order
	 */
	dptr = REAL(xmat);
	rp.xdata = (double **) ALLOC(rp.nvar, sizeof(double *));
	for (i = 0; i < rp.nvar; i++) {
		rp.xdata[i] = dptr;
		dptr += n;
	}

	iptr = INTEGER(opt);
	rp.maxsur = (int) iptr[4];

	/*
	 * create matrix of sort indices, for surrogate
	 * one for each continuous variable
	 * This sort is "once and for all".
	 * I don't have to sort the categoricals.
	 */
	rp.tempvec = (int *) ALLOC(n, sizeof(int));
	rp.xtemp = (double *) ALLOC(n, sizeof(double));
	rp.sorts = (int **) ALLOC(rp.nvar, sizeof(int *));
	rp.sorts[0] = (int *) ALLOC(n * rp.nvar, sizeof(int));
	maxcat = 0;
	for (i = 0; i < rp.nvar; i++) {
		rp.sorts[i] = rp.sorts[0] + i * n;
		for (k = 0; k < n; k++) {
			if (!R_FINITE(rp.xdata[i][k])) {
				/* this variable is missing (NA)*/
				rp.tempvec[k] = -(k + 1);
				rp.xtemp[k] = 0;
			} else {
				rp.tempvec[k] = k;
				rp.xtemp[k] = rp.xdata[i][k];
			}
		}
		if (ncat[i] == 0)
			sort_vec(0, n - 1, rp.xtemp, rp.tempvec);
		else
			if (ncat[i] > maxcat)
				maxcat = ncat[i];
		for (k = 0; k < n; k++)
			rp.sorts[i][k] = rp.tempvec[k];
	}

	/* finalize rp object */
	if (maxcat > 0) {
		rp.csplit = (int *) ALLOC(3 * maxcat, sizeof(int));
		rp.lwt = (double *) ALLOC(2 * maxcat, sizeof(double));
		rp.left = rp.csplit + maxcat;
		rp.right = rp.left + maxcat;
		rp.rwt = rp.lwt + maxcat;
	} else
		rp.csplit = (int *) ALLOC(1, sizeof(int));

	/* finalize tree object */
	nodesize = sizeof(Node);
	tree = (pNode) ALLOC(1, nodesize);
	memset(tree, 0, nodesize);

	dptr = REAL(node);
	/* the split structure is sized for 2 categories. */
	int splitsize = sizeof(Split);
	tree->primary = (pSplit) CALLOC(1, splitsize);
	memset(tree->primary, 0, splitsize);
	tree->primary->csplit[0] = 1;
	tree->primary->nextsplit = NULL;

	tree->primary->var_num = (int) dptr[0] - 1;
	tree->primary->spoint = dptr[1];

	if (rp.maxsur > 0)
		surrogate(tree, 0, n);

	/*
	 * Return the body of the tree
	 * For each component we first create a vector to hold the
	 * result, then a ragged array index into the vector.
	 * The rpmatrix routine then fills everything in.
	 */
	rpcountup(tree, &splitcount);

	dsplit3 = PROTECT(allocMatrix(REALSXP, splitcount, 3));
	dptr = REAL(dsplit3);
	for (i = 0; i < 3; i++) {
		ddsplit[i] = dptr;
		dptr += splitcount;
		for (j = 0; j < splitcount; j++)
			ddsplit[i][j] = 0.0;
	}

	isplit3 = PROTECT(allocMatrix(INTSXP, splitcount, 3));
	iptr = INTEGER(isplit3);
	for (i = 0; i < 3; i++) {
		iisplit[i] = iptr;
		iptr += splitcount;
	}

	rpmatrix(tree, rp.numcat, ddsplit, iisplit, NULL, 1);

	/* Create the output list */
	int nout = 2;
	SEXP rlist = PROTECT(allocVector(VECSXP, nout));
	SEXP rname = allocVector(STRSXP, nout);
	setAttrib(rlist, R_NamesSymbol, rname);
	SET_VECTOR_ELT(rlist, 0, dsplit3);
	SET_STRING_ELT(rname, 0, mkChar("dsplit"));
	SET_VECTOR_ELT(rlist, 1, isplit3);
	SET_STRING_ELT(rname, 1, mkChar("isplit"));

	UNPROTECT(1 + nout);

	/* let the memory go */
	free_tree(tree, 0);

	return rlist;
}
