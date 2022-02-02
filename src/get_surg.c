// This routine is originally from rpart.
//
// The main entry point for R to find surrogate variables for a node.
//
// Input variables:
//      wt      = vector of case weights
//      xmat    = matrix of continuous variables
//      opt     = vector of options.  Same order as get_surg.control, as a vector
//                of doubles.
//      node    = primary node
//
// Returned: a list with elements
//      dsplit = for each split, numeric variables (doubles)
//      isplit = for each split, integer variables
//      isplitcntSum =  for each node, integer variables
//
// Naming convention: ... = pointer to an integer vector, ...R = the
// input R object (SEXP) containing that vector

#define MAINRP
#include <math.h>
#include "get_surg.h"
#include "node.h"
#include "rpartproto.h"

#ifdef DEBUG
#include <unistd.h>	//for using the function sleep
#endif

SEXP getSurrogates(SEXP wt, SEXP testpar, SEXP xmat, SEXP opt, SEXP node) {

#ifdef DEBUG
	Rprintf("debug\n");
#endif
	double *dptr;
	int *iptr;
	int i, j, k, n;
	pNode tree;

	// return objects for R
	// end in "3" to avoid overlap with internal names

	SEXP dsplit3, isplit3;

	// output variables
	int nsplit[2];
	double *ddsplit[3];
	int *iisplit[3];

	// hand over arguments
	rp.n = nrows(xmat);
	n = rp.n;
	rp.nvar = ncols(xmat);

	rp.wt = REAL(wt);

	iptr = INTEGER(opt);
	rp.maxsur = (int) iptr[0];
	rp.sur_agree = (int) iptr[1];

#ifdef  OPENMP_ON
	rp.nthreads = (int) iptr[2];
	if (rp.nthreads < 1 || omp_get_num_procs() < rp.nthreads)
		rp.nthreads = omp_get_num_procs();
#endif
	// create pointers to the matrix
	// x and missmat are in column major order
	// y is in row major order
	dptr = REAL(xmat);
	rp.xdata = (double **) calloc(rp.nvar, sizeof(double *));
#ifdef OPENMP_ON
	#pragma omp parallel for private(i) schedule(dynamic) num_threads(rp.nthreads)
	// start parallel section
#endif
	for (i = 0; i < rp.nvar; i++) {
		//combined for parallel execution
		rp.xdata[i] = dptr + i * n;
	}
	// end parallel section

	// create matrix of sort indices, for surrogate
	// one for each continuous variable
	// This sort is "once and for all".
	// I don't have to sort the categories.

	rp.sorts = (int **) calloc(rp.nvar, sizeof(int *));
	rp.sorts[0] = (int *) calloc(n * rp.nvar, sizeof(int));

#ifdef OPENMP_ON
	#pragma omp parallel private(i,k) num_threads(rp.nthreads)
	{
// start parallel section
#endif
		int *tempvec = (int *) calloc(n, sizeof(int));
		double *xtemp = (double *) calloc(n, sizeof(double));
#ifdef OPENMP_ON
		#pragma omp for schedule(dynamic)
#endif
		for (i = 0; i < rp.nvar; i++) {
			rp.sorts[i] = rp.sorts[0] + i * n;
			for (k = 0; k < n; k++) {
				if (!R_FINITE(rp.xdata[i][k])) {
					// this variable is missing (NA)
					tempvec[k] = -(k + 1);
					xtemp[k] = 0;
				} else {
					tempvec[k] = k;
					xtemp[k] = rp.xdata[i][k];
				}
			}
			sort_vec(0, n - 1, xtemp, tempvec);
			for (k = 0; k < n; k++)
				rp.sorts[i][k] = tempvec[k];
		}
		//clear DMA
		free(tempvec);
		free(xtemp);
#ifdef OPENMP_ON
	}
// end parallel section
#endif
	// finalize tree object
	nodesize = sizeof(Node);
	tree = (pNode) calloc(1, nodesize);

	dptr = REAL(node);
	// the split structure is sized for 2 categories.
	int splitsize = sizeof(Split);
	tree->primary = (pSplit) calloc(1, splitsize);
	tree->primary->csplit = 1;
	tree->primary->nextsplit = NULL;

	tree->primary->var_num = (int) dptr[0] - 1;
	tree->primary->spoint = dptr[1];

	if (0 < rp.maxsur)
		surrogate(tree, 0, n);


	// Return the body of the tree
	// For each component we first create a vector to hold the
	// result, then a ragged array index into the vector.
	// The rpmatrix routine then fills everything in.

	rpcountup(tree, nsplit);
	int splitcnt = nsplit[0] + nsplit[1];

	dsplit3 = PROTECT(allocMatrix(REALSXP, splitcnt, 3));
	dptr = REAL(dsplit3);
#ifdef OPENMP_ON
	#pragma omp parallel for private(i,j) schedule(dynamic) collapse(2) num_threads(rp.nthreads)
	// start parallel section
#endif
	for (i = 0; i < 3; i++) {
		for (j = 0; j < splitcnt; j++) {
			//moved for perfectly nested loop | combined for parallel execution
			ddsplit[i] = dptr + i * splitcnt;
			ddsplit[i][j] = 0.0;
		}
	}
	// end parallel section

	isplit3 = PROTECT(allocMatrix(INTSXP, splitcnt, 3));
	iptr = INTEGER(isplit3);
// #ifdef OPENMP_ON
// 	#pragma omp parallel for private(i) schedule(dynamic) num_threads(4)
// 	// removed too short parallel section
// #endif
	for (i = 0; i < 3; i++) {
		//combined for parallel execution
		iisplit[i] = iptr + i * splitcnt;
	}

	rpmatrix(tree, ddsplit, iisplit, nsplit);

	// Create the output list
	int nout = 2;
	SEXP rlist = PROTECT(allocVector(VECSXP, nout));
	SEXP rname = allocVector(STRSXP, nout);
	setAttrib(rlist, R_NamesSymbol, rname);
	SET_VECTOR_ELT(rlist, 0, dsplit3);
	SET_STRING_ELT(rname, 0, mkChar("dsplit"));
	SET_VECTOR_ELT(rlist, 1, isplit3);
	SET_STRING_ELT(rname, 1, mkChar("isplit"));

	UNPROTECT(1 + nout);

	// let the memory go
	free_tree(tree, 0);

	// clear DMA
	free(rp.xdata);
	free(rp.sorts[0]);
	free(rp.sorts);
	free(tree);

	return rlist;
}