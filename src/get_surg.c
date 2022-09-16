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

SEXP getSurrogates(SEXP ncat2, SEXP wt, SEXP xmat, SEXP opt, SEXP var, SEXP split) {
	

	double *dptr;
	int *iptr;
	int i, j, k, n;
	pNode tree;

	// LCJ
	int *ncat;
	int maxcat;
	int *cat_direction;

	// return objects for R
	// end in "3" to avoid overlap with internal names

	SEXP dsplit3, isplit3;

	
	// output variables
	int nsplit[2];
	double *ddsplit[3];
	int *iisplit[3];


	//cat_direction = INTEGER(cat_direction2);

	ncat = INTEGER(ncat2);

	int test2 = 5;
	int *prim_var_num;
	prim_var_num = INTEGER(var);

	double *splitpoint;
	double *cat_dir_numeric;

	int *cat_dir_int = (int *) calloc(ncat[prim_var_num[0] - 1] + 1, sizeof(int));

	if (ncat[prim_var_num[0] - 1] == 0)
	{
		splitpoint = REAL(split);
	} 
	else
		{
			cat_dir_numeric = REAL(split);
			for (i = 0; i < cat_dir_numeric[0] + 1; i++)
			{
				cat_dir_int[i] = (int) cat_dir_numeric[i];
			}
		}

	// hand over arguments

	rp.numcat = INTEGER(ncat2);
	rp.n = nrows(xmat);
	n = rp.n;
	rp.nvar = ncols(xmat);
	rp.wt = REAL(wt);
	iptr = INTEGER(opt);
	rp.maxsur = (int) iptr[0];
	rp.sur_agree = (int) iptr[1];

	



	// create pointers to the matrix
	// x and missmat are in column major order
	// y is in row major order
	dptr = REAL(xmat);
	rp.xdata = (double **) calloc(rp.nvar, sizeof(double *));

	for (i = 0; i < rp.nvar; i++) {
		rp.xdata[i] = dptr;
		dptr += n;
	}

	// create matrix of sort indices, for surrogate
	// one for each continuous variable
	// This sort is "once and for all".
	// I don't have to sort the categories.

	rp.sorts = (int **) calloc(rp.nvar, sizeof(int *));
	rp.sorts[0] = (int *) calloc(n * rp.nvar, sizeof(int));

	maxcat = 0;


		int *tempvec = (int *) calloc(n, sizeof(int));
		double *xtemp = (double *) calloc(n, sizeof(double));

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
			if(ncat[i] == 0)
			{
				sort_vec(0, n - 1, xtemp, tempvec, i);
			} else if(ncat[i] > maxcat)
			{
				maxcat = ncat[i];
			}
			for (k = 0; k < n; k++)
			{
				rp.sorts[i][k] = tempvec[k];

				if (i == 8)
				{
					// printf("rp.sorts[i][k]: %d, xtemp[k]: %f\n", rp.sorts[i][k], xtemp[k]);
				}
			}	
			
		}
		//clear DMA
		free(tempvec);
		free(xtemp);


	if (maxcat > 0) {
		rp.csplit_for_cat = (int *) ALLOC(3 * maxcat, sizeof(int));
		rp.lwt = (double *) ALLOC(2 * maxcat, sizeof(double));
		rp.left = rp.csplit_for_cat + maxcat;
		rp.right = rp.left + maxcat;
		rp.rwt = rp.lwt + maxcat;
	} else {
		rp.csplit_for_cat = (int *) ALLOC(1, sizeof(int));
	}

	// finalize tree object
	nodesize = sizeof(Node);
	tree = (pNode) calloc(1, nodesize);

	//dptr = REAL(node);
	// the split structure is sized for 2 categories.
	int splitsize = sizeof(Split);
	tree->primary = (pSplit) calloc(1, splitsize);

	

	tree->primary->var_num = prim_var_num[0] - 1;


	if (rp.numcat[tree->primary->var_num] == 0)
	{
		tree->primary->csplit[0] = 1;
		tree->primary->spoint = splitpoint[0];
	} else 
	{
		for (i = 1; i < cat_dir_int[0]+1; i++)
		{
			tree->primary->csplit[i-1] = cat_dir_int[i];
		}
	}

	// free(cat_dir_int);

	tree->primary->nextsplit = NULL;


/*

	if (rp.numcat[tree->primary->var_num] == 0)
	{
		tree->primary->spoint = splitpoint[0];
	}
*/

	if (0 < rp.maxsur)
		surrogate(tree, 0, n);


	// Return the body of the tree
	// For each component we first create a vector to hold the
	// result, then a ragged array index into the vector.
	// The rpmatrix routine then fills everything in.

	rpcountup(tree, nsplit);
	int splitcnt = nsplit[0] + nsplit[1];

	dsplit3 = PROTECT(allocMatrix(REALSXP, splitcnt, 2));
	dptr = REAL(dsplit3);

	for (i = 0; i < 2; i++) {
		ddsplit[i] = dptr;
		dptr += splitcnt;
		for (j = 0; j < splitcnt; j++) {
			ddsplit[i][j] = 0.0;
		}
	}


	isplit3 = PROTECT(allocMatrix(INTSXP, splitcnt, 2));
	iptr = INTEGER(isplit3);

	for (i = 0; i < 2; i++) {
		iisplit[i] = iptr;
		iptr += splitcnt;
	}

	rpmatrix(tree, rp.numcat, ddsplit, iisplit, nsplit);

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