// This routine is from rpart
//
// Calculate the surrogate splits for a node and its primary
//    (This routine is an awful lot like bsplit)
//
// Input :      node
//              start and stop indices for the arrays (which obs apply)
//
// Output:      Fills in the node's
//                      surrogate splits
//                      lastsurrogate value
//
// Uses:        The global vector tempvec (integer) as a temporary, assumed
//                to be of length n.

#include "get_surg.h"
#include "node.h"
#include "rpartproto.h"

void surrogate(pNode me, int n1, int n2) {
	int i, j, k;
	// the primary split variable
	int var;
	double split;
	double improve;
	// weight sent left and right by primary
	double lcount, rcount;
	int extra;
	pSplit ss;
	pSplit *ss_list;
	int *index;
	int *tempy;
	int **sorts;
	double **xdata;
	double adj_agree;
	int ncat;

	// DMA 
	tempy = (int *) calloc(rp.n, sizeof(int));
	sorts = rp.sorts;
	xdata = rp.xdata;

	// First construct, in tempy, the "y" variable for this calculation.
	// It will be LEFT:goes left, 0:missing, RIGHT:goes right.
	// Count up the number of obs the primary sends to the left, as my
	// last surrogate (or to the right, if larger).

	var = (me->primary)->var_num;

	// LCJ
	// continuous variable
	if (rp.numcat[var] == 0)
	{
		split = (me->primary)->spoint;
		extra = (me->primary)->csplit[0];
		for (i = n1; i < n2; i++)
		{	
			// lcj: j is the variable num out of the sorted matrix of variables "sorts", whereby "sorts" is sorted in ascending order
			j = sorts[var][i];
			// printf("j (von primvar == 8): %d\n", j);
			// -> vgl. j mit logfile branch_2_sorts.txt: variable numbers stimmen überein
			if (j < 0)
			{
				tempy[-(j + 1)] = 0;
			} else
			{
				tempy[j] = (xdata[var][j] < split) ? extra : -extra;
			}
		}
	} else
	{ /* categorial variable */
	// lcj: index is a list of directions for the cats, whereby the index of "index" represents the categorial and the corresponding value the direction
	// example: index = [1, -1, -1, 1]: cat 0 goes 1 (right), cat 1 goes -1 (left), cat 2 goes -1 (left), cat 3 goes 1 (right)
	index = (me->primary)->csplit;
	
	for (i = n1; i < n2; i++) {
		j = sorts[var][i];
	    if (j < 0) {
		tempy[-(j + 1)] = 0;
		}
	    else {
		tempy[j] = index[(int) xdata[var][j] - 1];
		}
	
	}
	
	}


	lcount = 0;
	rcount = 0;
	// lcj: var: PrimVar (see line 48)
	for (i = n1; i < n2; i++) {
		j = sorts[var][i];
		//if (j < 0)
		//	j = -(j + 1);
		j = j < 0 ? -(j + 1) : j;
		switch (tempy[j]) {
		case LEFT:
			lcount += rp.wt[j];
			break;
		case RIGHT:
			rcount += rp.wt[j];
			break;
		default:
			break;
		}
	}
	// end parallel section

	if (lcount < rcount)
		me->lastsurrogate = RIGHT;
	else {
		if (lcount > rcount)
			me->lastsurrogate = LEFT;
		else
			// no default
			me->lastsurrogate = 0;
	}

	// Now walk through the variables
	me->surrogate = (pSplit) NULL;
	//int splitLR = rp.csplit;
	int nthreads;

		for (i = 0; i < rp.nvar; i++) {
			if (var == i)
				continue;
			ncat = rp.numcat[i];

			

			choose_surg(n1, n2, tempy, xdata[i], sorts[i], ncat, &improve, &split, rp.csplit_for_cat, lcount, rcount, &adj_agree);

			// org comment: was 0
			if (adj_agree <= 1e-10)
				// no better than default
				continue;
			

			// sort it onto the list of surrogates

			ss = insert_split(&(me->surrogate), ncat, improve, rp.maxsur);
			if (ss) {
				ss->improve = improve;
				ss->var_num = i;
				// corrected by nodesplit()
				ss->count = 0;
				ss->adj = adj_agree;
				if ( rp.numcat[i] == 0) {
					ss->spoint = split;
					ss->csplit[0] = rp.csplit_for_cat[0];
				} else {
					for (k = 0; k < rp.numcat[i]; k++) {
						ss->csplit[k] = rp.csplit_for_cat[k];
					}
				}
				
			}

		}
}
