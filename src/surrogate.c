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
	double rnd;
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

	// DMA 
	tempy = (int *) calloc(rp.n, sizeof(int));
	sorts = rp.sorts;
	xdata = rp.xdata;

	// First construct, in tempy, the "y" variable for this calculation.
	// It will be LEFT:goes left, 0:missing, RIGHT:goes right.
	// Count up the number of obs the primary sends to the left, as my
	// last surrogate (or to the right, if larger).

	var = (me->primary)->var_num;
	// continuous variable
	split = (me->primary)->spoint;
	extra = (me->primary)->csplit;
	for (i = n1; i < n2; i++) {
		j = sorts[var][i];
		if (j < 0)
			tempy[-(j + 1)] = 0;
		else
			tempy[j] = (xdata[var][j] < split) ? extra : -extra;
	}

	// count the total number sent left and right
	lcount = 0;
	rcount = 0;
#ifdef OPENMP_ON
	#pragma omp parallel for private(i,j) schedule(dynamic) reduction(+:lcount,rcount) //TODO num_threads(rp.nthreads)
	// start parallel section
#endif
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
	int splitLR = rp.csplit;
	int nthreads;
#ifdef OPENMP_ON
	#pragma omp parallel private(i,adj_agree,rnd,improve,split,splitLR,ss) num_threads(rp.nthreads)
	{
		// start parallel section
		#pragma omp single
		{
			nthreads = omp_get_num_threads();
			// DMA 
			ss_list = (pSplit*) calloc(nthreads, sizeof(pSplit));

		}
		#pragma omp for schedule(dynamic)
#endif
		for (i = 0; i < rp.nvar; i++) {
			if (var == i)
				continue;

			choose_surg(n1, n2, tempy, xdata[i], sorts[i], &improve, &split, &splitLR, lcount, rcount, &adj_agree);

			
			// org comment: was 0
			if (adj_agree <= 1e-10){
				// no better than default
				continue;
			}
			
				// sort it onto the list of surrogates
#ifdef  OPENMP_ON
				ss = insert_split(&(ss_list[omp_get_thread_num()]), &rnd, improve, rp.maxsur);
#else
				ss = insert_split(&(me->surrogate), &rnd, improve, rp.maxsur);
#endif
				if (ss) {
					ss->improve = improve;
					ss->rnd = rnd;
					ss->var_num = i;
					// corrected by nodesplit()
					ss->count = 0;
					ss->adj = adj_agree;
					ss->spoint = split;
					ss->csplit = splitLR;
				}


		}
#ifdef OPENMP_ON
	}
// end parallel section
	for (i = 0; i < nthreads; ++i) {
		if (ss_list[i] != NULL){
			for (pSplit lh = ss_list[i]; lh;lh = lh->nextsplit) {
					ss = insert_split(&(me->surrogate), &(lh->rnd), lh->improve, rp.maxsur);
					if (ss) {
						ss->improve = lh->improve;
						ss->var_num = lh->var_num;
						// corrected by nodesplit()
						ss->count = lh->count;
						ss->adj = lh->adj;
						ss->spoint = lh->spoint;
						ss->csplit = lh->csplit;
					}
			}
		}
	}
	// clear DMA
	for (i = 0; i < nthreads; ++i)
		free_split(ss_list[i]);
	free(ss_list);
	free(tempy);
#endif
}