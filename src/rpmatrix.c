// This routine is originally from rpart.
//
// For S's usage, convert the linked list data into matrix form

#include "get_surg.h"
#include "node.h"
#include "rpartproto.h"

void rpmatrix(pNode me, double **dsplit, int **isplit,  int *nsplit) {

	// dsplit  0: improvement
	//         1: split point for continuous
	//         2: surrogate: adjusted agreement,  primary: nothing
	// isplit  0: variable #
	//         1: count
	//         2: continuous: direction -1=left, 1=right


	int scnt;
	pSplit spl;

	scnt = 0;

#ifdef OPENMP_ON
	#pragma omp sections private(scnt,spl)
// start parallel section
	{
		#pragma omp section
		{
#endif
			spl = me->primary;
			for (scnt = 0; scnt < nsplit[0]; scnt++) {
				dsplit[0][scnt] = spl->improve;
				dsplit[1][scnt] = spl->spoint;

				isplit[0][scnt] = spl->var_num + 1;
				isplit[1][scnt] = spl->count;
				isplit[2][scnt] = spl->csplit;

				spl = spl->nextsplit;
			}
#ifdef OPENMP_ON
		}
		#pragma omp section
		{
#endif
			spl = me->surrogate;
			int splitcnt = nsplit[0] + nsplit[1];
			for (scnt = nsplit[0]; scnt < splitcnt; scnt++) {
				dsplit[0][scnt] = spl->improve;
				dsplit[1][scnt] = spl->spoint;
				dsplit[2][scnt] = spl->adj;


				isplit[0][scnt] = spl->var_num + 1;
				isplit[1][scnt] = spl->count;
				isplit[2][scnt] = spl->csplit;

				spl = spl->nextsplit;
			}
#ifdef OPENMP_ON
		}
	}
// end parallel section
#endif
}