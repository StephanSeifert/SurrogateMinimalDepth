// This routine is originally from rpart.
//
// For S's usage, convert the linked list data into matrix form

#include "get_surg.h"
#include "node.h"
#include "rpartproto.h"

void rpmatrix(pNode me, int *numcat, double **dsplit, int **isplit,  int *nsplit) {

	// dsplit  0: improvement
	//         1: split point for continuous
	//         2: surrogate: adjusted agreement,  primary: nothing
	// isplit  0: variable #
	//         1: count
	//         2: continuous: direction -1=left, 1=right


	int scnt;
	int j;
	pSplit spl;

	scnt = 0;


			spl = me->primary;
			for (scnt = 0; scnt < nsplit[0]; scnt++) {
				j = spl->var_num;
				
				if (numcat[j] == 0)
				{
					dsplit[1][scnt] = spl->spoint;
					isplit[1][scnt] = -1;
				} else
				{
					isplit[1][scnt] = numcat[j];
				}

				isplit[0][scnt] = spl->var_num + 1;
				

				spl = spl->nextsplit;
			}

			spl = me->surrogate;
			int splitcnt = nsplit[0] + nsplit[1];
			for (scnt = nsplit[0]; scnt < splitcnt; scnt++) {
				j = spl->var_num;
				
				dsplit[0][scnt] = spl->adj;
				if (numcat[j] == 0) {
					dsplit[1][scnt] = spl->spoint;
					isplit[1][scnt] = -1;
				} else {
					isplit[1][scnt] = numcat[j];
				}
				isplit[0][scnt] = j + 1;
				


				spl = spl->nextsplit;
			}

}