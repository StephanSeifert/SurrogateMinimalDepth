/*
 * This routine is orignaly from rpart.
 *
 *  For S's usage, convert the linked list data into matrix form
 */
#include "get_surg.h"
#include "node.h"
#include "rpartproto.h"

/* These four preserve from call to call */
static int ncnt, scnt, ccnt;
static double cp_scale;

void rpmatrix(pNode me, int *numcat, double **dsplit, int **isplit, int **csplit, int id) {
	/*
	 * dsplit  0: improvement
	 *         1: split point if continuous; index into csplit if not
	 *         2: surrogate: adjusted agreement,  primary: nothing
	 * isplit  0: variable #
	 *         1: count
	 *         2: if continuous: direction -1=left, 1=right
	 *            if categorical: # of categories
	 * csplit[i]: -1=left, 0=missing category, 1=right
	 */

	int i, j, k;
	pSplit spl;

	if (id == 1) {
		/* this is the top node */
		cp_scale = 1 / me->risk;
		scnt = 0;
		ncnt = 0;
		ccnt = 0;
	}

	i = 0;
	for (spl = me->primary; spl; spl = spl->nextsplit) {
		i++;
		j = spl->var_num;
		dsplit[0][scnt] = spl->improve;
		if (numcat[j] == 0) {
			dsplit[1][scnt] = spl->spoint;
			isplit[2][scnt] = spl->csplit[0];
		} else {
			/* categorical */
			dsplit[1][scnt] = ccnt + 1;
			isplit[2][scnt] = numcat[j];
			if (csplit != NULL)
				for (k = 0; k < numcat[j]; k++)
					csplit[k][ccnt] = spl->csplit[k];
			ccnt++;
		}
		/* use "1" based subscripts */
		isplit[0][scnt] = j + 1;
		isplit[1][scnt] = spl->count;
		scnt++;
	}

	i = 0;
	for (spl = me->surrogate; spl; spl = spl->nextsplit) {
		i++;
		j = spl->var_num;
		dsplit[0][scnt] = spl->improve;
		dsplit[2][scnt] = spl->adj;
		if (numcat[j] == 0) {
			dsplit[1][scnt] = spl->spoint;
			isplit[2][scnt] = spl->csplit[0];
		} else {
			dsplit[1][scnt] = ccnt + 1;
			isplit[2][scnt] = numcat[j];
			if (csplit != NULL)
				for (k = 0; k < numcat[j]; k++)
					csplit[k][ccnt] = spl->csplit[k];
			ccnt++;
		}
		isplit[0][scnt] = j + 1;
		isplit[1][scnt] = spl->count;
		scnt++;
	}

	ncnt++;

}
