// This routine is originally from rpart.
//
// count up the number of nodes and splits in the final result
//

#include "get_surg.h"
#include "node.h"
#include "rpartproto.h"

void rpcountup(pNode me, int *nsplit) {
	int i, j;
	pSplit ss;

	i = 0;
	j = 0;
#ifdef OPENMP_ON
	#pragma omp sections private (ss)
// start parallel section
	{
		#pragma omp section
#endif
		for (ss = me->primary; ss; ss = ss->nextsplit) {
			i++;
		}
#ifdef OPENMP_ON
		#pragma omp section
#endif
		for (ss = me->surrogate; ss; ss = ss->nextsplit) {
			j++;
		}
#ifdef OPENMP_ON
	}
// end parallel section
#endif
	nsplit[0] = i;
	nsplit[1] = j;
}