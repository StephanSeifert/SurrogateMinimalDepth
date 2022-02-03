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

		for (ss = me->primary; ss; ss = ss->nextsplit) {
			i++;
		}

		for (ss = me->surrogate; ss; ss = ss->nextsplit) {
			j++;
		}

	nsplit[0] = i;
	nsplit[1] = j;
}