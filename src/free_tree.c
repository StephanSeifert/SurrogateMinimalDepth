// This routine is from rpart
//
// free up all of the memory associated with a tree

#include "get_surg.h"
#include "node.h"
#include "rpartproto.h"

void free_split(pSplit spl) {
	if (spl) {
		free_split(spl->nextsplit);
		free(spl);
	}
}

// use freenode if the tree was CALLOC-ed, from xval.c
void free_tree(pNode node, int freenode) {
	free_split(node->surrogate);
	free_split(node->primary);

	if (freenode == 1)
		if (node)
			free(node);
		else {
			// don't point to things I just freed
			node->primary = (pSplit) NULL;
			node->surrogate = (pSplit) NULL;
		}
}