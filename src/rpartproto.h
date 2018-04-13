/*
 * This routine is from rpart
 *
 * prototypes for all of the functions
 *   This helps the ansi compiler do tight checking.
 *
 */
#include "node.h"

void choose_surg(int n1, int n2, int *y, double *x, int *order, int ncat, double *agreement, double *split, int *csplit, double ltot, double rtot, double *adj);

void free_tree(pNode node, int freenode);

pSplit insert_split(pSplit *listhead, int ncat, double improve, int max);

void sort_vec(int start, int stop, double *x, int *cvec);

void rpcountup(pNode me, int *nsplit);

void rpmatrix(pNode me, int *numcat, double **dsplit, int **isplit, int **csplit, int id);

void surrogate(pNode me, int n1, int n2);
