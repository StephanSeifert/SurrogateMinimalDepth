// This routine is originally from rpart
//
// prototypes for all of the functions
//   This helps the ansi compiler do tight checking.
//

#include "node.h"

void choose_surg(int n1, int n2, int *y, double *x, int *order, int ncat, double *agreement, double *split, int *csplit_for_cat, double ltot, double rtot, double *adj);

void free_split(pSplit spl);

void free_tree(pNode node, int freenode);

pSplit insert_split(pSplit *listhead, int ncat, double improve, int max);

void sort_vec(int start, int stop, double *x, int *cvec, int var);

void rpcountup(pNode me, int *nsplit);

void rpmatrix(pNode me, int *numcat, double **dsplit, int **isplit,  int *splits);

void surrogate(pNode me, int n1, int n2);