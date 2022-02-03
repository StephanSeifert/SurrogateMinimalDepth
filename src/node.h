#ifndef RPART_NODE_H
#define RPART_NODE_H

// This routine is originally from rpart.
//
// definition of a node in the tree
//
// The actual size of these structures when allocated in insert_split.c
// depends on the split.
// csplit[0] gets used even for continuous splits.

typedef struct split {
	double improve;
	// for surrogates only, adjusted agreement
	double adj;
	// only used if it is continuous
	double spoint;
	struct split *nextsplit;
	int var_num;
	int count;
	// the actual length depends on splitting rule
	int csplit[20];
} Split, *pSplit;

typedef struct node {
	// risk for the node
	double risk;
	pSplit primary, surrogate;
	int lastsurrogate;
} Node, *pNode;

//
//  Split:
//      variable number of the split; 0 = no more surrogates (or primaries)
//
//      split point: the actual split point for a continuous
//
//      improve:  For primary splits, the improvement index returned by the
//                bsplit routine. This is the measure that determines the
//                winning split.
//                For surrogate splits, this holds the error rate, i.e., the
//                % incorrect guesses of the primary by using this surrogate.
//
//      count: The number of observations split using this variable. For the
//             first primary, this will = the number of non-missing values.
//             For surrogates, it will be the number missing in the primary
//             and all earlier surrogates but not missing on this one. (For
//             all primaries but the first, the number is theoretical).
//
//   adj:  Let "maj" be the %agreement for going with the majority,
//         and "agree" the %agreement for this surrogate. The
//         adjusted value is (agree - maj)/(1-maj); the amount of
//         the potential improvement actually realized. The denominator
//         for both percents depends on the sur_agree option.
//
//      csplit[0]:   For a continuous variable, we also need to know the
//                   direction of the split. We use this "extra" variable
//                   as 1: <x to the left, -1: <x to the right.
//
//      csplit[]:    For a categorical, the labels are LEFT, RIGHT, and
//                   0=missing. (Even if a particular category is not empty,
//                   there may be no subjects from that category present
//                   at a particular split further down the tree).
//
//
//  Node:
//      num_obs: Number of observations in the node.
//
//      risk: From the eval routine. Estimate of risk, if this node were
//            terminal.
//
//      lastsurrogate: Which direction to send obs for which the primary and
//              all the surrogates are missing. (The child with the greatest
//             sum of weights).

#endif
