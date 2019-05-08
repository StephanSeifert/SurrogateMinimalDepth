/*
 * This routine is from rpart
 *
 * A particular split routine, optimized for the surrogate variable
 *  search.  The "goodness" of a split is the total weights of concordant
 *  observations between the surrogate and the primary split.
 *  Note that the CART folks use the %concordance, which factors missing
 *  values into the equations somewhat differently.
 *
 *  y is coded as  +1=left, -1=right, 0=missing
 *
 */
#include "get_surg.h"
#include "rpartproto.h"

void choose_surg(int n1, int n2, int *y, double *x, int *order, double *agreement, double *split, int *csplit, double tleft, double tright, double *adj) {
	/* sum of weights for each */
	double llwt, lrwt, rrwt, rlwt;
	double agree, majority, total_wt;
	/* set to 1 when something worthwhile is found */
	int success = 0;

	/*
	 * I enforce that at least 2 obs must go each way, to avoid having an
	 *  uncorrelated surrogate beat the "null" surrogate too easily
	 * Observations with 0 weight don't count in this total
	 */
	/* continuous case */
	/*
	 * ll = y's that go left that are also sent left by my split
	 * lr = y's that go left that I send right
	 * rl = y's that go right that I send to the left
	 * rr = y's that go right that I send to the right
	 *
	 * The agreement is max(ll+rr, lr+rl), if weights were = 1;
	 *   actually max(llwt + rrwt, lrwt + rlwt) / denominator
	 */
	double lastx = 0.0;
	int ll, lr, rr, rl;
	ll = rl = 0;
	llwt = 0;
	rlwt = 0;

	for (int i = n2 - 1; i >= n1; i--) {
		/* start with me sending all to the left */
		int j = order[i];
		if (j >= 0) {
			/* this is why I run the loop backwards */
			lastx = x[j];
			switch (y[j]) {
			case LEFT:
				if (rp.wt[j] > 0)
					ll++;
				llwt += rp.wt[j];
				break;
			case RIGHT:
				if (rp.wt[j] > 0)
					rl++;
				rlwt += rp.wt[j];
				break;
			default:
				;
			}
		}
	}

	if (llwt > rlwt)
		agree = llwt;
	else
		agree = rlwt;

	/*
	 *  Now we have the total agreement.  Calculate the %agreement and
	 *    the adjusted agreement
	 *  For both, do I use the total y vector as my denominator (my
	 *    preference), or only the y's for non-missing x (CART book)?
	 *    If the former, need to reset some totals.
	 */
	if (rp.sur_agree == 0) {
		/* use total table */
		total_wt = tleft + tright;
		if (tleft > tright)
			majority = tleft;
		else
			majority = tright;
	} else {
		/* worst possible agreement */
		majority = agree;
		total_wt = llwt + rlwt;
	}

	lr = rr = 0;
	lrwt = 0;
	rrwt = 0;
	/*
	 *  March across, moving things from the right to the left
	 *    the "lastx" code is caring for ties in the x var
	 *    (The loop above sets it to the first unique x value).
	 */
	/* NB: might never set csplit or split */
	*csplit = LEFT;
	/* a valid splitting value */
	*split = lastx;
	for (int i = n1; (ll + rl) >= 2; i++) {
		int j = order[i];
		if (j >= 0) {
			/* not a missing value */
			if ((lr + rr) >= 2 && x[j] != lastx) {
				/* new x found, evaluate the split */
				if ((llwt + rrwt) > agree) {
					success = 1;
					agree = llwt + rrwt;
					/* < goes to the right */
					*csplit = RIGHT;
					*split = (x[j] + lastx) / 2;
				} else if ((lrwt + rlwt) > agree) {
					success = 1;
					agree = lrwt + rlwt;
					*csplit = LEFT;
					*split = (x[j] + lastx) / 2;
				}
			}

			switch (y[j]) {
			/* update numbers */
			case LEFT:
				if (rp.wt[j] > 0) {
					ll--;
					lr++;
				}
				llwt -= rp.wt[j];
				lrwt += rp.wt[j];
				break;
			case RIGHT:
				if (rp.wt[j] > 0) {
					rl--;
					rr++;
				}
				rlwt -= rp.wt[j];
				rrwt += rp.wt[j];
				break;
			default:
				/* ignore missing y's */
				;
			}
			lastx = x[j];
		}
	}

	/*
	 * success = 0 means no split was found that had at least 2 sent each
	 *   way (not counting weights of zero), and for continuous splits also had
	 *   an improvement in agreement.
	 *  Due to round-off error such a split could still appear to have adj>0
	 *  in the calculation further below; avoid this.
	 */

	if (!success) {
		*agreement = 0.0;
		*adj = 0.0;
		return;
	}

	*agreement = agree / total_wt;
	majority /= total_wt;
	*adj = (*agreement - majority) / (1. - majority);
}
