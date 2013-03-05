
/* mexgspan.cpp - MEX interface to Taku Kudo's gspan implementation and
 *    additionally gspan for graph boosting.
 *
 * Author: Sebastian Nowozin <sebastian.nowozin@tuebingen.mpg.de>
 * Date: 31st August 2006
 *
 * Usage for gspan frequent subgraph mining:
 *
 * [graphs, count] = mexgspan (G, minsup, size_req, directed);
 *
 * or for graph boosting
 *
 * [graphs, ybase, GY] = mexgspan (G, minsup, size_req, directed,
 *    boostY, boostWeights, boostTau, boostN, boostMax, boostType);
 *
 * Find all subgraphs in the set of graphs G that appear at least minsup
 * times.  Return at most maxpat patterns.
 *
 ***
 * Input
 *
 * G: (1,n) cellarray of n graph structures with this layout
 *    g.nodelabels: (n,1) discrete integer labels [L_1 ; L_2 ; ... ; L_n];
 *    g.edges: (m,2) edges, [from to] at each line:
 *       [e_1_{from} e_1_{to} edgelabel_1 ; ... ; e_m_{from} e_m_{to} edgelabel_m]
 *       The node indices go from 1 to n.
 * minsup: The minimum number of times a frequent pattern has to appear in G.
 * size_req: (1,2) real vector containing lower/upper bounds on graph size.
 *    Use [0 0] for no limit, [n 0] for a lower bound, [0 n] for an upper bound.
 * directed: 1 in case the graphs are directed, 0 in case of undirected.
 *    Undirected graphs are working, directed graphs return the wrong number
 *    of substructes(?).
 *
 * (optional, only necessary for graph boosting)
 *
 * boostY: (n,1) array of labels in {-1, 1}.
 * boostWeights: (n,1) array of reals in [0 ; 1].
 * boostTau: tau parameter (normally set to 0).
 * boostN: integer > 0.  From the best patterns, only the best boostN patterns
 *    are returned.
 * boostMax: The maximum number of exploration steps to take or zero for
 *    infinity (1e6 is a good choice here).
 * boostType: 1 for 1-class and 1.5-class boosting, 2 for 2-class boosting.
 *
 ***
 * Output
 *
 * For normal frequent subgraph mining:
 *
 * graphs: (1,p) cellarray of p graph structures that appear frequently in G.
 * count: (1,p) uint32 array giving at element i the number of times graphs{i}
 *    appears in G (duplicates within one graph only counted once).
 *
 * For boosting:
 *
 * graphs: (1,p) cellarray of p graph structures that are discriminative in G.
 * ybase: (1,p) integers in {-1, 1} for 2-class boosting, {1} for 1-class.
 * GY: (n,p) response matrix
 */

#include <assert.h>
#include <limits.h>

#include <mex.h>

#include "gspan.h"

/* [graphs, count] = mexgspan (G, minsup, maxpat);
 */
void
mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	if (nlhs > 3) {
		mexPrintf ("Too many left hand side parameters.\n");
		return;
	}

	if (nrhs != 4 && nrhs != 10) {
		mexErrMsgTxt ("Wrong number of right hand side parameters (must be 4 or 10).\n");

		return;
	}

	if (mxIsCell (prhs[0]) == 0) {
		mexErrMsgTxt ("G must be a cellarray of graph structures.\n");

		return;
	}

	unsigned int minsup = (unsigned int) mxGetScalar (prhs[1]);
	if (mxGetM (prhs[2]) != 1 || mxGetN (prhs[2]) != 2) {
		mexErrMsgTxt ("size_req must be a (1,2) real vector.\n");

		return;
	}
	double* size_req = (double *) mxGetPr (prhs[2]);
	unsigned int maxpat_min = (unsigned int) size_req[0];
	unsigned int maxpat_max = (unsigned int) size_req[1];
	unsigned int directed = (unsigned int) mxGetScalar (prhs[3]);
	if (directed > 1) {
		mexErrMsgTxt ("Parameter 'directed' must be zero or 1.\n");
		return;
	}

	GSPAN::gSpan gspan;

	/* In case we should do graph boosting, we have to do some extra setup
	 * work, such as initializing the known labels of the graphs and their
	 * weights.
	 */
	if (nrhs == 10) {
		double boostTau = (double) mxGetScalar (prhs[6]);
		if (boostTau < 0.0) {
			mexErrMsgTxt ("Boost tau parameter must be >= 0.0.\n");
			return;
		}
		unsigned int boostN = (unsigned int) mxGetScalar (prhs[7]);

		double* boostYptr = (double *) mxGetPr (prhs[4]);
		double* boostWeightsPtr = (double *) mxGetPr (prhs[5]);
		unsigned int boostYlen = mxGetM (prhs[4]);
		std::vector<double> boostY;
		std::vector<double> boostWeights;
		for (unsigned int n = 0 ; n < boostYlen ; ++n) {
			boostY.push_back (boostYptr[n]);
			boostWeights.push_back (boostWeightsPtr[n]);
		}
		unsigned int boostMax = (unsigned int) mxGetScalar (prhs[8]);
		unsigned int boostType = (unsigned int) mxGetScalar (prhs[9]);

		gspan.boost_setup (boostN, boostTau, boostMax, boostY,
			boostWeights, boostType);
	}

	/* Start the actual gspan algorithm.  Before the graphs are converted from
	 * the Matlab structures to Taku's graph object format.
	 */
	gspan.run_graphs (prhs[0], nlhs, plhs,
		minsup, maxpat_min, maxpat_max,
		false, false, (directed == 1) ? true : false);
}


