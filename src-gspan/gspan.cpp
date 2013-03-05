/*
    $Id: gspan.cpp,v 1.8 2004/05/21 09:27:17 taku-ku Exp $;
 
   Copyright (C) 2004 Taku Kudo, All rights reserved.
     This is free software with ABSOLUTELY NO WARRANTY.
  
   This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.
    
   This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
    
   You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
     02111-1307, USA
*/
#include "gspan.h"
#include <iterator>

#include <stdlib.h>
#include <unistd.h>

#include <mex.h>

namespace GSPAN {

gSpan::gSpan (void)
{
	boost = false;
}

void gSpan::mexAppendGraph (Graph* g, unsigned int ID, double par,
	std::map<unsigned int, unsigned int>& counts)
{
	mexOutputG* ms = new mexOutputG ();
	ms->next = NULL;
	ms->graph = g->writemex ();
	ms->par = par;
	ms->counts = counts;	// copy

	// Chainlink at the end.
	if (mexOutputLast == NULL) {
		mexOutputRoot = mexOutputLast = ms;
	} else {
		mexOutputLast->next = ms;
		mexOutputLast = ms;
	}
}

void gSpan::read (const mxArray* graphs)
{
	unsigned int count = mxGetN (graphs);
	const mxArray** garr = (const mxArray**) mxGetPr (graphs);

	for (unsigned int n = 0 ; n < count ; ++n) {
		const mxArray* gcur = garr[n];

		Graph g(directed);
		g.read (gcur);
#ifdef	DEBUG
		g.check ();
#endif
		TRANS.push_back (g);
	}
}

std::istream &gSpan::read (std::istream &is)
{
	Graph g(directed);
	while (true) {
		g.read (is);
		if (g.empty()) break;
		TRANS.push_back (g);
	}
	return is;
}


/* 2-class LPBoosting gain function
 */
double
gSpan::gain2 (Projected &projected, double y)
{
	/* Compute actual gain (eq. (10) in [Dimiriz2002]) and Problem 1 in graph
	 * boost paper.
	 *
	 * For L samples x_i, calculate
	 *   gain = \sum_{i=1}^{L} d_i y_i h(x_i)
	 * where h is implicitly defined by this current subgraph pattern.
	 */
	unsigned int oid = 0xffffffff;
	double gainm = 0.0;

	Projected::iterator cur = projected.begin();
	for ( ; cur != projected.end() ; ++cur) {
		if (oid != cur->id) {
			/* A positive pattern: it appears in the graph.
			 *
			 * h_{<t,y>} = y, if t \subseteq x,
			 *           = -y, otherwise.
			 *
			 * t is implicitly defined by projected.
			 * First, process all skipped graphs.
			 */
			for (unsigned int skipn = oid+1 ; skipn < cur->id ; ++skipn)
				gainm += boostY[skipn]*boostWeights[skipn]*(-y);

			/* Now process positive pattern.
			 * This is also the gain function for 1.5 LP boosting, as we have
			 * the constraint
			 *   \sum_{n=1}^N \lambda_n h_j(x_{1,n})
			 *   - \sum_{m=1}^M \mu_m h_j(x_{2,m}) <= \gamma.
			 * In order to maximize the violation of the current ensemble, we
			 * search for the classifier that maximizes the left side.
			 */
			gainm += boostY[cur->id]*boostWeights[cur->id]*y;
		}
		oid = cur->id;
	}

	/* All skipped graphs to the end.
	 */
	for (unsigned int skipn = oid+1 ; skipn < boostY.size() ; ++skipn)
		gainm += boostY[skipn]*boostWeights[skipn]*(-y);

	return (gainm);
}


/* 1.5-class LPBoosting gain function
 */
double
gSpan::gain1d5 (Projected &projected, double y)
{
	/* Compute actual gain (eq. (10) in [Dimiriz2002]) and Problem 1 in graph
	 * boost paper.
	 *
	 * For L samples x_i, calculate
	 *   gain = \sum_{i=1}^{L} d_i y_i h(x_i)
	 * where h is implicitly defined by this current subgraph pattern.  In the
	 * 1.5-class case, h(x_i) is zero if the pattern does not appear, hence we
	 * can skip those instances silently.
	 */
	unsigned int oid = 0xffffffff;
	double gain = 0.0;

	Projected::iterator cur = projected.begin();
	for ( ; cur != projected.end() ; ++cur) {
		if (oid != cur->id)
			gain += boostY[cur->id]*boostWeights[cur->id]*y;

		oid = cur->id;
	}

	return (gain);
}

/* Gain bound for 1.5-class formulation
 */
double
gSpan::gainbound1d5 (Projected &projected)
{
	double gain = 0.0;
	unsigned int oid = 0xffffffff;

	for (Projected::iterator cur = projected.begin() ;
		cur != projected.end() ; ++cur)
	{
		if (oid != cur->id) {
			oid = cur->id;
			if (boostY[cur->id] <= 0.0)
				continue;

			gain += boostWeights[cur->id];
		}
	}

	return (gain);
}

/* Gain bound for 2-class formulation
 */
double
gSpan::gainbound2 (Projected &projected)
{
	/* Lemma 1, section 3.2 in the graph boosting paper.  We automatically
	 * only enumerate those patterns which contain the subgraph 'projected'.
	 *
	 * <t^,y^> = argmax_{t \in F,y \in {+/- 1}} d_i y_i h_{<t,y>}, where
	 * F = \unionset_{i=1}^{L} { t | t \subseteq x_i }.
	 */
	unsigned int oid = 0xffffffff;
	unsigned int size = 0;
	double gain_pos = 0.0;
	double gain_neg = 0.0;

	for (Projected::iterator cur = projected.begin() ;
		cur != projected.end() ; ++cur)
	{
		if (oid != cur->id) {
			if (boostY[cur->id] <= 0.0)
				gain_neg += boostWeights[cur->id];
			else
				gain_pos += boostWeights[cur->id];

			++size;
		}
		oid = cur->id;
	}

	/* The boostWeightSum has been precomputed as it is always the same for
	 * one run.
	 */
	gain_neg = 2.0*gain_neg + boostWeightSum;
	gain_pos = 2.0*gain_pos - boostWeightSum;

	if (gain_neg >= gain_pos)
		return (gain_neg);

	return (gain_pos);
}


std::map<unsigned int, unsigned int>
gSpan::support_counts (Projected &projected)
{
	std::map<unsigned int, unsigned int> counts;

	for (Projected::iterator cur = projected.begin() ;
		cur != projected.end() ; ++cur)
	{
		counts[cur->id] += 1;
	}

	return (counts);
}


unsigned int
gSpan::support (Projected &projected)
{
	unsigned int oid = 0xffffffff;
	unsigned int size = 0;

	for (Projected::iterator cur = projected.begin(); cur != projected.end(); ++cur) {
		if (oid != cur->id) {
			++size;
		}
		oid = cur->id;
	}

	return size;
}

void gSpan::report_boosting (Projected &projected, unsigned int sup,
	double gain, double yval)
{
	if (maxpat_max > maxpat_min && DFS_CODE.nodeCount () > maxpat_max)
		return;
	if (maxpat_min > 0 && DFS_CODE.nodeCount () < maxpat_min)
		return;

	Graph g(directed);
	DFS_CODE.toGraph (g);

	// insert individual counts
	std::map<unsigned int, unsigned int> GYcounts;
	unsigned int oid = 0xffffffff;
	for (Projected::iterator cur = projected.begin(); cur != projected.end(); ++cur) {
		if (oid != cur->id)
			GYcounts[cur->id] = 0;

		GYcounts[cur->id] += 1;
		oid = cur->id;
	}

	report_boosting_inter (g, sup, gain, yval, GYcounts);
}

void gSpan::report_boosting_inter (Graph &g, unsigned int sup,
	double gain, double yval, std::map<unsigned int, unsigned int>& GYcounts)
{
	g.mexprint ();

	/* Insertion sort based on the gain values.
	 */
	unsigned int insertidx;
	for (insertidx = 0 ; insertidx < bestGraphsGain.size () ; ++insertidx) {
		if (gain >= bestGraphsGain[insertidx])
			break;
	}

	bestGraphs.insert (bestGraphs.begin() + insertidx, g);
	bestGraphsY.insert (bestGraphsY.begin() + insertidx, yval);
	bestGraphsGain.insert (bestGraphsGain.begin() + insertidx, gain);
	bestGraphsCounts.insert (bestGraphsCounts.begin() + insertidx, GYcounts);

#ifdef DEBUG
	mexPrintf ("   inserted good output graph, gain %lf, sup %d, Y %lf, new size %d\n",
		gain, sup, yval, g.size());
#endif

	/* If the maximum allowed size is exceeded, prune it down.
	 */
	if (bestGraphsGain.size () > boostN) {
#ifdef DEBUG
		mexPrintf ("   one element removed at end with gain %lf\n",
			bestGraphsGain[bestGraphsGain.size()-1]);
#endif

		bestGraphs.pop_back ();
		bestGraphsY.pop_back ();
		bestGraphsGain.pop_back ();
		bestGraphsCounts.pop_back ();
	}

	++ID;
}


/* Special report function for single node graphs.
 */
void gSpan::report_single (Graph &g, std::map<unsigned int, unsigned int>& ncount)
{
	unsigned int sup = 0;
	for (std::map<unsigned int, unsigned int>::iterator it = ncount.begin () ;
		it != ncount.end () ; ++it)
	{
		sup += (*it).second;
	}

	if (maxpat_max > maxpat_min && g.size () > maxpat_max)
		return;
	if (maxpat_min > 0 && g.size () < maxpat_min)
		return;

	if (mex) {
		mexAppendGraph (&g, ID, sup, ncount);
		++ID;
	} else {
		if (enc == false) {
			if (where == false)
				*os << "t # " << ID << " * " << sup;
			*os << '\n';

			g.write (*os);
			*os << '\n';
		} else {
			std::cerr << "report_single not implemented for non-Matlab calls" << std::endl;
		}
	}
}


void gSpan::report (Projected &projected, unsigned int sup)
{
	/* Filter to small/too large graphs.
	 */
	if (maxpat_max > maxpat_min && DFS_CODE.nodeCount () > maxpat_max)
		return;
	if (maxpat_min > 0 && DFS_CODE.nodeCount () < maxpat_min)
		return;

	/* Output to matlab.
	 */
	if (mex) {
		Graph g(directed);
		DFS_CODE.toGraph (g);

		std::map<unsigned int, unsigned int> counts;
		unsigned int oid = 0xffffffff;
		for (Projected::iterator cur = projected.begin(); cur != projected.end(); ++cur) {
			if (oid != cur->id)
				counts[cur->id] = 0;

			counts[cur->id] += 1;
			//mexPrintf("counts[%d] = %d\n", cur->id, counts[cur->id]);
			oid = cur->id;
		}
		mexAppendGraph (&g, ID, sup, counts);
		++ID;
		return;
	}

	if (where) {
		*os << "<pattern>\n";
		*os << "<id>" << ID << "</id>\n";
		*os << "<support>" << sup << "</support>\n";
		*os << "<what>";
	}

	if (! enc) {
		Graph g(directed);
		DFS_CODE.toGraph (g);

		if (! where)
			*os << "t # " << ID << " * " << sup;

		*os << '\n';
		g.write (*os);
	} else {
		if (! where)
			*os << '<' << ID << ">    " << sup << " [";

		DFS_CODE.write (*os);
		if (! where) *os << ']';
	}

	if (where) {
		*os << "</what>\n<where>";
		unsigned int oid = 0xffffffff;
		for (Projected::iterator cur = projected.begin(); cur != projected.end(); ++cur) {
			if (oid != cur->id) {
				if (cur != projected.begin()) *os << ' ';
				*os << cur->id;
			}
			oid = cur->id;
		}
		*os << "</where>\n</pattern>";
	}

	*os << '\n';
	++ID;
}

/* Recursive subgraph mining function (similar to subprocedure 1
 * Subgraph_Mining in [Yan2002]).
 */
void gSpan::project (Projected &projected)
{
	if (boost && (boostmax > 0 && boostseen > boostmax))
		return;

	/* Check if the pattern is frequent enough.
	 */
	unsigned int sup = support (projected);
	if (sup < minsup)
		return;

	/* The minimal DFS code check is more expensive than the support check,
	 * hence it is done now, after checking the support.
	 */
	if (is_min () == false) {
		//      *os  << "NOT MIN [";  DFS_CODE.write (*os);  *os << "]" << std::endl;
		return;
	}

	if (boost) {
		/* Compute the gains for all classifiers we can build out of the current
		 * subgraph pattern (which is just two).
		 */
		double gain_ypos = 0.0;	// Gain if y=1 in the weak learner
		double gain_yneg = 0.0;	// Gain if y=-1 in the weak learner
		double gainmax = 0.0;	// Maximum of the negative/positive gain
		double yval = 1.0;	// y that achieves the maximum gain

		/* Compute gains.  This differs a little for the 1.5-class and the
		 * 2-class case.  (The gain function is replaced, as well).
		 */
		gain_ypos = (this->*gain) (projected, 1.0);
		if (boostType == 1) {
			gainmax = gain_ypos;
		} else if (boostType == 2) {
			gain_yneg = (this->*gain) (projected, -1.0);

			if (gain_ypos >= gain_yneg) {
				gainmax = gain_ypos;
			} else {
				gainmax = gain_yneg;
				yval = -1.0;
			}
		}

		/* Check if we have to replace the least best of our top performers.
		 */
		if (gainmax > boostTau) {
			/* Update best gain so far.  The best gain is also the tau bound.
			 */
#ifdef DEBUG
			mexPrintf ("boost gain observed with %lf, support %d\n", boostTau, sup);
#endif
			report_boosting (projected, sup, gainmax, yval);

			/* Update lower boosting gain bound in case we already exceeded
			 * the number of allowed patterns.
			 */
			if (bestGraphsGain.size () >= boostN) {
				boostTau = bestGraphsGain[bestGraphsGain.size() - 1];
				mexPrintf ("   top gains: %lf ... %lf\n",
					bestGraphsGain[0], boostTau);
			}
		}
		if (boostseen % 1000 == 0) {
			Graph g(directed);
			DFS_CODE.toGraph (g);

			mexPrintf ("%7u: gain %lf  this.bound %lf  "
				"global.bound %lf  DFS.size %u (gsize %d)\n", boostseen,
				gainmax, (this->*gainbound) (projected), boostTau,
				DFS_CODE.size(), g.size());
		}

		boostseen += 1;
	} else {
		// Output the frequent substructure
		report (projected, sup);
	}

	/* In case we have a valid upper bound and our graph already exceeds it,
	 * return.  Note: we do not check for equality as the DFS exploration may
	 * still add edges within an existing subgraph, without increasing the
	 * number of nodes.
	 */
	if (maxpat_max > maxpat_min && DFS_CODE.nodeCount () > maxpat_max)
		return;

	/* Compute gain bound for this pattern (projected) and only explore in
	 * case the bound allows this subgraph to be better.  That is, if the
	 * bound lets no supergraph be better than all reported ones so far, we
	 * just return.
	 */
	if (boost) {
		double gainb = (this->*gainbound) (projected);
#ifdef DEBUG
		mexPrintf ("   gain bound: %lf,  gain min required: %lf\n", gainb, boostTau);
		mexPrintf ("        actual gain for this subgraph: %lf\n", 
			(this->*gain) (projected, 1.0));
#endif

		if (gainb <= boostTau + 1e-8) {
#ifdef DEBUG
			mexPrintf ("      ==> return\n");
#endif
			return;
		}
	}

	/* We just outputted a frequent subgraph.  As it is frequent enough, so
	 * might be its (n+1)-extension-graphs, hence we enumerate them all.
	 */
	const RMPath &rmpath = DFS_CODE.buildRMPath ();
	int minlabel = DFS_CODE[0].fromlabel;
	int maxtoc = DFS_CODE[rmpath[0]].to;

	Projected_map3 new_fwd_root;
	Projected_map2 new_bck_root;
	EdgeList edges;

	/* Enumerate all possible one edge extensions of the current substructure.
	 */
	for (unsigned int n = 0; n < projected.size(); ++n) {

		unsigned int id = projected[n].id;
		PDFS *cur = &projected[n];
		History history (TRANS[id], cur);

		// XXX: do we have to change something here for directed edges?

		// backward
		for (int i = (int)rmpath.size()-1; i >= 1; --i) {
			Edge *e = get_backward (TRANS[id], history[rmpath[i]], history[rmpath[0]], history);
			if (e)
				new_bck_root[DFS_CODE[rmpath[i]].from][e->elabel].push (id, e, cur);
		}

		// pure forward
		// FIXME: here we pass a too large e->to (== history[rmpath[0]]->to
		// into get_forward_pure, such that the assertion fails.
		//
		// The problem is:
		// history[rmpath[0]]->to > TRANS[id].size()
		if (get_forward_pure (TRANS[id], history[rmpath[0]], minlabel, history, edges))
			for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
				new_fwd_root[maxtoc][(*it)->elabel][TRANS[id][(*it)->to].label].push (id, *it, cur);

		// backtracked forward
		for (int i = 0; i < (int)rmpath.size(); ++i)
			if (get_forward_rmpath (TRANS[id], history[rmpath[i]], minlabel, history, edges))
				for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it)
					new_fwd_root[DFS_CODE[rmpath[i]].from][(*it)->elabel][TRANS[id][(*it)->to].label].push (id, *it, cur);
	}

	/* Test all extended substructures.
	 */
	// backward
	for (Projected_iterator2 to = new_bck_root.begin(); to != new_bck_root.end(); ++to) {
		for (Projected_iterator1 elabel = to->second.begin(); elabel != to->second.end(); ++elabel) {
			DFS_CODE.push (maxtoc, to->first, -1, elabel->first, -1);
			project (elabel->second);
			DFS_CODE.pop();
		}
	}

	// forward
	for (Projected_riterator3 from = new_fwd_root.rbegin() ;
		from != new_fwd_root.rend() ; ++from)
	{
		for (Projected_iterator2 elabel = from->second.begin() ;
			elabel != from->second.end() ; ++elabel)
		{
			for (Projected_iterator1 tolabel = elabel->second.begin();
					tolabel != elabel->second.end(); ++tolabel)
			{
				DFS_CODE.push (from->first, maxtoc+1, -1, elabel->first, tolabel->first);
				project (tolabel->second);
				DFS_CODE.pop ();
			}
		}
	}

	return;
}

void gSpan::boost_setup (unsigned int _boostN, double _boostTau,
	unsigned int _boostmax,
	std::vector<double>& _boostY,
	std::vector<double>& _boostWeights,
	int _boostType)
{
	boost = true;
	boostseen = 0;
	boostmax = _boostmax;
	boostN = _boostN;
	boostTau = _boostTau;
	boostY = _boostY;
	boostWeights = _boostWeights;
	boostType = _boostType;

	/* Setup the corresponding function pointers for the gain function adn the
	 * gainbound function.
	 */
	if (boostType == 1) {
		gain = &gSpan::gain1d5;
		gainbound = &gSpan::gainbound1d5;
	} else if (boostType == 2) {
		gain = &gSpan::gain2;
		gainbound = &gSpan::gainbound2;
	}

	/* Precompute the overall label*weight sum-balance.
	 */
	boostWeightSum = 0.0;
	for (unsigned int n = 0 ; n < boostY.size() ; ++n)
		boostWeightSum += boostWeights[n] * boostY[n];

	bestGraphs.clear ();
}

void gSpan::run_graphs (const mxArray* graphs, int nlhs, mxArray* plhs[],
		 unsigned int _minsup,
		 unsigned int _maxpat_min, unsigned int _maxpat_max,
		 bool _enc,
		 bool _where,
		 bool _directed)
{
	ID = 0;
	minsup = _minsup;
	maxpat_min = _maxpat_min;
	maxpat_max = _maxpat_max;
	enc = _enc;
	where = _where;
	directed = _directed;

	mex = true;
	mexOutputRoot = NULL;
	mexOutputLast = NULL;

	unsigned int gcount = mxGetN (graphs);
	read (graphs);
	run_intern ();

	/* For the boosting mode, we only now convert the collected graphs to mex
	 * structures.
	 */
	if (boost) {
		for (unsigned int n = 0 ; n < bestGraphs.size() ; ++n)
			mexAppendGraph (&bestGraphs[n], 0, bestGraphsY[n], bestGraphsCounts[n]);
		// TODO: also output gains (bestGraphsGain)
	}

	/* Produce the desired mex structures.
	 */
	unsigned int count = 0;
	for (mexOutputG* cur = mexOutputRoot ; cur != NULL ; cur = cur->next)
		count += 1;

	/* So we have 'count' subgraphs now, allocate and fill in
	 *
	 * plhs[0]: subg
	 * plhs[1]: (count or ybase)
	 * plhs[2]: GY
	 */
	plhs[0] = mxCreateCellMatrix (1, count);
	assert (plhs[0] != NULL);
	//graphsP = (mxArray**) mxGetPr (plhs[0]);
	mxArray* graphsP = plhs[0];

	uint32_T* countP = NULL;
	double* subgY = NULL;
	double* countsGY = NULL;
	if (nlhs >= 2 && boost == false) {
		plhs[1] = mxCreateNumericMatrix (1, count, mxUINT32_CLASS, mxREAL);
		assert (plhs[1] != NULL);
		countP = (uint32_T *) mxGetPr (plhs[1]);
		assert (countP != NULL);
	} else if (nlhs >= 2 && boost) {
		plhs[1] = mxCreateNumericMatrix (1, count, mxDOUBLE_CLASS, mxREAL);
		subgY = (double *) mxGetPr (plhs[1]);
	}
	if (nlhs >= 3) {
		plhs[2] = mxCreateNumericMatrix (gcount, count, mxDOUBLE_CLASS, mxREAL);
		countsGY = (double*) mxGetPr (plhs[2]);
	}

	unsigned int idx = 0;
	for (mexOutputG* cur = mexOutputRoot ; cur != NULL ;
		cur = cur->next, ++idx)
	{
		//graphsP[idx] = cur->graph;
		mxSetCell (graphsP, idx, cur->graph);
		if (nlhs >= 2 && boost == false)
			countP[idx] = (unsigned int) cur->par;
		else if (nlhs >= 2 && boost)
			subgY[idx] = bestGraphsY[idx];

		/* Add individual graph counts.
		 */
		if (nlhs >= 3) {
			for (std::map<unsigned int, unsigned int>::iterator piter =
				cur->counts.begin() ; piter != cur->counts.end() ; ++piter)
			{
				countsGY[idx*gcount+(*piter).first] = (*piter).second;
			}
		}
	}

	/* Deallocate
	 */
	mexOutputG* cur = mexOutputRoot;
	while (cur != NULL) {
		mexOutputG* ccur = cur;
		cur = cur->next;

		delete ccur;
	}
}

void gSpan::run (std::istream &is, std::ostream &_os,
		 unsigned int _minsup,
		 unsigned int _maxpat_min, unsigned int _maxpat_max,
		 bool _enc,
		 bool _where,
		 bool _directed)
{
	os = &_os;
	ID = 0;
	minsup = _minsup;
	maxpat_min = _maxpat_min;
	maxpat_max = _maxpat_max;
	enc = _enc;
	where = _where;
	directed = _directed;
	mex = false;
	boost = false;

	read (is);
	run_intern ();
}

void gSpan::run_intern (void)
{
	/* In case 1 node subgraphs should also be mined for, do this as
	 * preprocessing step.
	 */
	if (maxpat_min <= 1) {
		/* Do single node handling, as the normal gspan DFS code based processing
		 * cannot find subgraphs of size |subg|==1.  Hence, we find frequent node
		 * labels explicitly.
		 */
		for (unsigned int id = 0; id < TRANS.size(); ++id) {
			for (unsigned int nid = 0 ; nid < TRANS[id].size() ; ++nid) {
				if (singleVertex[id][TRANS[id][nid].label] == 0) {
					// number of graphs it appears in
					singleVertexLabel[TRANS[id][nid].label] += 1;
				}

				singleVertex[id][TRANS[id][nid].label] += 1;
			}
		}
		/* All minimum support node labels are frequent 'subgraphs'.
		 * singleVertexLabel[nodelabel] gives the number of graphs it appears
		 * in.
		 *
		 * 1/1.5-class case: All nodelabels that do not appear at all have a
		 *    gain of zero, hence we do not need to consider them.
		 *
		 * 2-class case: Not appearing nodelabels are counted negatively.
		 */
		for (std::map<unsigned int, unsigned int>::iterator it =
			singleVertexLabel.begin () ; it != singleVertexLabel.end () ; ++it)
		{
			if ((*it).second < minsup)
				continue;

			unsigned int frequent_label = (*it).first;

			/* Found a frequent node label, report it.
			 */
			Graph g(directed);
			g.resize (1);
			g[0].label = frequent_label;

			/* [graph_id] = count for current substructure
			 */
			std::vector<unsigned int> counts (TRANS.size ());
			for (std::map<unsigned int, std::map<unsigned int, unsigned int> >::iterator it2 =
				singleVertex.begin () ; it2 != singleVertex.end () ; ++it2)
			{
				counts[(*it2).first] = (*it2).second[frequent_label];
			}

			if (boost) {
				/* Calculate gain and yval.  Here we do not use the normal
				 * gain function as there is no Projected/DFS_CODE there yet.
				 * Hence we need to make a distinction between the 1/1.5-class
				 * and the 2-class case here.
				 */
				double gainm = 0.0;
				double gainm_pos = 0.0;
				double gainm_neg = 0.0;

				for (unsigned int cid = 0 ; cid < counts.size () ; ++cid) {
					if (boostType == 1) {
						// Only consider the positive instances
						if (counts[cid] == 0)
							continue;

						gainm += boostY[cid]*boostWeights[cid];	// *1.0 (Y)
					} else if (boostType == 2) {
						double addfactor = 1.0;

						if (counts[cid] == 0)
							addfactor = -1.0;	// negation: pattern does not exist

						gainm_pos += addfactor*boostY[cid]*boostWeights[cid];
						gainm_neg += -addfactor*boostY[cid]*boostWeights[cid];
					}
				}
				double yval = 1.0;

				if (boostType == 2) {
					if (gainm_pos >= gainm_neg) {
						gainm = gainm_pos;
					} else {
						gainm = gainm_neg;
						yval = -1.0;
					}
				}

//#ifdef	DEBUG
				mexPrintf ("   single node graph, node label %d, gain %lf, yval %lf\n",
					frequent_label, gainm, yval);
//#endif
				if (gainm > boostTau) {
					/* Copy it into vector form
					 */
					std::map<unsigned int, unsigned int> gycounts;
					for (unsigned int n = 0 ; n < counts.size () ; ++n)
						gycounts[n] = counts[n];

					report_boosting_inter (g, (*it).second, gainm, yval, gycounts);
				}
			} else {
				std::map<unsigned int, unsigned int> gycounts;
				for (unsigned int n = 0 ; n < counts.size () ; ++n)
					gycounts[n] = counts[n];

				report_single (g, gycounts);
			}
		}
	}

	EdgeList edges;
	Projected_map3 root;

	for (unsigned int id = 0; id < TRANS.size(); ++id) {
		Graph &g = TRANS[id];
		for (unsigned int from = 0; from < g.size() ; ++from) {
			if (get_forward_root (g, g[from], edges)) {
				for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it)
					root[g[from].label][(*it)->elabel][g[(*it)->to].label].push (id, *it, 0);
			}
		}
	}

	for (Projected_iterator3 fromlabel = root.begin() ;
		fromlabel != root.end() ; ++fromlabel)
	{
		for (Projected_iterator2 elabel = fromlabel->second.begin() ;
			elabel != fromlabel->second.end() ; ++elabel)
		{
			for (Projected_iterator1 tolabel = elabel->second.begin();
				tolabel != elabel->second.end(); ++tolabel)
			{
				/* Build the initial two-node graph.  It will be grown
				 * recursively within project.
				 */
				DFS_CODE.push (0, 1, fromlabel->first, elabel->first, tolabel->first);
				project (tolabel->second);
				DFS_CODE.pop ();
			}
		}
	}
}

}
