
// TODO: fast matching, aborting after just one match (return true in visitor)

/* graphmatch.cpp - MEX interface to VFlib2 graph monomorphism matching
 *
 * Author: Sebastian Nowozin <sebastian.nowozin@tuebingen.mpg.de>
 * Date: 1st August 2006
 *
 * Usage:
 *
 * [count, matches] = graphmatch (subg1, g2, only_one, directed);
 *
 * Finds all exact occurences of subgraph subg1 in the larger graph g2.
 * Optionally stores all found matches.
 *
 ***
 * Input
 *
 * subg1, g2: structures defining graphs, of this layout:
 *    g.nodelabels: (n,1) discrete integer labels [L_1 ; L_2 ; ... ; L_n];
 *    g.edges: (m,2) edges, [from to] at each line:
 *       [e_1_{from} e_1_{to} edgelabel_1 ; ... ; e_m_{from} e_m_{to} edgelabel_m]
 *       The node indices go from 1 to n.
 *
 *    subg1 needs to be smaller or equally sized as g2, nodelabels and edges
 *    need to be uint32's.
 *
 * (optional) only_one: If given and non-zero, then only test whether the
 *    substructure appears once.  This aborts the search early in case a graph
 *    is found.
 *
 * (optional) directed: If zero, undirected graphs are matched, if non zero,
 *    the graphs are assumed to be directed graphs (default: 1).
 ***
 * Output
 *
 * count (double): The number of times the subg1 graph appears in g2.
 *
 * (optional) matches: (count,K) uint32 matrix, where count is the number of
 *    matches, one match per row; and K is the number of nodes in subg1.  In
 *    each row the node indices of g2 are stored such that the i'th column
 *    corresponds to the i'th node in subg1.  If there is a 5 as the 3rd
 *    element, then node number 5 from g2 is matched with node number 3 from
 *    subg1.
 */

#include <assert.h>

/* VFlib2 include files
 */
#include <argraph.h>
#include <argedit.h>
#include <vf2_mono_state.h>
#include <match.h>
#include <vector>
#include <inttypes.h>

/* Matlab MEX interface include files
 */
#include <mex.h>

/* Maximum number of labels at each edge.
 */
#define	ELABEL_MAX	256
#define	ELABEL_END	0xffffffff

/*** Local functions and classes
 */
typedef struct mv_chain {
	struct mv_chain* next;
	mxArray* row;
} mv_chain;

typedef struct {
	unsigned int count;	/* Count of the number of matches */

	int do_collection;		/* Collect all positive match permutations */
	int do_check_only;		/* Check only if at least one substructure exists */
	mv_chain* last;			/* Pointer to match permutations */
} mv_usr;

// If directed is 0: graph is undirected, directed is 1: directed graph
// (default)
static ARGraph<unsigned int, unsigned int> *
convert_mat_to_graph (const mxArray* par, int directed);

static bool
matched_visitor (int n, node_id ni1[], node_id ni2[], void *usr);

class UIntComparator
	: public AttrComparator
{
	public:
		virtual bool compatible (void* pa, void* pb)
		{
            // Cast to unsigned int ptr
            uintptr_t *a_ptr = (uintptr_t *)&pa;
            uintptr_t *b_ptr = (uintptr_t *)&pb;

            // Get value
            unsigned int a_val = *a_ptr;
            unsigned int b_val = *b_ptr;
            
            return (a_val == b_val);
		}
};

class UIntPtrComparator
	: public AttrComparator
{
	public:
		virtual bool compatible (void* pa, void* pb)
		{
			unsigned int* a = (unsigned int*) pa;
			unsigned int* b = (unsigned int*) pb;

			for (unsigned int n1 = 0 ; a[n1] != ELABEL_END ; ++n1)
				for (unsigned int n2 = 0 ; b[n2] != ELABEL_END ; ++n2)
					if (a[n1] == b[n2])
						return (true);

			return (false);
		}
};


static std::vector<unsigned int*> elist_ptr;


/* [count, associations] = graphmatch (subg1, g2);
 */
void
mexFunction (int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	elist_ptr.clear ();

	if (nlhs > 2) {
		mexPrintf ("Too many left hand side parameters.\n");

		return;
	}

	if (nrhs < 2 || nrhs > 4) {
		mexPrintf ("Wrong number of right hand side parameters.\n");

		return;
	}

	if (mxIsStruct (prhs[0]) == 0 || mxIsStruct (prhs[1]) == 0) {
		mexPrintf ("subg1 and g2 must be graph structures.\n");

		return;
	}

	unsigned int directed = 1;
	if (nrhs >= 4)
		directed = (unsigned int) mxGetScalar (prhs[3]);

	/* Convert Matlab structures to VFlib graphs.
	 */
	ARGraph<unsigned int, unsigned int>* subg1 = convert_mat_to_graph (prhs[0],
		directed);
	ARGraph<unsigned int, unsigned int>* g2 = convert_mat_to_graph (prhs[1],
		directed);
	if (subg1 == NULL || g2 == NULL) {
		mexPrintf ("Error parsing subg1 or g2 structure.\n");

		return;
	}

	/* Do the graph matching.
	 */
	subg1->SetNodeComparator (new UIntComparator ());
	subg1->SetEdgeComparator (new UIntPtrComparator ());
	g2->SetNodeComparator (new UIntComparator ());
	g2->SetEdgeComparator (new UIntPtrComparator ());

	/* Do the matching.
	 */
	//VF2SubState s0(subg1, g2);	// graph-subgraph isomorphism
	VF2MonoState s0(subg1, g2);	// monomorphism
	mv_usr mvs = { 0, (nlhs >= 2) ? 1 : 0,
		(nrhs > 2 && mxGetScalar (prhs[1]) > 0.0) ? 1 : 0, NULL };
	match (&s0, matched_visitor, &mvs);

	delete subg1;
	delete g2;

	plhs[0] = mxCreateDoubleScalar (mvs.count);

	if (nlhs >= 2 && mvs.count == 0) {
		plhs[1] = mxCreateNumericMatrix (0, 0, mxUINT32_CLASS, 0);
	} else if (nlhs >= 2) {
		unsigned int cols = mxGetN (mvs.last->row);
		plhs[1] = mxCreateNumericMatrix (mvs.count, cols, mxUINT32_CLASS, 0);
		uint32_T* ass = (uint32_T*) mxGetPr (plhs[1]);

		int drow = mvs.count - 1;
		for (mv_chain* cur = mvs.last ; cur != NULL ; ) {
			uint32_T* row = (uint32_T*) mxGetPr (cur->row);

			/* Copy from each single row into the large matrix, reverse row
			 * order.
			 */
			for (int i = 0 ; i < cols ; ++i)
				ass[i*mvs.count+drow] = row[i];

			mv_chain* tbd = cur;
			cur = cur->next;

			mxDestroyArray (tbd->row);
			free (tbd);

			drow -= 1;
		}
	}

	for (unsigned int n = 0 ; n < elist_ptr.size () ; ++n)
		free (elist_ptr[n]);
	elist_ptr.clear ();
}


/* Produce an attributed relational graph from a Matlab structure.
 */
static ARGraph<unsigned int, unsigned int> *
convert_mat_to_graph (const mxArray* par, int directed)
{
	mxArray* nodelabels = mxGetField (par, 0, "nodelabels");
	mxArray* edges = mxGetField (par, 0, "edges");

	if (nodelabels == NULL || edges == NULL) {
		mexPrintf ("Graph structures are missing the \"nodelabels\" "
			"or \"edges\" fields.\n");

		return (NULL);
	}

	if (mxIsUint32 (nodelabels) == 0 || mxIsUint32 (edges) == 0) {
		mexPrintf ("\"nodelabels\" or \"edges\" not of type Uint32.\n");

		return (NULL);
	}

	/* Create graph
	 */
	ARGEdit ed;
	uint32_T* nmat = (uint32_T*) mxGetPr (nodelabels);
	unsigned int nodes = mxGetM (nodelabels);
	for (unsigned int i = 0 ; i < nodes ; ++i) {
		ed.InsertNode ((void *) nmat[i]);
#ifdef	DEBUG
		mexPrintf ("inserting node %d with label %d\n", i, nmat[i]);
#endif
	}

	/* emat is column organised.
	 */
	uint32_T* emat = (uint32_T*) mxGetPr (edges);
	unsigned int edgecount = mxGetM (edges);

	for (unsigned int i = 0 ; i < edgecount ; ++i) {
		unsigned int from = emat[i+0*edgecount];
		unsigned int to = emat[i+1*edgecount];
		unsigned int elabel = emat[i+2*edgecount];

		/* Collect all labels
		 */
		unsigned int elabel_set[ELABEL_MAX];
		unsigned int n;
		for (n = 0 ; (i + n) < edgecount ; ++n) {
			if (from == emat[(i+n)+0*edgecount] &&
				to == emat[(i+n)+1*edgecount])
			{
				elabel_set[n] = emat[(i+n)+2*edgecount];
			} else
				break;
		}
		i += n-1;

		elabel_set[n] = ELABEL_END;
		n += 1;
		unsigned int* elabel_set_p = (unsigned int*) calloc (n, sizeof (unsigned int));
		for (unsigned int k = 0 ; k < n ; ++k)
			elabel_set_p[k] = elabel_set[k];

		/* Append to the global junkyard list.
		 */
		elist_ptr.push_back (elabel_set_p);

#ifdef	DEBUG
		mexPrintf ("inserting edge from %d to %d, label %d\n",
			from, to, elabel);
#endif
		ed.InsertEdge (from-1, to-1, (void *) elabel_set_p);
		if (directed == 0)
			ed.InsertEdge (to-1, from-1, (void *) elabel_set_p);
	}

	ARGraph<unsigned int, unsigned int>* gr =
		new ARGraph<unsigned int, unsigned int> (&ed);

	return (gr);
}


/* Callback function: called whenever the subgraph is successfully matched.
 * Return value false means "continue searching", true means "stop".
 */
static bool
matched_visitor (int n, node_id ni1[], node_id ni2[], void *usr)
{
	mv_usr* mv = (mv_usr*) usr;

#ifdef DEBUG
	for (int i = 0 ; i < n ; ++i)
		mexPrintf ("%d ", ni1[i]);
	mexPrintf ("\n");
	for (int i = 0 ; i < n ; ++i)
		mexPrintf ("%d ", ni2[i]);
	mexPrintf ("\n");
#endif

	mv->count += 1;
	if (mv->do_collection == 0)
		return (mv->do_check_only ? true : false);

	/* Collect
	 */
	mv_chain* c = (mv_chain*) calloc (1, sizeof (mv_chain));

	// XXX: Check that ni1 is always in the right order.
	for (int i = 0 ; i < (n-1) ; ++i) {
		assert (ni1[i] <= ni1[i+1]);
	}

	c->row = mxCreateNumericMatrix(1, n, mxUINT32_CLASS, 0);
	uint32_T* rowp = (uint32_T*) mxGetPr (c->row);
	for (int i = 0 ; i < n ; ++i)
		rowp[i] = ni2[i]+1;	// Fix index +1, conforming to Matlab convention

	/* Link.
	 */
	c->next = mv->last;
	mv->last = c;

#ifdef DEBUG
	printf ("matched\n");
#endif
	return (mv->do_check_only ? true : false);
}


