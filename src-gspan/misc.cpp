/*
    $Id: misc.cpp,v 1.6 2004/05/21 05:50:13 taku-ku Exp $;
 
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
#include <assert.h>

namespace GSPAN {

const RMPath &DFSCode::buildRMPath ()
{
	rmpath.clear ();

	int old_from = -1;

	for (int i = size() - 1 ; i >= 0 ; --i) {
		if ((*this)[i].from < (*this)[i].to && // forward
				(rmpath.empty() || old_from == (*this)[i].to))
		{
			rmpath.push_back (i);
			old_from = (*this)[i].from;
		}
	}

	return rmpath;
}

void History::build (Graph &graph, PDFS *e)
{
	// first build history
	clear ();
	edge.clear ();
	edge.resize (graph.edge_size());
	vertex.clear ();
	vertex.resize (graph.size());

	if (e) {
		push_back (e->edge);
		edge[e->edge->id] = vertex[e->edge->from] = vertex[e->edge->to] = 1;

		for (PDFS *p = e->prev ; p ; p = p->prev) {
			push_back (p->edge);	// this line eats 8% of overall instructions(!)
			edge[p->edge->id] = vertex[p->edge->from] = vertex[p->edge->to] = 1;
		}
		std::reverse (begin(), end());
	}
}

/* get_forward_pure ()
   e1 (from1, elabel1, to1)
   from から繋がる edge e2(from2, elabel2, to2) を返す.

   minlabel <= elabel2,
   (elabel1 < elabel2 ||
   (elabel == elabel2 && tolabel1 < tolabel2) の条件をみたす. 
   (elabel1, to1) のほうが先に探索されるべき
   また, いままで見た vertex には逝かない (backward のやくめ)
*/
bool get_forward_rmpath (Graph &graph, Edge *e, int minlabel, History& history, EdgeList &result)
{
	result.clear ();
	assert (e->to >= 0 && e->to < graph.size ());
	assert (e->from >= 0 && e->from < graph.size ());
	int tolabel = graph[e->to].label;

	for (Vertex::edge_iterator it = graph[e->from].edge.begin() ;
		it != graph[e->from].edge.end() ; ++it)
	{
		int tolabel2 = graph[it->to].label;
		if (e->to == it->to || minlabel > tolabel2 || history.hasVertex (it->to))
			continue;

		if (e->elabel < it->elabel || (e->elabel == it->elabel && tolabel <= tolabel2))
			result.push_back (&(*it));
	}

	return (! result.empty());
}

/* get_forward_pure ()
   e (from, elabel, to)
   to から繋がる edge を返す
   ただし, minlabel より大きいものにしかいかない (DFSの制約)
   また, いままで見た vertex には逝かない (backward のやくめ)
 */
bool get_forward_pure (Graph &graph, Edge *e, int minlabel, History& history, EdgeList &result)
{
	result.clear ();

	assert (e->to >= 0 && e->to < graph.size ());

	/* Walk all edges leaving from vertex e->to.
	 */
	for (Vertex::edge_iterator it = graph[e->to].edge.begin() ;
		it != graph[e->to].edge.end() ; ++it)
	{
		/* -e-> [e->to] -it-> [it->to]
		 */
		assert (it->to >= 0 && it->to < graph.size ());
		if (minlabel > graph[it->to].label || history.hasVertex (it->to))
			continue;

		result.push_back (&(*it));
	}

	return (! result.empty());
}

/* graph の vertex からはえる edge を探す
   ただし, fromlabel <= tolabel の性質を満たす.
*/
bool get_forward_root (Graph &g, Vertex &v, EdgeList &result)
{
	result.clear ();
	for (Vertex::edge_iterator it = v.edge.begin(); it != v.edge.end(); ++it) {
		assert (it->to >= 0 && it->to < g.size ());
		if (v.label <= g[it->to].label)
			result.push_back (&(*it));
	}

	return (! result.empty());
}

/* get_backward (graph, e1, e2, history);
   e1 (from1, elabel1, to1)
   e2 (from2, elabel2, to2)
   to2 -> from1 に繋がるかどうかしらべる.

   (elabel1 < elabel2 ||
   (elabel == elabel2 && tolabel1 < tolabel2) の条件をみたす. (elabel1, to1) のほうが先に探索されるべき
 */
Edge *get_backward (Graph &graph, Edge* e1, Edge* e2, History& history)
{
	if (e1 == e2)
		return 0;

	assert (e1->from >= 0 && e1->from < graph.size ());
	assert (e1->to >= 0 && e1->to < graph.size ());
	assert (e2->to >= 0 && e2->to < graph.size ());

	for (Vertex::edge_iterator it = graph[e2->to].edge.begin() ;
		it != graph[e2->to].edge.end() ; ++it)
	{
		if (history.hasEdge (it->id))
			continue;

		if ( (it->to == e1->from) &&
            (
              (e1->elabel < it->elabel) ||
                (
                  (e1->elabel == it->elabel) &&
                  (graph[e1->to].label <= graph[e2->to].label)
                )
            )
       )
		{
			return &(*it);
		}
	}

	return 0;
}
}

