/*
    $Id: ismin.cpp,v 1.5 2004/05/21 05:50:13 taku-ku Exp $;
 
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

namespace GSPAN {

bool gSpan::is_min ()
{
	if (DFS_CODE.size() == 1)
		return (true);

	DFS_CODE.toGraph (GRAPH_IS_MIN);
	DFS_CODE_IS_MIN.clear ();

	Projected_map3 root;
	EdgeList           edges;

	for (unsigned int from = 0; from < GRAPH_IS_MIN.size() ; ++from)
		if (get_forward_root (GRAPH_IS_MIN, GRAPH_IS_MIN[from], edges))
			for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it)
				root[GRAPH_IS_MIN[from].label][(*it)->elabel][GRAPH_IS_MIN[(*it)->to].label].push (0, *it, 0);

	Projected_iterator3 fromlabel = root.begin();
	Projected_iterator2 elabel    = fromlabel->second.begin();
	Projected_iterator1 tolabel   = elabel->second.begin();

	DFS_CODE_IS_MIN.push (0, 1, fromlabel->first, elabel->first, tolabel->first);

	return (project_is_min (tolabel->second));
}

bool gSpan::project_is_min (Projected &projected)
{
	const RMPath& rmpath = DFS_CODE_IS_MIN.buildRMPath ();
	int minlabel         = DFS_CODE_IS_MIN[0].fromlabel;
	int maxtoc           = DFS_CODE_IS_MIN[rmpath[0]].to;

	{
		Projected_map1 root;
		bool flg = false;
		int newto = 0;

		for (int i = rmpath.size()-1; ! flg  && i >= 1; --i) {
			for (unsigned int n = 0; n < projected.size(); ++n) {
				PDFS *cur = &projected[n];
				History history (GRAPH_IS_MIN, cur);
				Edge *e = get_backward (GRAPH_IS_MIN, history[rmpath[i]], history[rmpath[0]], history);
				if (e) {
					root[e->elabel].push (0, e, cur);
					newto = DFS_CODE_IS_MIN[rmpath[i]].from;
					flg = true;
				}
			}
		}

		if (flg) {
			Projected_iterator1 elabel = root.begin();
			DFS_CODE_IS_MIN.push (maxtoc, newto, -1, elabel->first, -1);
			if (DFS_CODE[DFS_CODE_IS_MIN.size()-1] != DFS_CODE_IS_MIN [DFS_CODE_IS_MIN.size()-1]) return false;
			return project_is_min (elabel->second);
		}
	}

	{
		bool flg = false;
		int newfrom = 0;
		Projected_map2 root;
		EdgeList edges;

		for (unsigned int n = 0; n < projected.size(); ++n) {
			PDFS *cur = &projected[n];
			History history (GRAPH_IS_MIN, cur);
			if (get_forward_pure (GRAPH_IS_MIN, history[rmpath[0]], minlabel, history, edges)) {
				flg = true;
				newfrom = maxtoc;
				for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it)
					root[(*it)->elabel][GRAPH_IS_MIN[(*it)->to].label].push (0, *it, cur);
			}
		}

		for (int i = 0; ! flg && i < (int)rmpath.size(); ++i) {
			for (unsigned int n = 0; n < projected.size(); ++n) {
				PDFS *cur = &projected[n];
				History history (GRAPH_IS_MIN, cur);
				if (get_forward_rmpath (GRAPH_IS_MIN, history[rmpath[i]], minlabel, history, edges)) {
					flg = true;
					newfrom = DFS_CODE_IS_MIN[rmpath[i]].from;
					for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it)
						root[(*it)->elabel][GRAPH_IS_MIN[(*it)->to].label].push (0, *it, cur);
				}
			}
		}

		if (flg) {
			Projected_iterator2 elabel  = root.begin();
			Projected_iterator1 tolabel = elabel->second.begin();
			DFS_CODE_IS_MIN.push (newfrom, maxtoc + 1, -1, elabel->first, tolabel->first);
			if (DFS_CODE[DFS_CODE_IS_MIN.size()-1] != DFS_CODE_IS_MIN [DFS_CODE_IS_MIN.size()-1]) return false;
			return project_is_min (tolabel->second);
		}
	}

	return true;
}
}
