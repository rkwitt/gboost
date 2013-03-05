/*
    $Id: dfs.cpp,v 1.3 2004/05/21 05:50:13 taku-ku Exp $;
 
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
#include <cstring>
#include <string>
#include <iterator>
#include <set>

namespace GSPAN {

/* Build a DFS code from a given graph.
 */
void
DFSCode::fromGraph (Graph &g)
{
	clear ();

	EdgeList edges;
	for (unsigned int from = 0 ; from < g.size () ; ++from) {
		if (get_forward_root (g, g[from], edges) == false)
			continue;

		for (EdgeList::iterator it = edges.begin () ; it != edges.end () ; ++it)
			push (from, (*it)->to, g[(*it)->from].label, (*it)->elabel, g[(*it)->to].label);
	}
}

bool DFSCode::toGraph (Graph &g)
{
	g.clear ();

	for (DFSCode::iterator it = begin(); it != end(); ++it) {
		g.resize (std::max (it->from, it->to) + 1);

		if (it->fromlabel != -1)
			g[it->from].label = it->fromlabel;
		if (it->tolabel != -1)
			g[it->to].label = it->tolabel;

		g[it->from].push (it->from, it->to, it->elabel);
		if (g.directed == false)
			g[it->to].push (it->to, it->from, it->elabel);
	}

	g.buildEdge ();

	return (true);
}

unsigned int
DFSCode::nodeCount (void)
{
	unsigned int nodecount = 0;

	for (DFSCode::iterator it = begin() ; it != end() ; ++it)
		nodecount = std::max (nodecount, (unsigned int) (std::max (it->from, it->to) + 1));

	return (nodecount);
}


std::ostream &DFSCode::write (std::ostream &os)
{
	if (size() == 0) return os;

	os << "(" << (*this)[0].fromlabel << ") " << (*this)[0].elabel << " (0f" << (*this)[0].tolabel << ")";

	for (unsigned int i = 1; i < size(); ++i) {
		if ((*this)[i].from < (*this)[i].to) {
			os << " " << (*this)[i].elabel << " (" << (*this)[i].from << "f" << (*this)[i].tolabel << ")";
		} else {
			os << " " << (*this)[i].elabel << " (b" << (*this)[i].to << ")";
		}
	}

	return os;
}
}

