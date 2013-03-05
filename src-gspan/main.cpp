/*
    $Id: main.cpp,v 1.4 2004/05/21 05:50:13 taku-ku Exp $;
 
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
#include <unistd.h>

#define OPT " [-m minsup] [-d] [-e] [-w] "

void usage (void)
{
	std::cout << "gspan implementation by Taku Kudo" << std::endl;
	std::cout << std::endl;
	std::cout << "usage: gspan [-m minsup] [-D] [-e] [-w] [-L maxpat]" << std::endl;
	std::cout << std::endl;
	std::cout << "options" << std::endl;
	std::cout << "  -h, show this usage help" << std::endl;
	std::cout << "  -m minsup, set the minimum support (absolute count)" << std::endl;
	std::cout << "  -D, use directed edges, default: undirected" << std::endl;
	std::cout << "  -e, output substructures in encoded form (?)" << std::endl;
	std::cout << "  -w, where (?)" << std::endl;
	std::cout << "  -L maxpat, the maximum number of outputted substructures" << std::endl;
	std::cout << "  -n minnodes, the minimum number of nodes in substructes (default: 0)" << std::endl;
	std::cout << std::endl;

	std::cout << "The graphs are read from stdin, and have to be in this format:" << std::endl;
	std::cout << "t" << std::endl;
	std::cout << "v <vertex-index> <vertex-label>" << std::endl;
	std::cout << "..." << std::endl;
	std::cout << "e <edge-from> <edge-to> <edge-label>" << std::endl;
	std::cout << "..." << std::endl;
	std::cout << "<empty line>" << std::endl;
	std::cout << std::endl;

	std::cout << "Indices start at zero, labels are arbitrary unsigned integers." << std::endl;
	std::cout << std::endl;
}

int main (int argc, char **argv)
{
	unsigned int minsup = 1;
	unsigned int maxpat = 0xffffffff;
	unsigned int minnodes = 0;
	bool where = false;
	bool enc = false;
	bool directed = false;

	int opt;
	while ((opt = getopt(argc, argv, "edws::m:L:Dhn:")) != -1) {
		switch(opt) {
		case 's':
		case 'm':
			minsup = atoi (optarg);
			break;
		case 'n':
			minnodes = atoi (optarg);
			break;
		case 'L':
			maxpat = atoi (optarg);
			break;
		case 'd': // same as original gSpan
		case 'e':
			enc = true;
			break;
		case 'w':
			where = true;
			break;
		case 'D':
			directed = true;
			break;
		case 'h':
		default:
			usage ();
			return -1;
		}
	}

	GSPAN::gSpan gspan;
	gspan.run (std::cin, std::cout, minsup, maxpat, minnodes, enc, where, directed);
}
