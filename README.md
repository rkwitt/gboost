gboost - Matlab graph boosting package
======================================

Version 0.1.1, 11th July 2007 (+ minor adjustments)

**Original website:** http://www.nowozin.net/sebastian/gboost/


Contents
--------

* [Description](#description)
* [Authors](#authors)
* [Licence](#lic)
* [Installation](#install)
* [Documentation](#doc)
* [Demonstration](#demo)
* [References](#ref)


Description
-----------
<a name="description"/>

This package contains a Matlab interface to various libraries in order to
perform graph boosting [\[Kudo2004\]](#ref) and frequent subgraph mining [\[Yan2002\]](#ref).
Graph boosting learns a classification function on discrete-labeled undirected
connected graphs.  Frequent subgraph mining determines subgraphs with a given
minimum support.

Discrete-labeled undirected connected graphs are a suitable representation for
problems where samples can be represented by individual parts (vertices) and
their relationships (edges).  Graph boosting has been particularly successful
for the classification of chemical molecules.

While graph boosting can also be used for regression this package does *not*
implement a regressor.


Authors
-------
<a name="authors"/>

* Sebastian Nowozin <sebastian.nowozin@tuebingen.mpg.de>
  * Matlab wrappers, LPBoost, modifications to gSpan implementation
* Taku Kudo <taku@google.com>
 * C++ gSpan implementation
* Intelligent Systems and Artificial Vision Lab, SIVALab of the University of Naples ''Federico II''.
 * VFLib graph matching library


License
-------
<a name="lic"/>

The software is dual licensed under the GNU General Public License version 2
and the Mozilla Public License, version 1.1.  This means that you can choose
any of the two licenses.

The licenses are included as LICENSE.txt (GPL version 2) and MPL-1.1.txt
(Mozilla Public License, version 1.1).

```bash
GNU General Public License

Copyright (C) 2006 Sebastian Nowozin,
Copyright (C) 2004 Taku Kudo,
Copyright (C) 2001 Dipartimento di Informatica e Sistemistica, Universit
    degli studi di Napoli ``Federico II'',
All rights reserved.

This is free software with ABSOLUTELY NO WARRANTY.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; version 2 of the License.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place - Suite 330, Boston, MA 02111-1307, USA
```

```
Mozilla Public License

``The contents of this distribution are subject to the Mozilla Public License
Version 1.1 (the "License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
http://www.mozilla.org/MPL/

Software distributed under the License is distributed on an "AS IS" basis,
WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License for
the specific language governing rights and limitations under the License.

     The Original Code is "gboost Matlab toolbox".

     The Initial Developer of the Original Code is Sebastian Nowozin.
     Portions created by Taku Kudo are Copyright (C) 2007.
     Portions created by Dipartimento di Informatica e Sistemistica,
        Universit degli studi di Napoli ``Federico II'' are
        Copyright (C) 2001.
     All Rights Reserved.

Alternatively, the contents of this distribution may be used under the terms
of the GNU General Public License, version 2 (the  "GPL 2 License"), in which
case the provisions of GNU General Public License are applicable instead of
those above.  If you wish to allow use of your version of this file only under
the terms of the GNU General Public License and not to allow others to use
your version of this file under the MPL, indicate your decision by deleting
the provisions above and replace  them with the notice and other provisions
required by the GNU General Public License.  If you do not delete the
provisions above, a recipient may use your version of this file under either
the MPL or the GNU General Public License, version 2.''
```

Installation
============
<a name="install"/>

Prequisites
     - Linux, x86, either 32 bit or 64 bit
     - Matlab R14 (also tested with 2012b)
     - CVX Matlab Optimization Kit, available under the GNU General Public
       License at http://www.stanford.edu/~boyd/cvx/

The program has been tested on Debian GNU/Linux, x86 32 bit and 64 bit
architectures with Matlab R14SP2, SP3, R2006a, and most recently on 
Ubuntu 12.04 (64bit) running R2012a.

If you want to modify the source code, you have to modify Makefile.options in
the root directory to setup your Matlab path.  Afterwards, a simple "make" in
the root directory should build both the graph matching wrapper as well as the
gSpan graph mining.


Documentation
-------------
<a name="doc"/>

The source code is well documented, but here is a list of the most important
parts.

`gspan.m` is the Matlab side interface of Taku's gSpan code.  It can perform
both frequent subgraph mining as well as weighted subgraph mining.  The first
is useful for data mining purposes, while the second is used in graph
boosting.

`findhypothesis_graph.m` is the interface between LPBoost and the weighted graph
mining algorithm (gSpan).  The duty of `findhypothesis_graph` is to create
decision stumps which correspond to the most violated constraint in the LP
dual (column-generation).  In fact, you can use the included `lpboost1d5.m` with
any other decision stump, you only need to write a suitable `findhypothesis_*.m`
function.

`lpboost.m` is an implementation of 1-class, 2-class and "1.5-class"
nu-LPBoosting.  The 1-class and 2-class formulations are explained in
[\[Demiriz2002\]](#ref), the 1.5-class formulation learns a 1-class classifier but also
takes into account negative samples.

`graphmatch.mex*` is a wrapper around VFLib to perform subgraph-graph
isomorphism matching.  It can output all matches and is used for the testing
on unlabeled samples.  It can match both directed and undirected graphs.

`mexgspan.mex*` is the gSpan Matlab wrapper.

`rocscore.m` is a simple function calculating the ROC AUC and ROC EER score as
well as the ROC curve itself.

Demonstration
-------------
<a name="demo"/>

Start Matlab and go to the bin/ directory.  Running the example.m script will
guide you through the training of a graph boosting classifier for a small
molecule example set. Assuming that your CVX install is at `/Software/cvx/`,
and you checked out gboost at `/Software/gboost`, run

```matlab
cd '/Software/cxv'
cvx_setup
cd '/Software/gboost'
example
```

If you have any questions, please feel free to email the first author.

References
----------
<a name="ref"/>

[Demiriz2002], Ayhan Demiriz, Kristin P. Bennett and John Shawe-Taylor,
   "Linear Programming Boosting via Column Generation", 2002, Journal of
   Machine Learning, Vol. 46, pages 225-254,
   http://www.rpi.edu/~bennek/rev_mlj6.ps

[Kudo2004], Taku Kudo, Eisaku Maeda and Yuji Matsumoto,
   "An Application of Boosting to Graph Classification", NIPS 2004,
   http://books.nips.cc/papers/files/nips17/NIPS2004_0369.pdf

[Yan2002], Xifeng Yan and Jiawei Han,
   "gSpan: Graph-Based Substructure Pattern Mining", ICDM 2002,
   http://computer.org/proceedings/icdm/1754/17540721abs.htm"

