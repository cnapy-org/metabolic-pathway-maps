2016-09-22
----------

General:
========

* tulip meta node functionality to collapse pathways/organelles

::

	tlp.Graph.createMetaNode(nodeSet[, multiEdges=True, delAllEdge=True])

	Closes a set of existing nodes into a metanode and returns it. Edges from nodes in the set to other nodes are replaced with edges from the metanode to the other nodes.

* click to highlight all edges connected to a reaction node
* alias a single connection
* have the option to export tulip graph as svg

Daily:
======
 Visualisation::

 - Mimoza: what creates the layout information in the two e2lo n2lo dictionaries in the first place or
   in other words,  what generates the layout.

   - Update 26.09.2016 on surfdrive


Python::

 - Start reading through Part One of "Learning Python, Fourth Edition (2009).pdf" in the literature folder.
   Skip the bits on IDLE and frozen binaries you should be using IPython as an interactive Python shell,
   see for example: http://ipython.readthedocs.io/en/stable/interactive/tutorial.html

   - Part One done (26.09.2016)


2016-09-28
----------

Mimoza::

 - Investigate transformation of the Tulip root (graph) object by "void" functions
 - Start working on script to load test SBML model (supplied by Brett) and extract/generate Tuplip; graph

Python 3 and CBMPy compatability, add to module::

 from __future__ import division, print_function
 from __future__ import absolute_import
 from __future__ import unicode_literals

New tools to possibly look into::

 - http://www.cs.bilkent.edu.tr/~ivis/ChiSE_daily/sample-app/ (can't handle large files)
 - https://github.com/persistentsystems/D3-SCoMoS.JS/ (can't handle large files)


2016-10-13
----------

Mimoza based layout, create SBML+Layout and SBGN

Read introduction and example of SBML layout package

 - http://co.mbine.org/specifications/sbml.level-3.version-1.layout.version-1.release-1

Look at w3schools XML tutorial (namespaces, elements, attributes, ...)

 - http://www.w3schools.com/xml/xml_whatis.asp

Brett - adding layout to existing file

Learning Python::

 - Chpt 4, 5, ...

2016-10-21
----------

 - Python modules for XML parsing: xml...minidom ...elementtree,

 - Manipulate Tulip graph e.g. line colour, dashes
 - reduce edge crossing
 - networkx for graph structure analysis, find cycles, subgraphs
 - different layouts on subsystems

2016-10-21
----------

Brett: install dependecies for running script (Wind64 Anaconda).

 - pip install geojson tarjan libsbgnpy tulip-python

General things.

 - Find tool for viewing SBGN to check output
 - SBML+Layout can be checked here: http://sysbioapps.dyndns.org/Layout
 - Brett: look into tarjan algorithm

Inverstigate the posibility of a pseudo-fixed layout method. Set specific metabolite nodes positions relative to each other and fix their positions.
Layout groups of reaction nodes around a fixed position metabolite node (in principle reaction nodes). Pseudo-code::

  - pathway1 = [S, R1, A, R2, B, R3, C]
  - pathway2 = [C, R4, D, R5, E, R6, F]
  - Put A at position XYa and layout subgraph/pathway1
  - Put D at position XYb and layout subgraph/pathway2
  - alias all non-connecting metabolites (in this case everything except C)

2016-11-11
----------

Brett: I was playing around with the code and some larger models, here are some notes:

 - created the file `TestSBML2TulipGraph_brett.py` to test autoaliasing on highly connected networks see the `autoSplitConnectedNodes()` function.
 - Mimoza can't seem to handle multiple compartments that are not called "Cell" find out what is going on.
 - have option to display reaction id's on final graph so, reaction and/or species ids/names
 - consider how we will deal with species we don't want to display at all (H+, H20) this should be an optional input, my feeling is that they should be there for the layout but be made invisible on display
 - (future) isolate and group all transport reactions to one side
 - (future) profile the Mimoza/tulip code with a big model to find bottlenecks
 - colour specific reactions different colours ... .
 - new CONDA installer for CBMPy: "conda install -c bgoli -c SBMLTeam cbmpy"











2016-11-22
----------

 - Check function definitions for global variables
 - create import library: MVLayout.py with header
 - clearly mark/group model specific code
 - https://sourceforge.net/p/auber/code/HEAD/tree/tulip/doc/python/tulippluginsdocumentation.rst
 - look at subgraph edge direction
 - look at subgraph properties
 - Brett: create test Tulip graphs, test Lactis model



"""
MVLayout: layout module
=======================

MetaToolkit genome scale visualisation module (http://cbmpy.sourceforge.net)
Copyright (C) 2016 Nicolas Huber and Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

Author: Nicolas Huber
Contact email: n.a.f.huber@vu.nl

"""

def applyCircularOGDF(g, params={}):
	"""
	# params['minDistCircle'] = ...
	# params['minDistLevel'] = ...
	# params['minDistSibling'] = ...
	# params['minDistCC'] = ...
	# params['pageRatio'] = ...

	"""

	sg_params = tlp.getDefaultPluginParameters('Circular (OGDF)', g)
	print(sg_params)
	for k in params:
		if k in sg_params:
			sg_params[k] = params[k]
	success = g.applyLayoutAlgorithm('Circular (OGDF)', sg_params)
	print('applyCircularOGDF: {}'.format(success))
