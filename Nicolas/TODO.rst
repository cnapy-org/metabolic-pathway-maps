2016-09-22
----------

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

Brett: install dependecies for running script (Wind64 Anaconda)

 - pip install geojson tarjan libsbgnpy tulip-python


2016-11-21
----------

 - Brett: look into tarjan algorithm
 - Find tool for viewing SBGN to check output
 - SBML+Layout can be checked here: http://sysbioapps.dyndns.org/Layout

Set specific metabolite nodes positions relative to each other and fix their positions
Layout groups of reaction nodes around a fixed position metabolite node (in principle reaction nodes)

  - pathway1 = [S, R1, A, R2, B, R3, C]
  - pathway2 = [C, R4, D, R5, E, R6, F]
  - Put A at position Xp and layout pathway1
  - Put E at position Yp and layout pathway2
  - alias all non-connecting metabolites (in this case everything except C)

