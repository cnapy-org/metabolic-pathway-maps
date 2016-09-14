..	include:: <isonum.txt>

############
MetDraw Test
############

* needs open source software Graphviz (Eclipse Public License)
	|rarr| PATH Variable needs to be changed and GRAPHVIZ_DOT environment variable that points to location of .../.../dot.exe needs to be added manually

* with default parameters MetDraw produces a lot of "messy hairballs" (see .png)
	|rarr| change parameters as described at http://csbl.bitbucket.org/metdraw/#install to obtain better result?

* if an SBML file is used which is not included as tested, default example the program causes an errorS (BML file: e_coli_core.xml, http://bigg.ucsd.edu/models/e_coli_core)

::

	Traceback (most recent call last):
	File "metdraw.py", line 206, in <module>
		metdraw(**metdraw_args)
	File "metdraw.py", line 88, in metdraw
		pieces = sbml.parse_sbml_file(file=sbml_filename)
	File "C:\Users\Nicolas\Desktop\LabR1_Teusink_Amsterdam\TestExistingModelViz\MetDraw\metdraw-master\src\sbml.py", line 39, in parse_sbml_file
		for elem in findall(listOfCompartments,'compartment'):
	File "C:\Users\Nicolas\Desktop\LabR1_Teusink_Amsterdam\TestExistingModelViz\MetDraw\metdraw-master\src\sbml.py", line 21, in findall
		return [e for e in elem.getiterator() if tagmatch(e,tag)]
	AttributeError: 'NoneType' object has no attribute 'getiterator'

*****************
How do they draw?
*****************

* python based software; takes a number of parameters to control the drawing process from an imput SBML file: "(...)The final layout uses the widely used Graphviz (Gansner and North, 2000) software. MetDraw converts the reactions to a series of edges and nodes using the Graphviz DOT language. Unlike other graph-drawing programs for biological networks, MetDraw inserts invisible nodes and edges and uses multiple graph elements to create a final layout that more closely resembles a classical biochemical map.(...)" 

*******
Opinion
*******

* does not do exactly what we want however can generate some kind of .svg just from SBML file
	|rarr| look at source code on github and see if there is something useful
	
