#################################################
CBMViz Documentation: dictionary graph_registry #
#################################################

Description
-----------

The dictionary 'graph_registry' contains all generated Tulip graph objects and their respective property dictionaries. Since the property objects are specific for each graph object, the property dictionary is saved together with the graph object.

Creation
--------

It is initially created after the Mimoza function 'tulip_graph2geojson' returns the initial graph (inside layoutModel). 

::

	# create the graph object
	graph = tulip_graph2geojson.graph2geojson(params)
	# set the name of the graph object
	graph.setName('initial graph')
	# create graph registry dict
	graph_registry = {}
	graph_registry['initial graph'] = [graph, graph_props]

Structure
---------

For each key (graph name) the dictionary contains a list (graph object, graph properties).
The second entry in this list (graph properties) is a dictionary itself with property names as keys and property objects as values: 

::

	# how graph_props is set up
	{'viewBorderColor': <Color property viewBorderColor>, 'viewLabelPosition': <Integer property viewLabelPosition>,...}
	

Usage in the Code
-----------------

Everytime a new graph object is created (output graphs like 'output_whole', etc.), the dictionary is updated:

graph_registry.update({'graph name': [graph object, dict of graph properties]})

