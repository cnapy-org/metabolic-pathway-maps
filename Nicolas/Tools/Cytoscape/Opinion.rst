..	include:: <isonum.txt>

##############
Cytoscape Test
##############

* "Cytoscape supports a variety of automated network layout algorithms, including spring-embedded layout, hierarchical layout, and circular layout (...) It models a mechanical system in which edges of the graph correspond to springs, creating an attractive force between nodes that are far apart, and a repulsive force between nodes that are close together."

* "Whereas layout determines the location of the nodes and edges in the window, an attribute-to-visual mapping allows data attributes to control the appearance of their associated nodes and edges. Cytoscape supports a wide variety of visual properties, such as node color, shape, and size; node border color and thickness; and edge color, thickness, and style; a data attribute is mapped to a visual property using either a lookup table or interpolation, depending on whether the attribute is discrete valued or continuous."

* can only load SBML files with additional plugin (http://apps.cytoscape.org/apps/cysbml)
	|rarr| but: no longer compatible with new Cytoscape version 3.4 (plugins where changed to apps, exact error message: "Jar is missing value for entry Cytoscape-App in its manifest file.")