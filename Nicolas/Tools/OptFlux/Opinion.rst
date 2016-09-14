############
OptFlux Test
############

* hairball maps! due to automatical layout (force directed layout, physics simulation where nodes repel each other)

* "Each node of the graph can have its 2D coordinates defined, but they are not mandatory. For the nodes that do not have a position, an automatic layout provided by Prefuse is applied, called Force Directed Layout, that positions graph elements based on a physics simulation of interacting forces; by default, nodes repel each other, edges act as springs, and drag forces are applied."

Supported Input:
================

"Currently, this layer supports two formats: CellDesigner SBML and XGMML. CellDesigner is a well-known tool for drawing cellular networks, based on the process diagram paradigm, with a graphical notation system proposed by Kitano. These are stored using an extension of the Systems Biology Markup Language (SBML). On the other hand, XGMML (eXtensible Graph Markup and Modeling Lan-guage) [2] is an extension to GML (Graph Markup Language), and it is used for graph descriptionNetwork Visualization Tools to Enhance Metabolic Engineering Platforms using XML tags to describe nodes and edges of a graph, being supported for instance by Cytoscape."

*****************
How do they draw?
*****************

* "Each metabolite has information on its identifier and coordinates. Therefore, there are three types of nodes: reaction, metabolite and information nodes. The visual representation of the network is handled by the Prefuse library, a graph visualization framework written in Java." 

*******
Opinion
*******

* we want defined lauout not automatical, force directed