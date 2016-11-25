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


from tulip import tlp
import graph_properties
import resize


def identifyConnectedNodes(model, spec_split, reac_split, split_biomass=True):
    """
    Identifies highly connected species and reactions::

     - *model* a CBMPy model object
     - *spec_split* species connected to this or more reactions are aliased
     - *reac_split* reactions connected to this or more non-split species are aliased
     - *split_biomass* [default=True] alias the objective function reaction reagents

    """
    sp_list = []
    for s in model.species:
        if len(s.isReagentOf()) > spec_split:
            sp_list.append(s.getId())

    # add ubiquitous reactions to a group
    re_list = []
    for r in model.reactions:
        if len((set(r.getSpeciesIds()).difference(set(sp_list)))) > reac_split:
            re_list.append(r.getId())
    if split_biomass:
        for r in model.getActiveObjectiveReactionIds():
            sp_list += [s for s in model.getReaction(r).getSpeciesIds() if s not in sp_list]
            if r in re_list:
                re_list.remove(r)

    return sp_list, re_list


def identifyBoundarySpecies(sps_list, re_list, cmod):
    """
    Identifies boundary species of given pathway by:
    1.) getting the nodes that have at least one reaction outside the pathway -> level0_nodes (still themselves inside
    the pathway - we don't want them!)
    2.) getting from the level0_nodes to the next connections (level1_nodes)
    => these are the first nodes outside the pathway => boundary nodes!
    ::

    - *sps_list* list containing all species of a pathway
    - *re_list* list containing all reactions of a pathway
    - *cmod* cbmpy model of the input sbml file

    """
    level0_sps = []
    level1_sps = []
    for sps in sps_list:
        sps = cmod.getSpecies(sps)
        for re in sps.isReagentOf():
            if re not in re_list and sps.getId() not in level0_sps:
                level0_sps.append(sps.getId())
                re = cmod.getReaction(re)
                # create a list of tuples (species, reaction): pair of pathway bound sps with connecting reaction
                level1_sps += [(s, re.getId()) for s in re.getSpeciesIds() if s not in level0_sps and s not in sps_list]

    bound_sps = level1_sps
    return list(set(bound_sps))


def create_props(graph):
    """
    Create tulip graph properties::

    - *graph* a tulip graph

    """
    graph.getStringProperty(graph_properties.ANCESTOR_TERM)
    graph.getStringProperty(graph_properties.ANCESTOR_ID)
    graph.getStringProperty(graph_properties.ANCESTOR_NAME)
    graph.getStringProperty(graph_properties.TERM)
    graph.getStringProperty(graph_properties.COMPARTMENT_ID)
    graph.getStringProperty(graph_properties.ID)
    graph.getStringProperty(graph_properties.NAME)
    graph.getIntegerProperty(graph_properties.TYPE)
    graph.getBooleanProperty(graph_properties.REVERSIBLE)
    graph.getDoubleProperty(graph_properties.STOICHIOMETRY)
    graph.getBooleanProperty(graph_properties.TRANSPORT)
    graph.getBooleanProperty(graph_properties.UBIQUITOUS)
    graph.getLayoutProperty(graph_properties.VIEW_LAYOUT)
    graph.getIntegerProperty(graph_properties.VIEW_SHAPE)
    graph.getSizeProperty(graph_properties.VIEW_SIZE)
    graph.getBooleanProperty(graph_properties.FAKE)
    graph.getStringProperty(graph_properties.CLONE_ID)
    graph.getStringProperty(graph_properties.FORMULA)
    graph.getStringVectorProperty(graph_properties.RELATED_COMPARTMENT_IDS)
    graph.getStringProperty(graph_properties.EDGE_DIRECTIONALITY)


def markUbiquitousNodes(g, ub_list):
    """
    Flag a list of nodes in a graph to be marked as ubiquitous (for aliasing)::

    - *g* a tulip/mimoza graph
    - *ub_list* a list of nodes

    """
    # set graph property for ub nodes
    #type = g.getIntegerProperty("type")
    ub_list = list(set(ub_list))
    ub = g.getBooleanProperty("ub")
    id_ = g.getStringProperty("id")
    for n in g.getEdges():
        if id_[n] in ub_list:
            #type[n] = 1
            ub[n] = True


def speciesToNodes(graph, input_model, ub_sps):
    """
    Create graph nodes for each species of a given libsbml model. Different properties such as size, shape,
    etc. are added to the nodes::

    - *graph* a tulip/mimoza graph
    - *input_model* a libsbml model
    - *ub_sps* a list of species to be flagged as ubiquitous (for aliasing)

    """
    id2n = {}
    for s in input_model.getListOfSpecies():
        # get current species id
        _id = s.getId()
        # get current species name
        name = s.getName()
        # get compartment of current species
        comp = input_model.getCompartment(s.getCompartment())
        # check if species is in ub list: return either true or false
        ub = (_id in ub_sps)

        # add new node to graph and return it (saved in n)
        n = graph.addNode()
        # add compartment id to new node of graph
        graph[graph_properties.COMPARTMENT_ID][n] = comp.getId()
        # add current species id to new node of graph
        graph[graph_properties.ID][n] = _id
        # if current species name exists:
        if name:
            # remove compartment from the name,
            # e.g. H2O [peroxisome] --> H2O
            if comp.getName():
                name = name.replace('[%s]' % comp.getName(), '').strip()
                name = name.replace('[%s]' % comp.getId(), '').strip()
            # add current species name to new node of graph
            graph[graph_properties.NAME][n] = name
        # add type of new node corresponding to the value that is defined for species nodes (1)
        # (reaction nodes have type = 2)
        graph[graph_properties.TYPE][n] = graph_properties.TYPE_SPECIES
        # add ub status (true or false) to new node of graph
        graph[graph_properties.UBIQUITOUS][n] = ub
        # add shape to new node of graph (set to circle shape, can be changed in graph_properties.py)
        graph[graph_properties.VIEW_SHAPE][n] = graph_properties.SPECIES_SHAPE
        # add size to new node of graph (calculated by function get_n_size in file resize.py)
        graph[graph_properties.VIEW_SIZE][n] = resize.get_n_size(graph, n)

        # assign the current species id (key) to the new tulip node n (value) in dict id2n
        id2n[_id] = n
    return id2n


def aliasNode(g, n_alias, target_string, x , y, n_new=None):
    """
    Alias a specific node (reaction nodes and/or metabolite nodes)::

    - *g* a tulip/mimoza graph
    - *n_alias* a tulip node to alias (format: tlp.node(index))
    - *target_string* a name string of the target node connected to n_alias
    - *x* x coordinate (int) of new node position
    - *y* y coordinate (int) of new node position
    - *n_new* newly created tulip node (cloned from n_alias)

    """
    viewName = g.getStringProperty("name")
    edge_direct = g.getStringProperty("edge_directionality")

    # select target node from target_string
    target_node = []
    for n in g.getNodes():
        if viewName[n] == target_string:
            target_node.append(n)

    """
    print("\n \n alias node: {0} {1}".format(n_alias, viewName[n_alias]))
    print("target node: {0} {1} \n".format(target_node[0], viewName[target_node[0]]))
    """

    # iterate over all incoming and outgoing edges of alias node
    for e in g.getInOutEdges(n_alias):
        # take the edge that begins at the target node
        if viewName[g.source(e)] == target_string:
            # if a node with the target name is the source of the edge, it's a product edge (outgoing, away from target)
            edge_direct[e] = "out" # set edge direct property to "out" = away from reaction node
            #print("{0} edge: {1} {2}".format(edge_direct[e], e, viewName[e]))
        if viewName[g.target(e)] == target_string:
            # if a node with the target name is the target of the edge, it's a educt edge (incoming, towards the target)
            edge_direct[e] = "in" # set edge direct property to "in" = towards reaction node
            #print("{0} edge: {1} {2}".format(edge_direct[e], e, viewName[e]))

    # add new node
    n_new = g.addNode()

    # copy properties from old node to new node
    for prop_string in g.getProperties():
        # convert property name (string) to property object
        prop_obj = g.getProperty(prop_string)
        # initially set all properties of new node to prop of old node
        prop_obj[n_new] = prop_obj[n_alias]
        # give new node some different coords to old node
        if prop_string == "viewLayout":
            prop_obj = g.getProperty(prop_string)
            prop_obj[n_new] = tlp.Coord(x, y, 0)
        # set clone id (node to which new node is connected)
        clone_id = g.getStringProperty("clone_id")
        clone_id[n_new] = target_string

    # add new edges and delete old ones
    for e in g.getInOutEdges(n_alias):
        # incoming edge (towards target)
        if edge_direct[e] == "in":
            e_new = g.addEdge(n_new, target_node[0])  # create new edge (from alias node to target)
            edge_direct[e_new] = edge_direct[e]
            viewName[e_new] = viewName[e]  # copy edge name to new edge
            g.delEdge(e)  # delete the old edge
        # outgoing edge (away from target)
        if edge_direct[e] == "out":
            e_new = g.addEdge(target_node[0], n_new)  # create new edge (from target to alias node)
            edge_direct[e_new] = edge_direct[e]
            viewName[e_new] = viewName[e]  # copy edge name to new edge
            g.delEdge(e)  # delete the edge

    return n_new


def delNodes(g, node_list):
    """
    Delete a set of nodes from a graph (by their id)::

    - *g* a tulip/mimoza graph
    - *node_list* a set of nodes

    """
    id_ = g.getStringProperty("id")
    for n in g.getNodes():
        if id_[n] in node_list:
            g.delNode(n)


def addLabels(g, node_type, label_type):
    """
    Add labels to nodes and/or edges of graph::

    - *g* a tulip/mimoza graph
    - *node_type* a predefined node type
    - *label_type* the label type, can be either "name" or "id"

    """
    viewLabel = g.getStringProperty("viewLabel")
    type = g.getIntegerProperty("type")
    view_label_type = g.getStringProperty(label_type)
    for n in g.getNodes():
        if type[n] == node_type:
            viewLabel[n] = view_label_type[n]


def setColours(g, item_list, colour):
    """
    Apply colour to a list of nodes or edges::

     - *g* a tulip/mimoza graph
     - *item_list* a tulip iterator (graph.getEdges or graph.getNodes)
     - *colour* a predefined colour

    """
    # apply colour to graph
    viewColor = g.getColorProperty("viewColor")
    for e in item_list:
        viewColor[e] = colour


def createSubgraph(graph, sg_name, n_list):
    """
    Creates a subgraph from a list of nodes::

    - *graph* tulip graph
    - *sg_name* subgraph name (string)
    - *n_list* a list of nodes

    """
    # define tulip graph properties
    id_ = graph.getStringProperty("id")
    selection = graph.getBooleanProperty("viewSelection")
    # create empty subgraph
    sg = graph.addSubGraph(sg_name)
    # reset selection
    selection.setAllNodeValue(False)
    selection.setAllEdgeValue(False)
    # populate subgraph
    for n in graph.getNodes():
        if id_[n] in n_list or n in n_list:
            selection[n] = True
            sg.addNode(n)
    for e in graph.getEdges():
        if selection[graph.source(e)] and selection[graph.target(e)]:
            selection[e] = True
            if e in graph.getInEdges(graph.source(e)):
                sg.addEdge(graph.target(e), graph.source(e))
            if e in graph.getOutEdges(graph.source(e)):
                sg.addEdge(graph.source(e), graph.target(e))
    # reset selection
    selection.setAllNodeValue(False)
    selection.setAllEdgeValue(False)

    return sg


def layoutSubgraph(sg, sg_layout = None, x_center = None, y_center = None, z_center = 0):
    """
    Apply layout to subgraph and center it around a given point::

    - *sg* a tulip/mimoza subgraph
    - *sg_layout* subgraph layout (string), ex: "Circular (OGDF)" (optional)
    - *x_center* x coordinate of layout center which defines subgraph position (optional)
    - *y_center* y coordinate of layout center which defines subgraph position (optional)
    - *z_center* z coordinate of layout center which defines subgraph position, default 0

    """
    # pseudo-fixed layout method
    # list of available algorithms: tlp.getLayoutAlgorithmPluginsList()
    # apply algorithm to subgraph
    # add default parameters and params as argument (TODO)
    if sg_layout is not None:
        sg.applyLayoutAlgorithm(sg_layout)
    # center layout: center(viewLayout, newCenter, subgraph)
    # newCenter can be given as coordinates (tlp.Coord(x,y,z)) or
    viewLayout = sg.getLayoutProperty("viewLayout")
    # pass a tuple and not three variables individually (TODO)
    if x_center is not None and y_center is not None:
        tlp.LayoutProperty.center(viewLayout, tlp.Coord(x_center, y_center, z_center), sg)

    return sg


