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

import random
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
    # assign all new graph properties to variables
    anc_term = graph.getStringProperty('anc_term')
    anc_id = graph.getStringProperty('anc_id')
    anc_name = graph.getStringProperty("anc_name")
    term = graph.getStringProperty('term')
    c_id = graph.getStringProperty("c_id")
    id_ = graph.getStringProperty('id')
    name = graph.getStringProperty('name')
    type = graph.getIntegerProperty('type')
    rev = graph.getBooleanProperty("rev")
    stoich = graph.getDoubleProperty('stoich')
    transport = graph.getBooleanProperty("tr")
    ub = graph.getBooleanProperty('ub')
    viewLabel = graph.getStringProperty("viewLabel")
    viewLayout = graph.getLayoutProperty("viewLayout")
    viewShape = graph.getIntegerProperty("viewShape")
    viewSize = graph.getSizeProperty("viewSize")
    viewColor = graph.getColorProperty("viewColor")
    fake = graph.getBooleanProperty("fake")
    clone_id = graph.getStringProperty("clone_id")
    formula = graph.getStringProperty("f")
    rel_c_ids = graph.getStringVectorProperty("rel_c_ids")
    edge_direct = graph.getStringProperty("edge_directionality")
    flux = graph.getDoubleProperty("flux")
    selection = graph.getBooleanProperty("viewSelection")

    return anc_term, anc_id, anc_name, term, c_id, id_, name, type, rev, stoich, transport, ub, viewLabel, \
           viewLayout, viewShape, viewSize, viewColor, fake, clone_id, formula, rel_c_ids, edge_direct, flux, selection


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


def aliasNode(g, n_alias, target_string, x=None , y=None, n_new=None):
    """
    Alias a specific node (reaction nodes and/or metabolite nodes)::

    - *g* a tulip/mimoza graph
    - *n_alias* a tulip node to alias (format: tlp.node(index))
    - *target_string* a name string of the target node connected to n_alias
    - *x* x coordinate (int) of new node position
    - *y* y coordinate (int) of new node position
    - *n_new* newly created tulip node (cloned from n_alias)

    """
    id_ = g.getStringProperty("id")
    edge_direct = g.getStringProperty("edge_directionality")
    # select target node from target_string
    target_node = []
    for n in g.getNodes():
        if id_[n] == target_string:
            target_node.append(n)
    # iterate over all incoming and outgoing edges of alias node
    for e in g.getInOutEdges(n_alias):
        # take the edge that begins at the target node
        if id_[g.source(e)] == target_string:
            # if a node with the target name is the source of the edge, it's a product edge (outgoing, away from target)
            edge_direct[e] = "out" # set edge direct property to "out" = away from reaction node
            #print("{0} edge: {1} {2}".format(edge_direct[e], e, viewName[e]))
        if id_[g.target(e)] == target_string:
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
        if prop_string == "viewLayout" and x is not None and y is not None:
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
            id_[e_new] = id_[e]  # copy edge name to new edge
            g.delEdge(e)  # delete the old edge
        # outgoing edge (away from target)
        if edge_direct[e] == "out":
            e_new = g.addEdge(target_node[0], n_new)  # create new edge (from target to alias node)
            edge_direct[e_new] = edge_direct[e]
            id_[e_new] = id_[e]  # copy edge name to new edge
            g.delEdge(e)  # delete the edge
    # delete old node
    # todo: add option (true/false) to either delete old node or keep it
    #g.delNode(n_alias)

    return n_new


def aliasObjFunc(g, obj_re_id, sg_layout, prop_dict):
    """
    Alias the objective reaction (with all reagents)::

    - *g* a tulip/mimoza graph
    - *obj_re_id* id (string) of objective reaction
    - *sg_layout* name of layout algorithm (string) for obj reaction subgraph
    - *prop_dict* dict with graph properties

    """
    # unwrap passed graph props
    id_ = prop_dict["id"]
    viewLayout = prop_dict["viewLayout"]
    viewSize = prop_dict["viewSize"]
    # select obj function nodes
    obj_node_list = []
    for n in g.getNodes():
        # select obj reaction node
        if id_[n] == obj_re_id:
            # add obj reaction node to list
            obj_node_list.append(n)
            # access obj reaction node coords (to layout connected nodes around it)
            x, y, z = viewLayout[n]
            # iterate over all nodes which are connected to obj reaction node
            for node in g.getInOutNodes(n):
                # alias node and move it to a random position close to obj reaction node (quick test hack)
                node_aliased = aliasNode(g, node, obj_re_id, x + random.random() * 50, y + random.random() * 50,
                                             None)
                # save newly created node in a list after aliasing
                obj_node_list.append(node_aliased)
    # move reaction node and all connected, aliased nodes 200 units to the right (x+200) and 200 units up (y+200)
    re_node_coords = []
    for node in obj_node_list:
        x, y, z = viewLayout[node]
        viewLayout[node] = tlp.Coord(x + 50, y + 15, z)  # shift coords here to move the biomass reaction independently
        if id_[node] == obj_re_id:
            # save coords of obj reaction node
            x, y, z = viewLayout[node]
            re_node_coords.append(x)
            re_node_coords.append(y)
            re_node_coords.append(z)
    # turn obj node group into subgraph (but without reaction node, map it seperately!)
    obj_node_list = [n for n in obj_node_list if id_[n] != obj_re_id]
    obj_sg = createSubgraph(g, "obj reaction", obj_node_list, None)
    # increasing circle radius of layout "Circular" only possible by increasing size of nodes
    for n in obj_sg.getNodes():
        viewSize[n] = tlp.Size(3, 3, 0)
    # apply layout to obj subgraph
    layoutSubgraph(obj_sg, sg_layout, None, re_node_coords[0], re_node_coords[1], 0)
    # set colour of new egdes to grey
    setColours(obj_sg, obj_sg.getEdges(), tlp.Color(220, 220, 220))


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


def createSubgraph(graph, sg_name, n_list, user_params=None):
    """
    Creates a subgraph from a list of nodes::

    - *graph* tulip graph
    - *sg_name* subgraph name (string)
    - *n_list* a list of nodes (ids or objects)
    - *user_params* parameter dictionary to overwrite default parameters of layout algorithm (optional)

    """
    # define tulip graph properties
    id_ = graph.getStringProperty("id")
    selection = graph.getBooleanProperty("viewSelection")
    # reset selection
    selection.setAllNodeValue(False)
    selection.setAllEdgeValue(False)
    # select nodes
    for n in graph.getNodes():
        if id_[n] in n_list or n in n_list:
            selection[n] = True
    # select reachable edges from node selection
    # get default parameters of Reachable Sub-Graph algorithm
    params = tlp.getDefaultPluginParameters('Reachable Sub-Graph', graph)
    # tlp.getDefaultPluginParameters("Reachable Sub-Graph")
    #    {'distance': 5L,
    #    'edges direction': 'output edges',
    #    'result': None,
    #    'starting nodes': None}
    # with default parameters
    if user_params is None:
        graph.applyBooleanAlgorithm('Reachable Sub-Graph', params)
    # with user parameters (overwrite default parameters)
    if user_params is not None:
        for key in user_params.keys():
            if key in params:
                params[key] = user_params[key]
                graph.applyBooleanAlgorithm('Reachable Sub-Graph', params)
    # create subgraph from selection
    sg = graph.addSubGraph(selection, sg_name)
    # reset selection
    selection.setAllNodeValue(False)
    selection.setAllEdgeValue(False)

    return sg


def layoutSubgraph(sg, sg_layout=None, user_params=None, x_center=None, y_center=None, z_center=0):
    """
    Apply layout to subgraph and center it around a given point::

    - *sg* a tulip/mimoza subgraph
    - *sg_layout* subgraph layout (string), ex: "Circular (OGDF)" (optional)
    - *user_params* parameter dictionary to overwrite default parameters of layout algorithm (optional)
    - *x_center* x coordinate of layout center which defines subgraph position (optional)
    - *y_center* y coordinate of layout center which defines subgraph position (optional)
    - *z_center* z coordinate of layout center which defines subgraph position, default 0

    """
    # pseudo-fixed layout method
    # list of available algorithms: tlp.getLayoutAlgorithmPluginsList()
    # get default parameters of layout algorithm
    params = tlp.getDefaultPluginParameters(sg_layout, sg)
    # apply algorithm to subgraph
    # with default parameters
    if sg_layout is not None and user_params is None:
        sg.applyLayoutAlgorithm(sg_layout, params)
    # with user parameters (overwrite default parameters)
    if sg_layout is not None and user_params is not None:
        for key in user_params.keys():
            if key in params:
                params[key] = user_params[key]
                sg.applyLayoutAlgorithm(sg_layout, params)
    # center layout: center(viewLayout, newCenter, subgraph)
    # newCenter can be given as coordinates (tlp.Coord(x,y,z)) or
    viewLayout = sg.getLayoutProperty("viewLayout")
    # todo: pass a tuple and not three variables individually
    if x_center is not None and y_center is not None:
        tlp.LayoutProperty.center(viewLayout, tlp.Coord(x_center, y_center, z_center), sg)

    return sg


def colourFBA(g, flux_dict, color_pos, color_neg, color_zero):
    """
    Colour reaction nodes and adjacent edges dependent on flux value::

    - *g* a tulip/mimoza graph
    - *flux_dict* a dictionary with reaction id keys and flux values
    - *color_pos* a color value (RGB) for maximum flux
    - *color_neg* a color value (RGB) for minimum flux
    - *color_zero* a color value (RGB) for zero flux

    """
    # save flux results as graph integer property "flux"
    flux = g.getDoubleProperty("flux")
    id_ = g.getStringProperty("id")
    node_ids = [id_[r] for r in g.getNodes()] # create a list of reaction node ids for entire network (keys)
    node_objs = [r for r in g.getNodes()] # create a list of reaction node objects for entire network (values)
    network_nodes = {}
    for i in range(len(node_ids)):
        network_nodes[node_ids[i]] = node_objs[i]
    for r_id, f in flux_dict.items():
        if r_id in network_nodes.keys():
            flux[network_nodes[r_id]] = f
    # TODO: add flux labels to reaction nodes
    #mvl.addLabels(graph, 2, "flux")
    # -> addLabel only accepts string label properties so far!
    # extract reaction nodes from dict
    pos_flux = []
    neg_flux = []
    zero_flux = []
    for key in flux_dict:
        if flux_dict[key] > 0:
            pos_flux.append(key)
        elif flux_dict[key] < 0:
            neg_flux.append(key)
        else:
            zero_flux.append(key)
    # color total graph (zero fluxes will remain grey)
    setColours(g, g.getNodes(), color_zero)
    setColours(g, g.getEdges(), color_zero)
    # tulip provides a colour scale function!
    # 1. initialize empty color scales
    colorScale_pos = tlp.ColorScale([])
    colorScale_neg = tlp.ColorScale([])
    # 2. fill colorScales (gradiant from white to chosen pos/neg color)
    colorScale_pos.setColorAtPos(0.0, tlp.Color(255, 255, 255))
    colorScale_pos.setColorAtPos(1.0, color_pos)
    colorScale_neg.setColorAtPos(0.0, tlp.Color(255, 255, 255))
    colorScale_neg.setColorAtPos(1.0, color_neg)
    # 3. get colors for flux values
    # color in and out edges of flux nodes and reaction nodes
    # TODO: use numpy to generate ranges for flux ranges numpy.linspace(start, end, number of intervalls)
    viewColor = g.getColorProperty("viewColor")
    for n in network_nodes.values():
        if id_[n] in pos_flux:
            if flux[n] > 0 and flux[n] <= 0.5:
                viewColor[n] = colorScale_pos.getColorAtPos(0.25)
                for e in g.getInOutEdges(n):
                    viewColor[e] = colorScale_pos.getColorAtPos(0.25)
            if flux[n] > 0.5 and flux[n] <= 1:
                viewColor[n] = colorScale_pos.getColorAtPos(0.5)
                for e in g.getInOutEdges(n):
                    viewColor[e] = colorScale_pos.getColorAtPos(0.5)
            if flux[n] > 1 and flux[n] <= 3:
                viewColor[n] = colorScale_pos.getColorAtPos(0.75)
                for e in g.getInOutEdges(n):
                    viewColor[e] = colorScale_pos.getColorAtPos(0.75)
            if flux[n] > 3 :
                viewColor[n] = colorScale_pos.getColorAtPos(1)
                for e in g.getInOutEdges(n):
                    viewColor[e] = colorScale_pos.getColorAtPos(1)
        if id_[n] in neg_flux:
            if flux[n] < 0 and flux[n] <= -0.5:
                viewColor[n] = colorScale_neg.getColorAtPos(0.25)
                for e in g.getInOutEdges(n):
                    viewColor[e] = colorScale_neg.getColorAtPos(0.25)
            if flux[n] < -0.5 and flux[n] <= -1:
                viewColor[n] = colorScale_neg.getColorAtPos(0.5)
                for e in g.getInOutEdges(n):
                    viewColor[e] = colorScale_neg.getColorAtPos(0.5)
            if flux[n] < -1 and flux[n] <= -3:
                viewColor[n] = colorScale_neg.getColorAtPos(0.75)
                for e in g.getInOutEdges(n):
                    viewColor[e] = colorScale_neg.getColorAtPos(0.75)
            if flux[n] < -3 :
                viewColor[n] = colorScale_neg.getColorAtPos(1)
                for e in g.getInOutEdges(n):
                    viewColor[e] = colorScale_neg.getColorAtPos(1)
