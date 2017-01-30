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

import os
import random
from tulip import tlp
from tulipgui import *
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


def createPathwayLists(cmod, path_name, del_ub_sps=True, small_threshold=30, ub_sps=None):
    """
    Creates a dictionary containing group id, metabolite ids and reaction ids of a given pathway (if it exists as
    a group in the SBML file) and returns dict content as lists path_sps (metabolite ids) and path_re (reaction ids)::

    - *cmod* a cbmpy model object
    - *path_name* a string of the pathway/group name
    - *del_ub_sps* [default=True] delete ids from highly connected species from pathway lists
    - *small_threshold* an (optional) integer threshold level for filtering highly connected metabolites (default
    level: 30)
    - *ub_sps* a (optional) list of highly connected species ids to remove

    """
    # create a dictionary to store model groups (key: group name - value: group id)
    gids = {}
    for g in cmod.groups:
        if g.getName() == path_name:
            gids[path_name] = {}
            gids[path_name]['id'] = g.getId()
            gids[path_name]['reactions'] = g.getMemberIDs()
            # leave it blank for now - usually groups only contain reaction ids (not metabolites)
            gids[path_name]['metabolites'] = []
    # get all member species per group for lists of reaction nodes
    path_re = gids[path_name]['reactions']
    # for each reaction get the associated species for lists of metabolite nodes
    path_sps = []
    for r in path_re:
        path_sps.extend(cmod.getReaction(r).getSpeciesIds())
    path_sps = list(set(path_sps))
    # add species lists to dict
    gids[path_name]['metabolites'] = path_sps
    # get ubiquitous small metabolites by filtering for high connectivity (default: > 30)
    # get stoich matrix (rows = species, columns = reactions)
    if del_ub_sps:
        n = cmod.N.array
        # convert all matrix entries which are not 0 to 1 (so that one involvement in a reaction is represented by a 1)
        n[n!=0] = 1 # convert stoich matrix to adjacency matrix
        # save number of connections (= involvment in a reaction) per species (= row in stoich matrix) in a dict
        # dict: key = species id of a row; value = integer value of reactions that species is involved in (connections )
        n_r_conn = {}
        for row in range(n.shape[0]):
            n_r_conn[cmod.N.row[row]] = round(n[row, :].sum(), 0)
        # if the number of connections for a species is bigger than the given threshold
        # -> put the repsective species id in the list small_sps
        if ub_sps:
            small_sps = ub_sps
        # if no user input is provided, try to automatically determine ubiquitous species based on connectivity
        else:
            # take user threshold if provided
            if small_threshold:
                small_sps = [sps_id for sps_id, num_re in n_r_conn.items() if num_re > small_threshold]
            # if not user threshold is provided make automatic threshold level (dependent on model size)
            else:
                small_threshold = len(n_r_conn.values()) * 0.2
                small_sps = [sps_id for sps_id, num_re in n_r_conn.items() if num_re > small_threshold]
        # remove the small metabolites from the pathways
        for ubs in small_sps:
            if ubs in path_sps:
                path_sps.remove(ubs)

    return path_sps, path_re


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
    t = graph.getIntegerProperty('type')
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
    path = graph.getStringProperty("pathway")
    font_size = graph.getIntegerProperty("viewFontSize")

    return anc_term, anc_id, anc_name, term, c_id, id_, name, t, rev, stoich, transport, ub, viewLabel, \
           viewLayout, viewShape, viewSize, viewColor, fake, clone_id, formula, rel_c_ids, edge_direct, flux, \
           selection, path, font_size


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


def speciesToNodes(graph, input_model, ub_sps, graph_props):
    """
    Create graph nodes for each species of a given libsbml model. Different properties such as size, shape,
    etc. are added to the nodes::

    - *graph* a tulip/mimoza graph
    - *input_model* a libsbml model
    - *ub_sps* a list of species to be flagged as ubiquitous (for aliasing)
    - *graph_props* a dict of graph properties

    """
    id2n = {}
    for s in input_model.species:
        # get current species id
        _id = s.getId()
        # get current species name
        name = s.getName()
        # get compartment id of current species
        comp_id = s.getCompartmentId()
        # check if species is in ub list: return either true or false
        ub = (_id in ub_sps)
        # add new node to graph and return it (saved in n)
        n = graph.addNode()
        # add compartment id to new node of graph
        if comp_id:
            graph_props["c_id"][n] = comp_id
        # add current species id to new node of graph
        graph_props["id"][n] = _id
        # add current species name to new node
        graph_props["name"][n] = name
        # add type of new node
        #  (species nodes = 1, reaction nodes = 2)
        graph_props["type"][n] = graph_properties.TYPE_SPECIES
        # add ub status (true or false) to new node of graph
        graph_props["ub"][n] = ub
        # add shape to new node of graph
        # (species nodes = circle, reaction nodes = square)
        graph_props["viewShape"][n] = graph_properties.SPECIES_SHAPE
        # add size to new node of graph (calculated by function get_n_size in file resize.py)
        graph_props["viewSize"][n] = resize.get_n_size(graph, n)
        # assign the current species id (key) to the new tulip node n (value) in dict id2n
        id2n[_id] = n

    return id2n


def reactionsToNodes(get_r_comp, graph, id2n, input_model):
    get_sp_comp = lambda _id: input_model.getSpecies(_id).getCompartmentId()

    def link_reaction_to_species(reaction_node, sp_ref, all_comps, is_substrate=None, reversible=None):
        s_id = sp_ref.getId()
        all_comps.add(get_sp_comp(s_id))
        species_node = id2n[s_id]
        # create edge
        e = graph.addEdge(species_node, reaction_node) if is_substrate else graph.addEdge(reaction_node, species_node)
        graph[graph_properties.NAME][e] = input_model.getSpecies(s_id).getName()
        ub = graph[graph_properties.UBIQUITOUS][species_node]
        graph[graph_properties.UBIQUITOUS][e] = ub
        graph[graph_properties.REVERSIBLE][e] = reversible and not ub

    for r in input_model.reactions:
        name = r.getName()
        # create reaction node and add basic graph properties
        n = graph.addNode()
        graph[graph_properties.ID][n] = r.getId()
        graph[graph_properties.NAME][n] = name
        graph[graph_properties.TYPE][n] = graph_properties.TYPE_REACTION
        # get reversible and irreversible reactions
        rev_re = [re for re in input_model.getReversibleReactionIds()]
        irrev_re = [reac for reac in input_model.getIrreversibleReactionIds()]
        # add reversibility graph property
        if r.getId() in irrev_re:
            graph[graph_properties.REVERSIBLE][n] = False
        elif r.getId() in rev_re:
            graph[graph_properties.REVERSIBLE][n] = True
        # add shape and size graph properties
        graph[graph_properties.VIEW_SHAPE][n] = graph_properties.REACTION_SHAPE
        graph[graph_properties.VIEW_SIZE][n] = resize.get_n_size(graph, n)
        # add edges with function link_reaction_to_species
        # -> preserve type of species and directionality
        all_comps = set()
        for sp_ref in r.getSpeciesObj():
            # species is substrate of irreversible reaction
            if sp_ref.getId() in r.getSubstrateIds() and r.getId() in irrev_re:
                link_reaction_to_species(n, sp_ref, all_comps, is_substrate=True, reversible=False)
            # species is substrate of reversible reaction
            if sp_ref.getId() in r.getSubstrateIds() and r.getId() in rev_re:
                link_reaction_to_species(n, sp_ref, all_comps, is_substrate=True, reversible=True)
            # species is product of irreversible reaction
            if sp_ref.getId() in r.getProductIds() and r.getId() in irrev_re:
                link_reaction_to_species(n, sp_ref, all_comps, is_substrate=False, reversible=False)
            # species is product of reversible reaction
            if sp_ref.getId() in r.getProductIds() and r.getId() in rev_re:
                link_reaction_to_species(n, sp_ref, all_comps, is_substrate=False, reversible=True)
        # add graph property transport to reaction node (boolean value)
        graph[graph_properties.TRANSPORT][n] = len(all_comps) > 1
        # add compartment id to reaction node
        graph[graph_properties.COMPARTMENT_ID][n] = get_r_comp(all_comps)


# aliasNode() function description:
# ---------------------------------
# the function aliasNode does only a very specific thing!
# -> it cuts exactly one edge between a node and a target and then adds a new egde with a new node to the target
# this new node therefore has only 1 edge attached to it, namely the new edge to the target
# => this is only helpful for the obj reaction because only there does the new node not need any other edges
#
# -> you can't use the node created by aliasNode for a proper pathway like Glycolysis because all the other edges to
# other reaction and metabolite nodes are not created for the new node, only 1 to the specified target
def aliasNode(g, n_alias, n_target, x=None, y=None, n_new=None):
    """
    Alias a specific node (reaction nodes and/or metabolite nodes)::

    - *g* a tulip/mimoza graph
    - *n_alias* a tulip node to alias (format: tlp.node(index))
    - *n_target* a tulip node connected to n_alias (the target)
    - *x* x coordinate (int) of new node position
    - *y* y coordinate (int) of new node position
    - *n_new* newly created tulip node (cloned from n_alias)

    """
    id_ = g.getStringProperty("id")
    edge_direct = g.getStringProperty("edge_directionality")
    # select target node from target_string
    target_node = []
    for n in g.getNodes():

        # TODO: use node object instead of id for target selection (because id is not unique sometimes!)

        if n == n_target:
            target_node.append(n)
    # iterate over all incoming and outgoing edges of alias node
    for e in g.getInOutEdges(n_alias):
        # take the edge that begins at the target node
        if g.source(e) == n_target:
            # if a node with the target name is the source of the edge, it's a product edge (outgoing, away from target)
            edge_direct[e] = "out" # set edge direct property to "out" = away from reaction node
            #print("{0} edge: {1} {2}".format(edge_direct[e], e, viewName[e]))
        if g.target(e) == n_target:
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
        id_ = g.getStringProperty("id")
        clone_id[n_new] = id_[n_target]
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


def layoutObjFunc(g, obj_re_id, sg_layout, prop_dict, obj_name="objective reaction"):
    """
    Layout the objective reaction (with all reagents)::

    - *g* a tulip/mimoza graph
    - *obj_re_id* id (string) of objective reaction
    - *sg_layout* name of layout algorithm (string) for obj reaction subgraph
    - *prop_dict* dict with graph properties
    - *obj_name*

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
            obj_re_node = n
            # add obj reaction node to list
            obj_node_list.append(n)
            # access obj reaction node coords (to layout connected nodes around it)
            x, y, z = viewLayout[n]
            # iterate over all nodes which are connected to obj reaction node
            for node in g.getInOutNodes(n):
                # alias node and move it to a random position close to obj reaction node (quick test hack)
                node_aliased = aliasNode(g, node, obj_re_node, x + random.random() * 50, y + random.random() * 50,
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
    # apply layout to obj subgraph
    layoutSubgraph(obj_sg, sg_layout, None, re_node_coords[0], re_node_coords[1], 0)
    # set colour of new egdes to grey
    setColours(obj_sg, obj_sg.getEdges(), tlp.Color(220, 220, 220))
    # set subgraph name
    obj_sg.setName(obj_name)

    return obj_sg


def layoutSeparateReactions(g, re_list, re_node_obj, sg_name, layout_alg, graph_props):
    """
    Separate reactions, then create a subgraph with the isolated reactions next to each other and apply a
    layout algorithm on the group of non-connected reactions::

    - *g* a tulip/mimoza graph
    - *re_list* a list of reaction ids (strings)
    - *re_node_obj* a list of tulip node objects
    - *sg_name* the name of the new subgraph
    - *layout_alg* a tulip layout algorithm for the new subgraph
    - *graph_props* a dict of graph properties

    """
    # alias a group of reactions
    # I.) construct first dictionary with key: reaction id and value = reaction node object
    # re_node_obj contains list of ex re tulip nodes (value list)
    # (without the ones excluded for having been already layouted)
    # re_list contains the ids (key list)
    # merge the key and value lists to a dictionary
    re_id_dict = {}
    for i in range(len(re_list)):
        re_id_dict[re_list[i]] = re_node_obj[i]

    # II.) construct second dictionary with key = ex reaction id and value = list of connected out nodes
    # re_list already contains the ids (key list)
    # create a dict with reaction ids as keys and a list of connected in and out nodes as value
    re_con_dict = dict((graph_props["id"][r], [n for n in g.getInOutNodes(r) if graph_props["id"][r] in
                                               re_list]) for r in re_node_obj)

    # III.) call aliasNode function for every node connected to an reaction node
    re_alias_result_nodes = []
    for r_id in re_con_dict:
        for node in re_con_dict[r_id]:
            n_new = aliasNode(g, node, re_id_dict[r_id], None, None)
            re_alias_result_nodes.append(n_new)
    # IV.) add newly created, aliased nodes to a subgraph and layout it
    re_sg = createSubgraph(g, sg_name, re_alias_result_nodes, None)
    # createSubgraph with only list of connected metabolite nodes (re_alias_result_nodes) does the same as:
    # re_sg = g.inducedSubGraph(re_alias_result_nodes+re_node_obj)
    #  => inducedSubGraph with reaction AND metabolite nodes
    setColours(re_sg, re_sg.getEdges(), tlp.Color(220, 220, 220))
    layoutSubgraph(re_sg, layout_alg, None, 0, 0, 0)

    return re_sg


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


# TODO: write a function that takes in model ids and returns a list of the associated nodes and/or egdes
# TODO: use graph_elements dict

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
    Creates a reachable subgraph from a list of nodes (may include more elements than are stored in n_list!).
    For creation of a subgraph strictly from a set of nodes, use graph.inducedSubGraph(node_list)::

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


def layoutPathway(g, graph_elements, graph_props, id_list, sg_name, sg_layout):
    """
    Create and layout a subgraph based initially on a list of model ids::

    - *g* the initial tulip graph generated by Mimoza
    - *graph_elements* a dict of all initial graph elements (nodes and their layout status)
    - *graph_props* a dict of all graph properties
    - *id_list* a list of model ids (strings)
    - *sg_name* a name for the subgraph (string)
    - *sg_layout* a tulip layout algorithm to be applied on the subgraph (string)

    """
    # generate (initial) list of nodes for pathway elements by going through lists containing selected ids
    # this list will also include nodes that already have been used for layouting
    # only take nodes that were initially generated by Mimoza (not new ones that are created by calling this function!)
    nodes = [n for n in g.getNodes() if graph_props["id"][n] in id_list and n == graph_elements[graph_props["id"][n]][
        "initial_node"]]
    # check for nodes which have already been used for layouting (with the flag in dict graph_elements[node][
    # 'layout'] - boolean value)
    # save the nodes which are excluded because they have already been used for layouting in a list for later use
    ex_nodes = [n for n in nodes if graph_elements[graph_props["id"][n]]['layout'] is True]
    # save the connectivity of the excluded edges
    edge_dict = {}
    for n in ex_nodes:
        edge_dict[graph_props["id"][n]] = {}
        edge_dict[graph_props["id"][n]]["in"] = []
        edge_dict[graph_props["id"][n]]["out"] = []
        for e in g.getInOutEdges(n):
            if e in list(g.getInEdges(n)) and g.source(e) in nodes and g.target(e) in nodes:
                edge_dict[graph_props["id"][n]]["in"].append([g.source(e),g.target(e)])
            if e in list(g.getOutEdges(n)) and g.source(e) in nodes and g.target(e) in nodes:
                edge_dict[graph_props["id"][n]]["out"].append([g.source(e), g.target(e)])
    # exclude the already layouted nodes (ex_nodes) from the main node list (nodes)
    for n in list(nodes):
        if n in ex_nodes:
            nodes.remove(n)
    # create subgraph with inducedSubgraph function (automatically add all existing edges between seletced nodes aswell)
    # this graph will not have any nodes that have been used for layouting before
    sg = g.inducedSubGraph(nodes)
    # name subgraph
    sg.setName(sg_name)
    # if any nodes have been excluded due to the fact that they already have been used for layouting
    # -> recreate these nodes
    # => they also need to have all the edges connected to them! (as it was the case before)
    # 1.) add new nodes (as many as were excluded)
    new_nodes = sg.addNodes(len(ex_nodes))
    # give new nodes the same properties as old ones
    for i in range(len(new_nodes)):
        for prop_string in g.getProperties():
        # convert property name (string) to property object
            prop_obj = g.getProperty(prop_string)
            # set all properties of new node to prop of old node
            prop_obj[new_nodes[i]] = prop_obj[ex_nodes[i]]
    # 2.) re-establish connectivity between the nodes (is stored in model file) and create the necessary edges
    for n in new_nodes:
        for e_type, conn_list in edge_dict[graph_props["id"][n]].items():
            # replace old nodes in edge dict (which were excluded) with newly created nodes that take their place
            # conn_list = a list of anchor node lists (anchor_pair_list)
            # anchor_pair_list = a list of 2 elements: source and target node of edge
            for anchor_pair_list in list(conn_list):
                if len(anchor_pair_list) != 0:
                    for node in list(anchor_pair_list):
                        if node in ex_nodes:
                            node_index = anchor_pair_list.index(node)
                            anchor_pair_list.insert(node_index, n)
                            anchor_pair_list.remove(node)
            try:
                if e_type == "out":
                    for anchor_pair_list in conn_list:
                        sg.addEdge(anchor_pair_list[0], anchor_pair_list[1])
                if e_type == "in":
                    for anchor_pair_list in conn_list:
                        sg.addEdge(anchor_pair_list[0], anchor_pair_list[1])
            except Exception as ex:
                print(ex)
                print("layoutPathway", sg_name)
                print(n)
                print(e_type)
                print(conn_list)
    # layout subgraph
    sg = layoutSubgraph(sg, sg_layout)
    # flag all the nodes of the glycolysis subgraph as "already used" and add pathway info ( = subgraph name) to their
    # graph properties
    for n in sg.getNodes():
        graph_elements[graph_props["id"][n]]['layout'] = True
        graph_props["pathway"][n] = sg_name
    # add pathway info to egdes aswell
    for e in sg.getEdges():
        graph_props["pathway"][e] = sg_name

    return sg


def dumpToGraph(input_graph, input_graph_props, new_graph_name, sg_selection, sg_labels=True, transport_radial=True):
    """
    Take subgraphs from one graph object and copy them to another graph object::

    - *input_graph* a tulip graph object to take subgraphs from
    - *input_graph_props* a dict of input graph properties
    - *new_graph_name* a name of the newly created tulip graph that recieves the subgraphs (string)
    - *sg_selection* a list of ids of subgraphs which are attached to the input_graph and are to be copied
    - *sg_labels* [default=True] a boolean parameter controlling the addition of label nodes to the subgraphs after
    they have been copied
    - *transport_radial* [default=True] a boolean parameter controlling whether or not the Transport subgraph is
    layouted using the 'Tree Radial' Algorithm

    """
    # create new output graph
    new_graph = tlp.newGraph()
    # create props for new combined output graph
    # if a completely new graph object is created, the properties need to be generated again because they are tied
    # to a graph object (properties of a graph can be used for its subgraphs but thats just because subgraphs are only partial
    # views of a whole graph, not new objects)!
    new_graph_tup = create_props(new_graph)
    # save graph props in dictionary
    new_graph_props = {}
    for prop in new_graph_tup:
        new_graph_props[prop.getName()] = prop
    # add name to new graph
    new_graph.setName(new_graph_name)
    # copy subgraphs to output graph 'new_graph' and recreate subgraphs from selection
    graph_sg_list_new_graph = [sg for sg in input_graph.getSubGraphs()
                               if sg.getName() in sg_selection]
    for sg in graph_sg_list_new_graph:
        # reset selection
        input_graph_props["viewSelection"].setAllNodeValue(False)
        input_graph_props["viewSelection"].setAllEdgeValue(False)
        # add pathway to output graph
        tlp.copyToGraph(new_graph, sg, None, input_graph_props["viewSelection"])
        new_graph_sg = new_graph.addSubGraph(input_graph_props["viewSelection"], sg.getName())
        # reset selection
        input_graph_props["viewSelection"].setAllNodeValue(False)
        input_graph_props["viewSelection"].setAllEdgeValue(False)
    # separate overlapping individual pathways
    new_graph.applyLayoutAlgorithm("Connected Component Packing")
    # layout Transport subgraph
    transport_radial = True
    if transport_radial:
        for sg in new_graph.getSubGraphs():
            if sg.getName() == "Transport":
                user_params = {"layer spacing": 60, "node spacing": 80}
                layoutSubgraph(sg, "Tree Radial", user_params, 75, 100, 0)
    # TODO: change type of labels between id or name for all nodes of the entire graph (control with use_names=Fales)
    # add pathway labels (invisible nodes with visible label)
    # TODO: link the pathway label to a nodeon the border of the pathway (so that it doesn't overlap and is visible
    # TODO: in the whole output graph)
    # TODO: to identify border label -> has only In or only Out edges?
    if sg_labels:
        for sg in new_graph.getSubGraphs():
            n = sg.addNode()
            # TODO: labels have type 0 for identfication!
            # TODO: add graph property to be able to flag labels so that we can exclude them from FBA coloring later
            # TODO: (so that the white color for background blending is not overwritten by FBA Coloring)
            new_graph_props["viewLabel"][n] = sg.getName()
            new_graph_props["viewColor"][n] = tlp.Color(255, 255, 255) # set node color to white (invis on white bg)
            new_graph_props["viewSize"][n] = tlp.Size(5, 5, 1) # normal node = (2.5, 2.5, 1); bigger -> (5, 5, 1)
            new_graph_props["viewFontSize"][n] = 16 # default is 18

    # reset all selection
    new_graph_props["viewSelection"].setAllNodeValue(False)
    new_graph_props["viewSelection"].setAllEdgeValue(False)

    return new_graph, new_graph_props


def colourFBA(g, flux_dict, color_pos, color_neg, color_zero):
    """
    Colour reaction nodes and adjacent edges dependent on flux value::

    - *g* a tulip/mimoza graph
    - *flux_dict* a dictionary with reaction id keys and flux values
    - *color_pos* a color value (RGB) for maximum flux
    - *color_neg* a color value (RGB) for minimum flux
    - *color_zero* a color value (RGB) for zero flux

    """
    # get necessary graph properties
    # save flux results as graph integer property "flux"
    flux = g.getDoubleProperty("flux")
    id_ = g.getStringProperty("id")
    viewColor = g.getColorProperty("viewColor")
    node_ids = [id_[r] for r in g.getNodes()] # create a list of reaction node ids for entire network (keys)
    node_objs = [r for r in g.getNodes()] # create a list of reaction node objects for entire network (values)
    network_nodes = {}
    for i in range(len(node_ids)):
        network_nodes[node_ids[i]] = node_objs[i]
    for r_id, f in flux_dict.items():
        if r_id in network_nodes.keys():
            flux[network_nodes[r_id]] = f
    # creates a regular color scale going from neg color to zero color to pos color
    #colorScale = [color_neg, color_zero, color_pos]
    # new FAME rainbow coloring: only apply color to nodes with flux
    # don't fade into zero grey! -> have a clear difference between zero flux and flux != 0
    colorScale = tlpgui.ColorScalesManager.getColorScale('rainbow')
    # a dictionary can also be used (allow to define non regular scales)
    # get min and max flux values
    # colorMap = {0: tlp.Color.Blue, 0.5: tlp.Color.Yellow, 1.0: tlp.Color.Red}
    # set the parameters for the 'Color Mapping' color algorithm
    colorMappingParams = tlp.getDefaultPluginParameters('Color Mapping', g)
    colorMappingParams['input property'] = flux
    colorMappingParams['color scale'] = colorScale
    # overwrite min and max values to make color scale symmetric -> color zero needs to be applied to zero and not a
    # non zero middle value between real min and max -> mid value between -1 and 1 is 0
    colorMappingParams['override minimum value'] = True
    colorMappingParams['minimum value'] = -1
    colorMappingParams['override maximum value'] = True
    colorMappingParams['maximum value'] = 1
    # compute the color mapping
    g.applyColorAlgorithm('Color Mapping', colorMappingParams)
    # color all nodes with flux == 0 grey
    setColours(g, [n for n in g.getNodes() if flux[n] == 0], color_zero)
    # color all edges grey (color of edges connected to reactions with flux will be overwritten by node color)
    setColours(g, g.getEdges(), color_zero)
    # apply reaction node color to in and out egdes if the reaction node has a flux value != 0
    for n in g.getNodes():
        if flux[n] != 0:
            for e in g.getInOutEdges(n):
                viewColor[e] = viewColor[n]


def exportGML(graph, model_id, output_graph_name, gml_path):
    """
    Exports a tulip graph as a GML file::

    - *graph* a tulip graph (input graph to get elements and properties from)
    - *model_id* id of the libSBML model
    - *output_graph_name* a name for the output graph G (has unique node ids), file type: .tlp
    - *gml_path* a path for the output GML file (string), file type: .gml

    """
    # create new output graph 'GML'
    GML = tlp.newGraph()
    GML.setName("GML output graph")
    # if a completely new graph object is created, the properties need to be generated again because they are tied to a
    # graph object (properties of a graph can be used for its subgraphs but thats just because subgraphs are only partial
    # views of a whole graph, not new objects)
    GML_props_tup = create_props(GML)
    GML_props = {}
    for prop in GML_props_tup:
        GML_props[prop.getName()] = prop
    # copy elements from the graph 'graph' to the graph 'GML'
    tlp.copyToGraph(GML, graph)
    # change nodes to be unique (add integers to node ids and labels)
    i = 0
    for n in GML.getNodes():
        GML_props["id"][n] = GML_props["id"][n] + str(i)
        GML_props["viewLabel"][n] = GML_props["viewLabel"][n] + str(i)
        i += 1
    # save new output graph 'GML'
    tlp.saveGraph(GML, '{}'.format(model_id) + '_' + output_graph_name)
    # export new output graph 'GML' as .gml file (graphic markup language)
    params = tlp.getDefaultPluginParameters('GML Export', graph)
    tlp.exportGraph('GML Export', GML, gml_path, params)


def containsSublist(l):
    """
    Check if a list contains sublists. Return 'True' if list contains sublists::

    - *l* a list

    """
    # only add string 'list' to l_item_types if any of the items of l is a list, otherwise the l_item_types os empty
    l_item_types = ['list' for item in l if isinstance(item, list)]
    if len(l_item_types) == 0: #if l_item_types is epmty, no list has been found inside l
        return False
    else: #at least one list has been found inside l
        return True