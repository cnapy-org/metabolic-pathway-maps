# Test transition from sbml file to structured tulip graph, sbml with layout and sbgn file
# necessary  dependencies (Wind64 Anaconda): pip install geojson tarjan libsbgnpy tulip-python
# for SQL database: install sqlite 3 (https://www.tutorialspoint.com/sqlite/sqlite_installation.htm)

import os
# get current directory
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))

# tulip imports
import sys
sys.path.append('C:\\Program Files\\Tulip-4.10.0\\bin')
from tulipgui import * # for SVG export
from tulip import tlp

# modules for SBML model interaction
import libsbml
import cbmpy

# SQL database API
import sqlite3

# necessary Mimoza modules for graph visualization
import sbml2tlp
import tulip_graph2geojson
# sbml and sbgn export
import sbml_manager
import sbgn_helper

# module with own functions
import MVLayout as mvl



############################
# global control variables #
############################

# apply aliasing to every node of the graph with a certain connectivity + option to alias biomass (result in maingraph)
# if auto alias is active the automatically generated layout will still be overwritten if any other layouting
# functions are active aswell
AUTO_ALIAS = False

# layout annotated groups of sbml file as seperate pathways
PATHWAY_LISTS = True
LAYOUT_PATHWAYS = True # needs PATHWAY_LISTS true, so that necessary node lists are generated (for pathways)

# layout the exchange reactions and biomass seperately (results in subgraphs: only ex re, only biomass and both)
LAYOUT_EX_RE_AND_BIOMASS = True

# export generated graph as .gml file
EXPORT_GML = True
# export generated graph as .svg file
EXPORT_SVG = True

# colour the reaction nodes + in and outcoming edges based on flux values obtained from FBA
COLOUR_FBA = False



#############################################
# read SBML with libsbml and cbmpy, do FBA  #
#############################################

# read SBML file (with libsbml) # Auto alias settings
#sbml = os.path.join(cDir, "TestModels", "toy_model_biofuel.xml") # 4, 4, False
#sbml = os.path.join(cDir, "TestModels", "e_coli_core.xml") # 6, 4, True
#sbml = os.path.join(cDir, "TestModels", "Ecoli_Carlson2003.xml") # 6, 4, True
#sbml = os.path.join(cDir, "TestModels", "Hpylori_iIT341.xml") # 6, 4, True
#sbml = os.path.join(cDir, "TestModels", "iJO1366.xml") # 6, 4, True
#sbml = os.path.join(cDir, "TestModels", "Lactococcus_lactis_MG1363.xml") # ?, ?, set it to 'True' probably
sbml = os.path.join(cDir, "TestModels", "iIT341.xml") # BiGG2 Helicobacter pylori 26695

# create libsbml model
reader = libsbml.SBMLReader()
doc = reader.readSBML(sbml)
model = doc.getModel()
model_id = model.getId()

# create CBMPy model
cmod = cbmpy.readSBML3FBC(sbml)
if sbml == os.path.join(cDir, "TestModels", "iJO1366.xml"):
    # this is a manual hack for now, but there should be a way to do it on the display level
    cmod.deleteSpecies('M_h_c', also_delete='reagents')
    cmod.deleteSpecies('M_h2_c', also_delete='reagents')
    cmod.deleteSpecies('M_h2o_c', also_delete='reagents')
    cmod.deleteReactionAndBounds('R_BIOMASS_Ec_iJO1366_WT_53p95M')
    cbmpy.writeSBML3FBCV2(cmod, '_temp.xml')
    sbml = '_temp.xml'

# call stoichiometric matrix
stM = cmod.N
# stM = cmod.buildStoichMatrix(matrix_type = 'numpy', only_return=True)

# get objective reaction id
objective_re_id = str(cmod.getActiveObjectiveReactionIds()[0])

# analyze model (FBA)
cbmpy.doFBA(cmod)
fba_results = {}
for r in cmod.reactions:
    fba_results[r.getId()] = round(r.getValue(), 6)



#################################################
# select nodes for aliasing (different methods) #
#################################################

# 1.) create alias node list based on connectivity <-> user defined thresholds (aliasing by Mimoza)
ub_sp_list = []
ub_re_list = []

if AUTO_ALIAS:
    sp_list, re_list = mvl.identifyConnectedNodes(cmod, 6, 4, True)
    ub_sp_list.extend(sp_list)
    ub_re_list.extend(re_list)

# combine alias species and alias reaction lists (for Mimoza node cloning functions)
ub_list = ub_sp_list + ub_re_list

# add ubiquitous (alias) species group to sbml file and save it as "..._with_groups.xml"
# create group for ubiquitous species
cmod.createGroup('g_ubiquitous_sps')
g1 = cmod.getGroup('g_ubiquitous_sps')
g1.setName('ubiquitous species')
ub_sps = [cmod.getSpecies(i) for i in ub_sp_list]
ub_react = [cmod.getReaction(i) for i in ub_re_list]
# add reactions and species to group
g1.addMember(ub_sps+ub_react)
# add RDF annotation to group:
# MIRIAM annotation and SBO term (important for drawing ub species seperately)
g1.setSBOterm('SBO:0000248')
g1.addMIRIAMannotation('isDescribedBy', 'ChEBI', 'CHEBI:000000') #dummy annotation
g1.miriam.isDescribedBy[0] = "http://identifiers.org/ubiquitous" #change dummy to fit mimoza
cbmpy.writeSBML3FBC(cmod, '{}_with_groups.xml'.format(model_id), add_groups=True)
groups_sbml = os.path.join(cDir, "{}_with_groups.xml".format(model_id))


# 2.) create alias node list based on biological pathways <-> group info in sbml - if available (aliasing by own
# functions)

# -> IF GROUPS DON'T EXIST: provide pathway lists manually (example: e coli carlson model):
if PATHWAY_LISTS and sbml == os.path.join(cDir, "TestModels", "Ecoli_Carlson2003.xml"):
    # create pathway dictionary (key: pathway name - value: tupel of two lists)
    # path_dict = {"glycolysis": ([list of metabolite nodes],[list of reaction nodes]), "TCA Cycle": ([metabolite
    # nodes],[reaction nodes]), ...}
    glyc_sps = []
    tca_sps = []
    glyc_re = []
    tca_re = []
    small_sps = []

    path_dict = {
        "Glycolysis":
            {
            "metabolites": ["GLU_6_P", "FRU_6_P", "FRU_BIS_P", "DHAP", "GA_3P", "PG", "PEP"],
            "reactions": ["R2r", "R3", "R4", "R5r", "R6r", "R7r", "R8r", "R9", "RR9"]
            },
        "TCA Cycle":
            {
            "metabolites": ["CITRATE", "ISOCIT", "AKG", "SUCC_CoA", "SUCC", "FUMARATE", "MALATE", "OXALO"],
            "reactions": ["R23r", "R24", "R25", "R26r", "R27r", "R28r", "R29r", "R22"]
            }
        }
    # ubiquitous small metabolites
    small_sps.extend(["ATP", "ADP", "NAD", "NADH", "CO2", "NH3"])
    # extract pathway lists (for LAYOUT_GLYC_TCA)
    # lists of metabolite nodes
    glyc_sps.extend(path_dict["Glycolysis"]["metabolites"])
    tca_sps.extend(path_dict["TCA Cycle"]["metabolites"])
    # lists of reaction nodes
    glyc_re.extend(path_dict["Glycolysis"]["reactions"])
    tca_re.extend(path_dict["TCA Cycle"]["reactions"])
    # create list of all reaction and species nodes that are not part of any pathway
    not_spec_sps = [sps for sps in cmod.getSpeciesIds() if sps not in glyc_sps or sps not in tca_sps]
    not_spec_re = [re for re in cmod.getReactionIds() if re not in glyc_re or re not in tca_re]
    # create groups for evey pathway dict entry (not necessary if provided by sbml)
    for path in path_dict:
        cmod.createGroup(path)
        group = cmod.getGroup(path)
        group.setName(path)
        group.addMember([cmod.getSpecies(s) for s in path_dict[path]["metabolites"]] + [cmod.getReaction(r) for r in
                                                                                        path_dict[path]["reactions"]])
    # create group for all other species and reactions
    cmod.createGroup('not_specified')
    path0 = cmod.getGroup('not_specified')
    path0.setName('Not specified')
    path0.addMember([cmod.getSpecies(s) for s in not_spec_sps] + [cmod.getReaction(r) for r in not_spec_re])
    cmod.createGroup('glycolysis')


# -> IF GROUPS EXIST: fetch groups from sbml (this doesn't use ub_sps_list or ub_re_list!, specific to iIT341)
# create dictionary of groups where each group name is a key and the values are again dicts with three keys: 'id',
# 'reactions' and 'metabolites'
path_dict = {}
if PATHWAY_LISTS and sbml == os.path.join(cDir, "TestModels", "iIT341.xml"):

    # create id lists for a few selected pathways of model iIT341
    #glyc_sps, glyc_re = mvl.createPathwayLists(cmod, 'Glycolysis')
    #tca_sps, tca_re = mvl.createPathwayLists(cmod, 'TCA cycle')
    #fas_sps, fas_re = mvl.createPathwayLists(cmod, 'Fatty Acid Synthesis')
    #uc_sps, uc_re = mvl.createPathwayLists(cmod, 'Urea Cycle')
    #am_sps, am_re = mvl.createPathwayLists(cmod, 'Aminosugar Metabolism')

    # loop over all groups in iIT341 and automatically create id lists for all of them
    # save all these lists in a new dictionary
    for group in cmod.groups:
        path_sps, path_re = mvl.createPathwayLists(cmod, group.name)
        path_dict.update({group.name: {"metabolites": path_sps, "reactions": path_re}})


print("Pathway lists created.")



##################################################
# create Tulip graph and generate graph objects  #
##################################################

# create new (empty) tulip graph
graph = tlp.newGraph()
graph.setName(model.getId())

# create and store graph properties for new tulip graph elements (nodes and egdes)
props_tup = mvl.create_props(graph)
# save graph props in dictionary
graph_props = {}
for prop in props_tup:
    graph_props[prop.getName()] = prop

# flag (true/false) a list of nodes (already saved as ub_list) so that they can be aliased later
mvl.markUbiquitousNodes(graph, ub_list)

# create nodes for all species
id2n = mvl.speciesToNodes(graph, cmod, ub_list, graph_props)

# TODO: initial tulip graph creation is slow because of Mimoza Black Box Code!
# TODO: reaction node creation still relies on copied code (problem: reaction node compartments need helper
# TODO: functions for proper assignment -> understand what the copied code does and rewrite!)
########################################################################################################################
# save compartment hierarchy to dictionary "c_id2level" for function import_sbml
# required syntax: c_id2level = {'Cell': (0, ''), 'Boundary': (1, 'Cell')}
comp_list = [''] + cmod.getCompartmentIds()
level_list = [(x, comp_list[x]) for x in range(len(comp_list))]
c_id2level = dict(zip(comp_list[1:], level_list))
# function copied from Mimoza
def compute_c_id2info(c_id2level, input_model):
    c_id2info = {}
    for comp in input_model.compartments:
        c_id = comp.getId()
        c_name = comp.getName()
        if not c_name:
            c_name = c_id
        c_id2info[c_id] = (c_name, '', c_id2level[c_id])
    c_id2outs = {}
    for c_id in c_id2info.keys():
        _, _, (_, out_c_id) = c_id2info[c_id]
        outs = []
        while out_c_id:
            outs.append(out_c_id)
            _, _, (_, out_c_id) = c_id2info[out_c_id]
        c_id2outs[c_id] = outs
    return c_id2info, c_id2outs
c_id2info, c_id2outs = compute_c_id2info(c_id2level, cmod)

# helper function for reaction compartment assignment (copied from Mimoza)
def get_r_comp(all_comp_ids):
    if len(all_comp_ids) == 1:
        return all_comp_ids.pop()
    get_level = lambda c_id: c_id2level[c_id][0]
    outer_most = min(all_comp_ids, key=get_level)
    inner_most = max(all_comp_ids, key=get_level)
    outer_level, inner_level = get_level(outer_most), get_level(inner_most)
    if outer_level == inner_level or (outer_most not in c_id2outs[inner_most]):
        candidates = set(c_id2outs[inner_most]) & set(c_id2outs[outer_most])
        if candidates:
            return max(candidates, key=get_level)
        else:
            return outer_most
    if inner_level - outer_level > 1:
        return max(c_id2outs[inner_most], key=get_level)
    return outer_most

# create nodes for all reactions
mvl.reactionsToNodes(get_r_comp, graph, id2n, cmod)
########################################################################################################################

# clone/alias nodes selected above (in ub_list)
# -> chebi seems to be only important for extracting ubiquitous species ids!
# => provide them by passing a list (ubs), converted to a set
root = sbml2tlp.cloneNodes(set(ub_list), graph)



##############################
#run Mimoza layout functions #
##############################

# actual layouting, creation of graph structure
fc, (n2lo, e2lo), hidden_c_ids, c_id_hidden_ubs, graph = tulip_graph2geojson.graph2geojson(c_id2info, c_id2outs, root)

# save raw Mimoza graph output as backup
tlp.saveGraph(graph, "{}_raw.tlp".format(model_id))


print("Initial Tulip graph created.")



#######################################################
# map sbml element ids to graph nodes in a dictionary #
#######################################################

# Section Description:
# --------------------
# example: iIT341 initial graph (returned by Mimoza): 1039 nodes created
# at this point only 1 node per id if auto aliasing is off!
# check: len(cmod.getSpeciesIds()+cmod.getReactionIds()) returns 1039 as well
# since in the raw graph only 1 node per id exists (if auto aliasing is off)
#
# => create mapping between ids (keys) and nodes (lists of 1) in a dict that reflects the status of the raw graph
# => each new node that gets linked to the same id (via aliasing) will get a new entry with a new key (= pathway name)

# use only if 1 model id = 1 graph node is true (only if auto aliasing is off and only initial, unchanged graph is
# the current graph)
graph_elements = {}

if AUTO_ALIAS == False:
    # all species and all reactions of sbml model = total set of model elements
    model_elements = cmod.species + cmod.reactions
    #dict with keys: model ids and values = dict with keys: pathway names and values = lists of associated graph nodes
    for m_e in model_elements:
        graph_elements[m_e.getId()] = {}
        # use [0] to get first entry of that list so that the returned object is the nodes itself and not a list again
        graph_elements[m_e.getId()]["initial_node"] = [n for n in graph.getNodes() if graph_props["id"][n] ==
                                                       m_e.getId()][0]
        # edges do not exist in SBML model! (and created graph edges have no ids, just names)
        # -> edges are not included in dict!
        # if a node was copied for layouting the obj reaction, it gets saved here (default value: None, gets replaces
        #  with tulip node object that is created for obj re subgraph layout)
        # TODO: we probably don't need the obj re tracking anymore!
        graph_elements[m_e.getId()]["obj_re_copied_node"] = None
        # track if graph element was used for layouting already (default value: False, gets changed to True once the
        # respective node was layouted in a subgraph)
        graph_elements[m_e.getId()]["layout"] = False



########################################################################################################################
# Problem: if a node is shared between pathways, the last valid layouting call will determine its position meaning that
# previous positions of that node in other layouted subgraphs that also contain it are overwritten!!
#
# example in Carlsson model: since PEP and R1 are part of glyc AND exch re, the assigned positions of the PEP
# and R1 nodes after layouting glyc are overwritten by layouting exch re!
#
# -> give subgraphs their own versions of a shared node by copying it from the initial graph created by Mimoza
# -> also: track what graph elements have already been used for layouting to determine what nodes and edges need to
#          be copied
#
# 1.) create graph with initial graph object creation (at this point: 1 id = 1 node only)
# 2.) after graph creation: see if some nodes need to be copied so that each pathway has its own version
#     (two pathways can have the same id entry "PEP" but one will have e.g. PEP 39 and the other one PEP 90 etc.)
# 3.) then when all tulip nodes (for species and reactions) are assembled
#     -> inducedSubGraph, name, layout, save file
########################################################################################################################


# create objective reaction nodes list:
obj_re_nodes = []

for n in graph.getNodes():
    # select obj reaction node
    if graph_props["id"][n] == objective_re_id:
        # add obj reaction node to the list
        obj_re_nodes.append(n)
    # iterate over all nodes which are connected to obj reaction node and add them to the list
    # the new PEP node will be added instead of the old(raw) node because the connection to the raw node was cut and a
        # new connection to the new PEP node was created
        for node in graph.getInOutNodes(n):
            obj_re_nodes.append(node)



#################################
# change basic graph properties #
#################################

# add labels to nodes
# define node types for labeling
NODE_TYPES = {'metabolites': 1,
              'reactions'  : 2
              }

mvl.addLabels(graph, NODE_TYPES['metabolites'], "id")
mvl.addLabels(graph, NODE_TYPES['reactions'], "id")

# TODO: add labels to edges (possible?)
#for e in graph.getEdges():
#    viewLabel[e] = "test"
# -> edge labels are set in "viewLabel" but not visible!

# change node or edge colour
# define colours
COLOURS = {'red'  : tlp.Color(255, 0, 0),
           'green': tlp.Color(0, 255, 0),
           'blue' : tlp.Color(0, 0, 255),
           'white': tlp.Color(255, 255, 255),
           'grey' : tlp.Color(220, 220, 220),
           'black': tlp.Color(0, 0, 0)
           }

# set colour of all edges (complete graph)
mvl.setColours(graph, graph.getEdges(), COLOURS['grey'])

# change line types (Polyline, BezierCurve, CatmullRomCurve, CubicSplineCurve)
for e in graph.getEdges():
    graph_props["viewShape"][e] = tlp.EdgeShape.BezierCurve



##############################################################
# layout pathways in subgraphs (Glycolysis, TCA Cycle, etc.) #
##############################################################

# define function to get subsystem node list from different biological pathways (e.g.: glycolysis, tca) and apply layout
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
        for e in g.getInOutEdges(n):
            if g.source(e) in ex_nodes and g.target(e) in nodes:
                edge_dict[graph_props["id"][n]].update({"out": g.target(e)})
            if g.target(e) in ex_nodes and g.source(e) in nodes:
                edge_dict[graph_props["id"][n]].update({"in": g.source(e)})
    # exclude the already layouted nodes from the main node list
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
        for e_type, node in edge_dict[graph_props["id"][n]].items():
            if e_type == "out":
                sg.addEdge(n, node)
            if e_type == "in":
                sg.addEdge(node, n)
    # layout subgraph
    sg = mvl.layoutSubgraph(sg, sg_layout)
    # flag all the nodes of the glycolysis subgraph as "already used" and add pathway info ( = subgraph name) to their
    # graph properties
    for n in sg.getNodes():
        graph_elements[graph_props["id"][n]]['layout'] = True
        graph_props["pathway"][n] = sg_name
    # add pathway info to egdes aswell
    for e in sg.getEdges():
        graph_props["pathway"][e] = sg_name

    return sg

if LAYOUT_PATHWAYS:

    # TODO: generalize this!

    """
    for path_name, path_ids in path_dict.items():
        path_list = path_ids["metabolites"] + path_ids["reactions"]
        if "ycle" in path_name:
            path_layout = "Circular (OGDF)"
        else:
            path_layout = "Sugiyama (OGDF)"
        path_sg = layoutPathway(graph, graph_elements, graph_props, path_list, path_name, path_layout)
    """

    # layout Glycolysis
    glyc_list = path_dict["Glycolysis"]["metabolites"] + path_dict["Glycolysis"]["reactions"]
    glyc_sg = layoutPathway(graph, graph_elements, graph_props, glyc_list,
                                                        "Glycolysis", "Sugiyama (OGDF)")
    # layout TCA cycle
    tca_list = path_dict["TCA cycle"]["metabolites"] + path_dict["TCA cycle"]["reactions"]
    tca_sg = layoutPathway(graph, graph_elements, graph_props, tca_list,
                                                        "TCA cycle", "Circular (OGDF)")
    # layout Fatty Acid Synthesis
    fas_list = path_dict["Fatty Acid Synthesis"]["metabolites"] + path_dict["Fatty Acid Synthesis"]["reactions"]
    fas_sg = layoutPathway(graph, graph_elements, graph_props, fas_list,
                                                         "Fatty Acid Synthesis", "Sugiyama (OGDF)")
    # layout Urea Cycle
    uc_list = path_dict["Urea Cycle"]["metabolites"] + path_dict["Urea Cycle"]["reactions"]
    uc_sg = layoutPathway(graph, graph_elements, graph_props, uc_list,
                                                         "Urea Cycle", "Circular (OGDF)")
    # layout Aminosugar Metabolism
    am_list = path_dict["Aminosugar Metabolism"]["metabolites"] + path_dict["Aminosugar Metabolism"]["reactions"]
    am_sg = layoutPathway(graph, graph_elements, graph_props, am_list,
                                                         "Aminosugar Metabolism", "Sugiyama (OGDF)")

    # combine pathways in output graph
    # ! if a completely new graph object is created, the properties need to be generated again because they are tied
    # to a graph object (properties of a graph can be used for its subgraphs but thats just because subgraphs are only partial
    # views of a whole graph, not new objects) !

    # create output graph 'G'
    G = tlp.newGraph()
    # create props for new combined output graph 'G'
    G_tup = mvl.create_props(G)
    # save graph props in dictionary
    G_props = {}
    for prop in G_tup:
        G_props[prop.getName()] = prop
    # add name to graph 'G'
    G.setName("output_combined")

    # copy subgraphs to output graph 'G' and recreate subgraphs from selection in graph 'G'
    # TODO: why do i need to use the properties of graph here if i want to change selection status of G nodes?
    graph_props["viewSelection"].setAllNodeValue(False)
    graph_props["viewSelection"].setAllEdgeValue(False)
    #layout Glycolysis
    tlp.copyToGraph(G, glyc_sg, None, graph_props["viewSelection"])
    G_glyc = G.addSubGraph(graph_props["viewSelection"], "Glycolysis")
    G_glyc.setName("Glycolysis")
    graph_props["viewSelection"].setAllNodeValue(False)
    graph_props["viewSelection"].setAllEdgeValue(False)
    # layout TCA cycle
    tlp.copyToGraph(G, tca_sg, None, graph_props["viewSelection"])
    G_tca = G.addSubGraph(graph_props["viewSelection"], "TCA cycle")
    G_tca.setName("TCA cycle")
    graph_props["viewSelection"].setAllNodeValue(False)
    graph_props["viewSelection"].setAllEdgeValue(False)
    #
    tlp.copyToGraph(G, fas_sg, None, graph_props["viewSelection"])
    G_fas = G.addSubGraph(graph_props["viewSelection"], "Fatty Acid Synthesis")
    G_fas.setName("Fatty Acid Synthesis")
    graph_props["viewSelection"].setAllNodeValue(False)
    graph_props["viewSelection"].setAllEdgeValue(False)

    tlp.copyToGraph(G, uc_sg, None, graph_props["viewSelection"])
    G_uc = G.addSubGraph(graph_props["viewSelection"], "Urea Cycle")
    G_uc.setName("Urea Cycle")
    graph_props["viewSelection"].setAllNodeValue(False)
    graph_props["viewSelection"].setAllEdgeValue(False)

    tlp.copyToGraph(G, am_sg, None, graph_props["viewSelection"])
    G_am = G.addSubGraph(graph_props["viewSelection"], "Aminosugar Metabolism")
    G_am.setName("Aminosugar Metabolism")
    graph_props["viewSelection"].setAllNodeValue(False)
    graph_props["viewSelection"].setAllEdgeValue(False)

    #reset selection of G elements
    G_props["viewSelection"].setAllNodeValue(False)
    G_props["viewSelection"].setAllEdgeValue(False)
    # separate overlapping individual pathways
    G.applyLayoutAlgorithm("Connected Component Packing")
    # add pathway labels (invisible nodes with visible label)
    for sg in G.getSubGraphs():
        n = sg.addNode()
        G_props["viewLabel"][n] = sg.getName()
        G_props["viewColor"][n] = tlp.Color(255, 255, 255) # set node color to white (invis on white bg)
        G_props["viewSize"][n] = tlp.Size(5, 5, 1) # normal node = (2.5, 2.5, 1); bigger -> (5, 5, 1)
        G_props["viewFontSize"][n] = 24 # default is 18
    # TODO: control subgraph ( = pathway) positions

    # save output graph
    tlp.saveGraph(G, "{0}_output_combined.tlp".format(model_id))


    print("Pathway layouting finished, results saved in {0}_output_combined.tlp".format(model_id))



#########################################################################
# individually layout the objective reaction and all exchange reactions #
#########################################################################

if LAYOUT_EX_RE_AND_BIOMASS:

    # I.) Exchange Reactions Layouting
    # get all exchange reactions (from ex_re_list, extracted from model with cbmpy if info available)

    #create id list of all boundary (exchange) reactions
    ex_re_list = []
    ex_re_list.extend(cmod.getExchangeReactionIds())

    # if exchange reactions are not annotated in sbml -> look for a exchange reaction group
    if len(ex_re_list) == 0:
        for g_id in cmod.getGroupIds():
            group = cmod.getGroup(g_id)
            if group.name in ["Transport", "Transport, Extracellular", "Exchange", "Exchange Reactions"]:
                ex_re_list.extend(group.getMemberIDs())

    # if objective reaction (e.g. biomass reaction) is in the list of exchange reactions, remove it
    for re in ex_re_list:
        if re == objective_re_id:
            ex_re_list.remove(re)

    # create an initial list of reaction node objects from ids
    ex_re_node_obj = [r for r in graph.getNodes() if graph_props["id"][r] in ex_re_list]

    # check for nodes which have already been used for layouting
    # (with the flag in dict graph_elements[node]['layout'] - boolean value)
    # save the nodes which are excluded because they have already been used for layouting in a list for later use
    ex_re_exclude_nodes = [n for n in ex_re_node_obj if graph_elements[graph_props["id"][n]]['layout'] is True]

    # exclude the already layouted nodes from the initial ex re node list
    for n in list(ex_re_node_obj):
        if n in ex_re_exclude_nodes:
            ex_re_node_obj.remove(n)

    # create subgraph with all exchange reactions
    ex_re_sg = mvl.aliasSeparateReactions(graph, ex_re_list, ex_re_node_obj, "exchange reactions", "Tree Radial",
                                          graph_props)

    # TODO: add new nodes with properties of the excluded ones to the subgraph

    # flag all nodes in exchange reaction subgraph as having been layouted
    for n in ex_re_sg.getNodes():
        graph_elements[graph_props["id"][n]]["layout"] = True

    # get nodes of ex reaction subgraph
    ex_node_list = []
    for node in ex_re_sg.getNodes():
        ex_node_list.append(node)

    #------------------------------------------------------------------------------------------------------------------#

    # II.) Objective Reaction Layouting
    # check obj_re_nodes list for nodes that have already been used for layouting

    # save the nodes which are excluded because they have already been used for layouting in a list for later use
    obj_re_exclude_nodes = [n for n in obj_re_nodes if graph_elements[graph_props["id"][n]]['layout'] is True]

    # exclude the already layouted nodes from the initial ex re node list
    for n in list(obj_re_nodes):
        if n in obj_re_exclude_nodes:
            obj_re_nodes.remove(n)

    # create subgraph for objective reaction
    obj_re_sg = graph.inducedSubGraph(obj_re_nodes)
    obj_re_sg.setName("only objective reaction")

    # TODO: add new nodes with properties of the excluded ones to the subgraph

    # flag all nodes in exchange reaction subgraph as having been layouted
    for n in obj_re_sg.getNodes():
        graph_elements[graph_props["id"][n]]["layout"] = True

    # create combined subgraph with all ex reactions and objective reaction
    obj_re_and_ex_re_list = obj_re_nodes + ex_node_list
    obj_re_and_ex_re_sg = graph.inducedSubGraph(obj_re_and_ex_re_list)
    obj_re_and_ex_re_sg.setName("objective and ex reactions")

    #------------------------------------------------------------------------------------------------------------------#

    # III.) Combine Exchange and Objective Reaction Layouts
    # apply layout to obj reaction and ex reations individually and save combined result
    mvl.layoutSubgraph(ex_re_sg, "Tree Radial", {"layer spacing": 10}, 0, 0, 0)
    # Sugiyama: change node distance (default: 3) and layer distance (default: 3) to scale up the biomass reaction
    mvl.layoutSubgraph(obj_re_sg, "Sugiyama (OGDF)", {"node distance": 8, "layer distance": 8}, 0, 0, 0)
    # in combined graph: the biomass function is not in the center of the inner circle of exchange reactions
    # -> reasons: although both subgraphs are centered around the point (0,0), the center of the ex re tree
    # => radial is not the same point as the center of the inner circle due to appendages on the outside!

    # save combined subgraph
    tlp.saveGraph(obj_re_and_ex_re_sg, "{}_biomass_ex_re.tlp".format(model_id))
    # TODO: put biomass and ex re in new output graph
    # TODO: generalize the copyToGraph function calls so that it can get changed quickly (at the very end)

    print("Biomass and exchange reaction layouting finished, results saved in {}_biomass_ex_re.tlp".format(model_id))



########################
# visualize FBA result #
########################

# color nodes and edges depending on flux value
# input format: dict = {"M1": -3, ...}
# TODO: test if flux overlay also works on the new output graph G, and others
if COLOUR_FBA:
    mvl.colourFBA(graph, fba_results, COLOURS["green"], COLOURS["red"], COLOURS["grey"])



############################################################################
# save final graph, generate SBGN and SBML+Layout Files (Mimoza functions) #
############################################################################

# save total graph with all changes
tlp.saveGraph(graph, "{}_graph.tlp".format(model_id))

# add layout to SBML file
groups_document = reader.readSBML(groups_sbml)
groups_model = groups_document.getModel()
gen_model = False
#sbml_manager.save_as_layout_sbml(groups_model, gen_model, groups_sbml, n2lo, ub_sps)

# create SBGN file
# TODO: check SBGN: is it valid?, what information is stored?
groups_sbgn = os.path.join(cDir, '{}.sbgn'.format(model_id))
# TODO: can we link Mimoza sbgn export function to our graphs?
#sbgn_helper.save_as_sbgn(n2lo, e2lo, groups_model, groups_sbgn)



######################################
# export the graph as GML and/or SVG #
######################################

if EXPORT_GML:
    mvl.exportGML(graph, model_id, "GML_test.tlp", "test.gml")
    print("Graph exported as GML file.")

if EXPORT_SVG:
    # TODO: creates an .svg file but it opening it with Inkscape doesn't show its content
    # get default parameters of SVG export plugin
    params = tlp.getDefaultPluginParameters('SVG Export', graph)
    # export graph as .svg file
    tlp.exportGraph('SVG Export', graph, '{}.svg'.format(model_id), params)
    print("Graph exported as SVG file.")



######################################
# create SQL table from graph object #
######################################


# 1. create necessary data structures

# create a dictionary of dicts for the graph properties to link types of properties to their data:
# dict setup: keys = property names with values = dict with 1. key 'type' => type of property
#                                                       and 2. key 'data' => list of data
#                                                                              (same order as columns in graph object)
#
# for SQLite: convert boolean values to integers (True = 1, False = 0) since sqlite has no boolean data type!
graph_props_node_dict = {}
graph_props_edge_dict = {}

graph_props_ordered_n_list = []
graph_props_ordered_e_list = []

graph_prop_ordered_types = []

for prop in graph.getObjectProperties():
    # add types of props to a list
    graph_prop_ordered_types.append(prop.getTypename())
    # add property values of nodes to list
    prop_node_sublist = []
    for n in graph.getNodes():
        if prop.getTypename() == "bool":
            prop_node_sublist.append(prop[n] * 1) # convert True/False to 1/0
        else:
            prop_node_sublist.append(prop[n])
    # add graph node prop info to ordered list
    graph_props_ordered_n_list.append(prop_node_sublist)
    # fill dict with graph node prop info
    graph_props_node_dict[prop.getName()] = {}
    graph_props_node_dict[prop.getName()]["type"] = prop.getTypename()
    graph_props_node_dict[prop.getName()]["data"] = prop_node_sublist

    # add property values of edges to list
    prop_edge_sublist = []
    for e in graph.getEdges():
        if prop.getTypename() == "bool":
            prop_edge_sublist.append(prop[e] * 1) # convert True/False to 1/0
        else:
            prop_edge_sublist.append(prop[e])
    # add graph edge prop info to ordered list
    graph_props_ordered_e_list.append(prop_edge_sublist)
    # fill dict with graph edge prop info
    graph_props_edge_dict[prop.getName()] = {}
    graph_props_edge_dict[prop.getName()]["type"] = prop.getTypename()
    graph_props_edge_dict[prop.getName()]["data"] = prop_edge_sublist


# 2. encode graph property data to SQLite compatible types

# go through the list of graph property lists:
# since a sublist contains only one graph property for all nodes/edges -> the first element of the sublist is
# exemplary for the whole list
graph_props_node_tuplelist = [] # after encoding is finished: save encoded prop lists in tuples to prepare for INSERT
# add node ids to beginning of tuple list to serve as primary key in table
graph_props_node_tuplelist.append([i for i in range(len(list(graph.getNodes())))])

for prop_type, prop_data_list in zip(graph_prop_ordered_types, graph_props_ordered_n_list):
    # allowed data types, no encoding needed
    if prop_type == 'string':
        pass
    if prop_type == 'int':
        pass
    if prop_type == 'double': # = float (in SQL -> REAL)
        pass
    if prop_type == 'bool': # was converted to INTEGER
        pass
    # data types that require encoding
    # every item of properties with these types is a list itself!
    enc_types = ['size', 'color', 'layout', 'vector<string>']
    # convert data of a type that requires encoding to string
    if prop_type in enc_types:
        for list_item in list(prop_data_list):
            # convert list item to string
            str_item = ','.join(str(item) for item in list_item)
            # remove old list_item
            prop_data_list.remove(list_item)
            # add new string item
            prop_data_list.append(str_item)
    # exlucde props with type 'graph' from tuplelist (viewMetaGraph)
    included = True
    if prop_type == 'graph':
        included = False
    # final step: save included, encoded sublist as tuple
    if included is True:
        graph_props_node_tuplelist.append(tuple(prop_data_list))

graph_props_edge_tuplelist = []
# add edge ids to beginning of tuple list to serve as primary key in table
graph_props_edge_tuplelist.append([i for i in range(len(list(graph.getEdges())))])
for prop_type, prop_data_list in zip(graph_prop_ordered_types, graph_props_ordered_e_list):
    # allowed data types, no encoding needed
    if prop_type == 'string':
        pass
    if prop_type == 'int':
        pass
    if prop_type == 'double': # = float (in SQL -> REAL)
        pass
    if prop_type == 'bool': # was converted to INTEGER
        pass
    # data types that require encoding
    # every item of properties with these types is a list itself!
    enc_types = ['size', 'color', 'layout', 'vector<string>']
    # convert data of a type that requires encoding to string
    if prop_type in enc_types:
        for list_item in list(prop_data_list):
            # convert list item to string
            str_item = ','.join(str(item) for item in list_item)
            # remove old list_item
            prop_data_list.remove(list_item)
            # add new string item
            prop_data_list.append(str_item)
    # exlucde props with type 'graph' from tuplelist (viewMetaGraph)
    included = True
    if prop_type == 'graph':
        included = False
    # final step: save included, encoded sublist as tuple
    if included is True:
        graph_props_edge_tuplelist.append(tuple(prop_data_list))


# 3. create an SQL database with sqlite3
# connect to database
db = sqlite3.connect('{0}_graph.db'.format(model_id))
# get a cursor object
c = db.cursor()
# pass SQL statements into cursor
# create empty graph node table (skips if table already exists)
c.execute('''CREATE TABLE IF NOT EXISTS graph_nodes(node_id INTEGER PRIMARY KEY,
anc_id TEXT,
anc_name TEXT,
anc_term TEXT,
c_id TEXT,
clone_id TEXT,
edge_directionality TEXT,
f TEXT,
fake INTEGER,
flux REAL,
id TEXT,
name TEXT,
pathway TEXT,
rel_c_ids TEXT,
rev INTEGER,
stoich REAL,
term TEXT,
tr INTEGER,
type INTEGER,
ub INTEGER,
viewBorderColor TEXT,
viewBorderWidth REAL,
viewColor TEXT,
viewFont TEXT,
viewFontAwesomeIcon TEXT,
viewFontSize INTEGER,
viewLabel TEXT,
viewLabelBorderColor TEXT,
viewLabelBorderWidth REAL,
viewLabelColor TEXT,
viewLabelPosition INTEGER,
viewLayout TEXT,
viewMetric REAL,
viewRotation REAL,
viewSelection INTEGER,
viewShape INTEGER,
viewSize INTEGER,
viewSrcAnchorShape INTEGER,
viewSrcAnchorSize TEXT,
viewTexture TEXT,
viewTgtAnchorShape INTEGER,
viewTgtAnchorSize TEXT)
''')
# create empty graph edge table (skips if table already exists)
c.execute('''CREATE TABLE IF NOT EXISTS graph_edges(edge_id INTEGER PRIMARY KEY,
anc_id TEXT,
anc_name TEXT,
anc_term TEXT,
c_id TEXT,
clone_id TEXT,
edge_directionality TEXT,
f TEXT,
fake INTEGER,
flux REAL,
id TEXT,
name TEXT,
pathway TEXT,
rel_c_ids TEXT,
rev INTEGER,
stoich REAL,
term TEXT,
tr INTEGER,
type INTEGER,
ub INTEGER,
viewBorderColor TEXT,
viewBorderWidth REAL,
viewColor TEXT,
viewFont TEXT,
viewFontAwesomeIcon TEXT,
viewFontSize INTEGER,
viewLabel TEXT,
viewLabelBorderColor TEXT,
viewLabelBorderWidth REAL,
viewLabelColor TEXT,
viewLabelPosition INTEGER,
viewLayout TEXT,
viewMetric REAL,
viewRotation REAL,
viewSelection INTEGER,
viewShape INTEGER,
viewSize INTEGER,
viewSrcAnchorShape INTEGER,
viewSrcAnchorSize TEXT,
viewTexture TEXT,
viewTgtAnchorShape INTEGER,
viewTgtAnchorSize TEXT)
''')
# reorder tuple list for INSERT INTO statement:
# instead of [(all prop1), (all prop2), ...] I want [(prop1[1], prop2[1], ...), (prop1[2], prop2[2], ...), ...]
graph_props_node_insertlist = []
graph_props_edge_insertlist = []
for i in range(len(list(graph.getNodes()))):
    graph_props_node_insertlist.append([tuple[i] for tuple in graph_props_node_tuplelist])
for i in range(len(list(graph.getEdges()))):
    graph_props_edge_insertlist.append([tuple[i] for tuple in graph_props_edge_tuplelist])
# fill tables with encoded values
for tup in graph_props_node_insertlist:
    c.execute('''INSERT INTO graph_nodes VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,
    ?,?,?,?,?,?)''', tup)
for tup in graph_props_edge_insertlist:
    c.execute('''INSERT INTO graph_edges VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,
    ?,?,?,?,?,?)''', tup)
# commit changes to database
db.commit()
# close connection to database
db.close()


print("Graph data has been stored in SQL database.")



#####################################
# retrieve graphs from SQL database #
#####################################
# TODO: after putting in graph info try to retireve it again and save it as new tlp (see if layout stay constant)
# TODO: use biomass layout for testing purposes


# fetch graph properties of nodes and columns from database
db_return_node_props = []

# TODO: decode graph props again after getting them from the database












