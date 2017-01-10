# Test transition from sbml file to structured tulip graph, sbml with layout and sbgn file
# necessary  dependencies (Wind64 Anaconda): pip install geojson tarjan libsbgnpy tulip-python
# for SQL database: install sqlite 3 (https://www.tutorialspoint.com/sqlite/sqlite_installation.htm)

import os
# get current directory
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))

import libsbml
import cbmpy
from tulip import tlp
import sqlite3

# necessary Mimoza modules
# graph visualization
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

# TODO: let this refer to the generalized version of the layouting function (instead of glyc + tca specific)
PATHWAY_ALIAS = True
LAYOUT_GLYC_TCA = True # needs PATHWAY_ALIAS true, so that necessary node lists are generated (for pathways)

# layout the exchange reactions and biomass seperately (results in subgraphs: only ex re, only biomass and both)
LAYOUT_EX_RE_AND_BIOMASS = True

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
glyc_sps = []
tca_sps = []
glyc_re = []
tca_re = []
small_sps = []

# -> IF GROUPS DON'T EXIST: provide pathway lists manually (example: e coli carlson model):
if PATHWAY_ALIAS and sbml == os.path.join(cDir, "TestModels", "Ecoli_Carlson2003.xml"):
    # create pathway dictionary (key: pathway name - value: tupel of two lists)
    # path_dict = {"glycolysis": ([list of metabolite nodes],[list of reaction nodes]), "TCA Cycle": ([metabolite
    # nodes],[reaction nodes]), ...}
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

# -> IF GROUPS EXIST: fetch groups from sbml (this doesn't use neither ub_sps_list nor ub_re_list!, specific to iIT341)
# create dictionary of groups where each group name is a key and the values are again dicts with three keys: 'id',
# 'reactions' and 'metabolites'
if PATHWAY_ALIAS and sbml == os.path.join(cDir, "TestModels", "iIT341.xml"):
    # create a dictionary to store model groups (key: group name - value: group id)
    gids = {}
    for g in cmod.groups:
        if g.getName() == 'Glycolysis':
            gids['Glycolysis'] = {}
            gids['Glycolysis']['id'] = g.getId()
            gids['Glycolysis']['metabolites'] = []
            # in iIT341 (and many other models): groups only contain reaction ids (not metabolites)
            gids['Glycolysis']['reactions'] = g.getMemberIDs()
        elif g.getName() == 'TCA cycle':
            gids['TCA cycle'] = {}
            gids['TCA cycle']['id'] = g.getId()
            gids['TCA cycle']['metabolites'] = []
            gids['TCA cycle']['reactions'] = g.getMemberIDs()
    # get all member species per group for lists of reaction nodes
    glyc_re = gids['Glycolysis']['reactions']
    tca_re = gids['TCA cycle']['reactions']
    # for each reaction get the associated species for lists of metabolite nodes
    glyc_sps = []
    tca_sps = []
    for r in glyc_re:
        glyc_sps.extend(cmod.getReaction(r).getSpeciesIds())
    for r in tca_re:
        tca_sps.extend(cmod.getReaction(r).getSpeciesIds())
    glyc_sps = list(set(glyc_sps))
    tca_sps = list(set(tca_sps))
    # add species lists to dict
    # TODO: do that with update function instead?
    gids['Glycolysis']['metabolites'] = glyc_sps
    gids['TCA cycle']['metabolites'] = tca_sps
    # get ubiquitous small metabolites by filtering for high connectivity (> 30)
    small_sps = [s.getId() for s in cmod.species if len(s.isReagentOf()) > 30]
    # remove the small metabolites from the pathways
    for ubs in small_sps:
        if ubs in glyc_sps:
            glyc_sps.remove(ubs)
        if ubs in tca_sps:
            tca_sps.remove(ubs)



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
# special manipulation of Glycolysis and TCA Cycle subgraphs #
##############################################################

# TODO: generalize this!

# get subsystem node list from different biological pathways (e.g.: glycolysis, tca) and apply layout
if LAYOUT_GLYC_TCA and len(glyc_re) > 1 and len(tca_re) > 1:

    # I.) initial subgraph generation
    # generate (initial) list of nodes for pathway elements by going through lists containing selected ids
    # this list will also include nodes that already have been used for layouting
    glyc_nodes = [n for n in graph.getNodes() if graph_props["id"][n] in glyc_sps + glyc_re]

    # check for nodes which have already been used for layouting (with the flag in dict graph_elements[node][
    # 'layout'] - boolean value)
    # save the nodes which are excluded because they have already been used for layouting in a list for later use
    glyc_ex_nodes = [n for n in glyc_nodes if graph_elements[graph_props["id"][n]]['layout'] is True]

    # exclude the already layouted nodes from the main node list: glyc_nodes
    for n in glyc_nodes:
        if n in glyc_ex_nodes:
            glyc_nodes.remove(n)

    # create subgraph with inducedSubgraph function (automatically add all existing edges between seletced nodes aswell)
    # this graph will not have any nodes that have been used for layouting before
    glyc_sg = graph.inducedSubGraph(glyc_nodes)
    # name subgraph
    glyc_sg.setName("glycolysis")

    # if any nodes have been excluded due to the fact that they already have been used for layouting
    # -> we need to fetch copies of these nodes (create new nodes and equip them with all of the properties of the
    # excluded ones)
    # => they also need to have all the edges connected to them! (as it was the case before)
    # use tulip function coptyToGraph, that takes in a selection of graph elements to add to another graph
    # first, reset selection
    graph_props["viewSelection"].setAllNodeValue(False)
    graph_props["viewSelection"].setAllEdgeValue(False)

    # select nodes that have been excluded
    for n in glyc_ex_nodes:
        graph_props["viewSelection"][n] = True

    # select all edges of the excluded nodes that are connected to a glyc node
    # TODO: adding the copied nodes, edges... from the initial graph to the subgraph just makes them float in
    # TODO: empty space next to the network! -> connect the copied parts to the rest
    # TODO: maybe write a function that detects equal reaction nodes (copied nodes) and merges them so that in the end
    # TODO: there is only one reaction node with the combined sum of all In and Out edges
    for e in graph.getEdges():
        if graph.source(e) in glyc_ex_nodes and graph.target(e) in glyc_nodes or graph.source(e) in glyc_nodes and \
                graph.target(e) in glyc_ex_nodes:
            graph_props["viewSelection"][e] = True

    # copy nodes that have been excluded from graph to glyc subgraph + their edges aswell
    # tlp.copyToGraph(outGraph, inGraph[, inSelection = None, outSelection = None])
    tlp.copyToGraph(glyc_sg, graph, graph_props["viewSelection"])

    # now, reset selection again
    graph_props["viewSelection"].setAllNodeValue(False)
    graph_props["viewSelection"].setAllEdgeValue(False)

    # layout subgraph
    mvl.layoutSubgraph(glyc_sg, "Sugiyama (OGDF)", None, 0, 0, 0)

    # save graph object with new subgraphs
    tlp.saveGraph(graph, "{}_graph_test.tlp".format(model_id))

    # flag all the nodes of the glycolysis subgraph as "already used"
    for n in glyc_sg.getNodes():
        graph_elements[graph_props["id"][n]]['layout'] = True



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

    # X.) LAYOUTING STATUS CHECK 1: ex re nodes -> identify and exclude already used nodes
    # create an initial list of reaction node objects from ids
    ex_re_node_obj = [r for r in graph.getNodes() if graph_props["id"][r] in ex_re_list]
    # check for nodes which have already been used for layouting
    # (with the flag in dict graph_elements[node]['layout'] - boolean value)
    # save the nodes which are excluded because they have already been used for layouting in a list for later use
    ex_re_exclude_nodes = [n for n in ex_re_node_obj if graph_elements[graph_props["id"][n]]['layout'] is True]
    # exclude the already layouted nodes from the initial ex re node list
    for n in ex_re_node_obj:
        if n in ex_re_exclude_nodes:
            ex_re_node_obj.remove(n)

    # if any nodes have been excluded due to the fact that they already have been used for layouting
    # -> we need to fetch copies of these nodes (create new nodes and equip them with all of the properties of the
    # excluded ones)
    # => they also need to have all the edges connected to them! (as it was the case before)
    # use tulip function coptyToGraph, that takes in a selection of graph elements to add to another graph
    # first, reset selection
    graph_props["viewSelection"].setAllNodeValue(False)
    graph_props["viewSelection"].setAllEdgeValue(False)
    # select nodes that have been excluded
    for n in ex_re_exclude_nodes:
        graph_props["viewSelection"][n] = True
    # select all edges of the excluded nodes that are connected to a ex_re node
    for e in graph.getEdges():
        if graph.source(e) in ex_re_exclude_nodes and graph.target(e) in ex_re_node_obj \
        or graph.source(e) in ex_re_node_obj and graph.target(e) in ex_re_exclude_nodes:
            graph_props["viewSelection"][e] = True

    # create a subgraph with all the ex reaction nodes and the all biomass reaction nodes together
    # create subgraph with all exchange reactions
    ex_re_sg = mvl.aliasSeparateReactions(graph, ex_re_list, ex_re_node_obj, "exchange reactions", "Tree Radial",
                                          graph_props)

    # X.) LAYOUTING STATUS CHECK 2: ex re nodes -> fetch copies of excluded nodes from initial graph and add them
    # copy nodes that have been excluded from graph to ex re subgraph + their edges aswell
    # tlp.copyToGraph(outGraph, inGraph[, inSelection = None, outSelection = None])
    tlp.copyToGraph(ex_re_sg, graph, graph_props["viewSelection"])
    # now, reset selection again
    graph_props["viewSelection"].setAllNodeValue(False)
    graph_props["viewSelection"].setAllEdgeValue(False)
        # flag all nodes in exchange reaction subgraph as having been layouted
    for n in ex_re_sg.getNodes():
        graph_elements[graph_props["id"][n]]["layout"] = True

    # get nodes of ex reaction subgraph
    ex_node_list = []
    for node in ex_re_sg.getNodes():
        ex_node_list.append(node)

    #------------------------------------------------------------------------------------------------------------------#

    # II.) Objective Reaction Layouting
    # X.) LAYOUTING STATUS CHECK 1: obj re nodes -> identify and exclude already used nodes
    # check obj_re_nodes list for nodes that have already been used for layouting
    # save the nodes which are excluded because they have already been used for layouting in a list for later use
    obj_re_exclude_nodes = [n for n in obj_re_nodes if graph_elements[graph_props["id"][n]]['layout'] is True]
    # exclude the already layouted nodes from the initial ex re node list
    for n in list(obj_re_nodes):
        if n in obj_re_exclude_nodes:
            obj_re_nodes.remove(n)
    # if any nodes have been excluded due to the fact that they already have been used for layouting
    # -> we need to fetch copies of these nodes (create new nodes and equip them with all of the properties of the
    # excluded ones)
    # => they also need to have all the edges connected to them! (as it was the case before)
    # use tulip function coptyToGraph, that takes in a selection of graph elements to add to another graph
    # first, reset selection
    graph_props["viewSelection"].setAllNodeValue(False)
    graph_props["viewSelection"].setAllEdgeValue(False)
    # select nodes that have been excluded
    for n in obj_re_exclude_nodes:
        graph_props["viewSelection"][n] = True
    # select all edges of the excluded nodes that are connected to a ex_re node
    for e in graph.getEdges():
        if graph.source(e) in obj_re_exclude_nodes and graph.target(e) in obj_re_nodes \
        or graph.source(e) in obj_re_nodes and graph.target(e) in obj_re_exclude_nodes:
            graph_props["viewSelection"][e] = True

    # create subgraph for objective reaction
    obj_re_sg = graph.inducedSubGraph(obj_re_nodes)
    obj_re_sg.setName("only objective reaction")

    # X.) LAYOUTING STATUS CHECK 2: obj re nodes -> fetch copies of excluded nodes from initial graph and add them
    # copy nodes that have been excluded from graph to ex re subgraph + their edges aswell
    # tlp.copyToGraph(outGraph, inGraph[, inSelection = None, outSelection = None])
    tlp.copyToGraph(obj_re_sg, graph, graph_props["viewSelection"])
    # now, reset selection again
    graph_props["viewSelection"].setAllNodeValue(False)
    graph_props["viewSelection"].setAllEdgeValue(False)
    # flag all nodes in exchange reaction subgraph as having been layouted
    for n in obj_re_sg.getNodes():
        graph_elements[graph_props["id"][n]]["layout"] = True

    # create combined subgraph with all ex reactions and objective reaction
    obj_re_and_ex_re_list = obj_re_nodes + ex_node_list
    obj_re_and_ex_re_sg = graph.inducedSubGraph(obj_re_and_ex_re_list)
    obj_re_and_ex_re_sg.setName("objective and ex reactions")

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



########################
# visualize FBA result #
########################

# color nodes and edges depending on flux value
# input format: dict = {"M1": -3, ...}
if COLOUR_FBA:
    mvl.colourFBA(graph, fba_results, COLOURS["green"], COLOURS["red"], COLOURS["grey"])



############################################################################
# save final graph, generate SBGN and SBML+Layout Files (Mimoza functions) #
############################################################################

# save total graph with all changes
tlp.saveGraph(graph, "{}_graph.tlp".format(model_id))

# TODO: export tulip graph as SVG file
#params = tlp.getDefaultPluginParameters('SVG Export', graph)
#tlp.exportGraph('SVG Export', graph, '{}.svg'.format(model_id), params
# NOT WORKING! -> so far: Exception: No plugin named  SVG Export.
# emailed tulip staff about it

# add layout to SBML file
groups_document = reader.readSBML(groups_sbml)
groups_model = groups_document.getModel()
gen_model = False
#sbml_manager.save_as_layout_sbml(groups_model, gen_model, groups_sbml, n2lo, ub_sps)

# create SBGN file
groups_sbgn = os.path.join(cDir, '{}.sbgn'.format(model_id))
#sbgn_helper.save_as_sbgn(n2lo, e2lo, groups_model, groups_sbgn)



######################################
# create SQL table from graph object #
######################################

# create a list with lists of graph property data which carry the same data as one column in tulip graph object
# order of list = order of nodes in tulip graph object table
# for SQLite: convert boolean values to integers (True = 1, False = 0) since sqlite has no boolean data type!

graph_props_node_data = []
graph_props_edge_data = []

# add node and edge tulip ids to lists
graph_props_node_data.append([str(n) for n in graph.getNodes()])
graph_props_edge_data.append([str(e) for e in graph.getEdges()])

for prop in graph.getObjectProperties():
    # add property values to node list
    prop_node_sublist = []
    for n in graph.getNodes():
        if prop.getTypename() == "bool":
            prop_node_sublist.append(prop[n] * 1)
        else:
            prop_node_sublist.append(prop[n])
    graph_props_node_data.append(prop_node_sublist)

    # add property values to edge list
    prop_edge_sublist = []
    for e in graph.getEdges():
        if prop.getTypename() == "bool":
            prop_edge_sublist.append(prop[e] * 1)
        else:
            prop_edge_sublist.append(prop[e])
    graph_props_edge_data.append(prop_edge_sublist)

# create an SQL database with sqlite3
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
rel_c_ids TEXT,
rev INTEGER,
stoich REAL,
term TEXT,
tr INTEGER,
type INTEGER,
ub INTEGER,
viewBorderColor BLOB,
viewBorderWidth REAL,
viewColor BLOB,
viewFont TEXT,
viewFontAwesomeIcon TEXT,
viewFontColor INTEGER,
viewLabel TEXT,
viewLabelBorderColor BLOB,
viewLabelBorderWidth REAL,
viewLabelColor BLOB,
viewLabelPosition INTEGER,
viewLayout BLOB,
viewMetaGraph BLOB,
viewMetric REAL,
viewRotation REAL,
viewSelection INTEGER,
viewShape INTEGER,
viewSize BLOB,
viewSrcAnchorShape INTEGER,
viewSrcAnchorSize BLOB,
viewTexture TEXT,
viewTgtAnchorShape INTEGER,
viewTgtAnchorSize BLOB)
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
rel_c_ids TEXT,
rev INTEGER,
stoich REAL,
term TEXT,
tr INTEGER,
type INTEGER,
ub INTEGER,
viewBorderColor BLOB,
viewBorderWidth REAL,
viewColor BLOB,
viewFont TEXT,
viewFontAwesomeIcon TEXT,
viewFontColor INTEGER,
viewLabel TEXT,
viewLabelBorderColor BLOB,
viewLabelBorderWidth REAL,
viewLabelColor BLOB,
viewLabelPosition INTEGER,
viewLayout BLOB,
viewMetaGraph BLOB,
viewMetric REAL,
viewRotation REAL,
viewSelection INTEGER,
viewShape INTEGER,
viewSize BLOB,
viewSrcAnchorShape INTEGER,
viewSrcAnchorSize BLOB,
viewTexture TEXT,
viewTgtAnchorShape INTEGER,
viewTgtAnchorSize BLOB)
''')

# TODO: fill tables with values in graph_props_node_data and graph_props_edge_data

# commit changes to database
db.commit()
# close connection to database
db.close()
