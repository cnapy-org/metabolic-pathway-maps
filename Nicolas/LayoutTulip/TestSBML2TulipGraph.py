# Test transition from sbml file to structured tulip graph, sbml with layout and sbgn file
# necessary  dependencies (Wind64 Anaconda): pip install geojson tarjan libsbgnpy tulip-python

import os
# get current directory
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))

import random
import libsbml
import cbmpy
from tulip import tlp

# necessary Mimoza modules
import sbml2tlp
import sbml_manager
import tulip_graph2geojson
import sbgn_helper

# module with own functions
import MVLayout as mvl



#############################################
# read SBML with libsbml and cbmpy, do FBA  #
#############################################

# read SBML file (with libsbml)
#sbml = os.path.join(cDir, "TestModels", "toy_model_biofuel.xml") # 4, 4, False
#sbml = os.path.join(cDir, "TestModels", "e_coli_core.xml") # 6, 4, True
sbml = os.path.join(cDir, "TestModels", "Ecoli_Carlson2003.xml") # 6, 4, True
#sbml = os.path.join(cDir, "TestModels", "Hpylori_iIT341.xml") # 6, 4, True
#sbml = os.path.join(cDir, "TestModels", "iJO1366.xml") # 6, 4, True
#sbml = os.path.join(cDir, "TestModels", "Lactococcus_lactis_MG1363") # ?, ?, set it to 'True' probably

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

# analyze model (FBA)
cbmpy.doFBA(cmod)
fba_results = {}
for r in cmod.reactions:
    fba_results[r.getId()] = round(r.getValue(), 6)



################################################################################
# extract compartment hierarchy, select nodes for aliasing (different methods) #
################################################################################

# save compartment hierarchy to dictionary "c_id2level" for function import_sbml
# c_id2level needs to be automatically extracted from input SBML file
# required form: c_id2level = {'Cell': (0, ''), 'Boundary': (1, 'Cell')}
comp_list = [''] + cmod.getCompartmentIds()
level_list = [(x, comp_list[x]) for x in range(len(comp_list))]
c_id2level = dict(zip(comp_list[1:], level_list))


# 1.) create alias node list based on connectivity <-> user defined thresholds (strategy 1)

#ub_sp_list, ub_re_list = identifyConnectedNodes(cmod, 4, 4, False)


# 2.) create alias node list based on biological pathways <-> group info in sbml - if available (strategy 2)

# create lists of species ids and reaction ids for each pathway (not necessary if provided by sbml)

# specific to carlsson model
# lists of metabolite nodes
glyc_sps = ["GLU_6_P", "FRU_6_P", "FRU_BIS_P", "DHAP", "GA_3P", "PG", "PEP"]
tca_sps = ["CITRATE", "ISOCIT", "AKG", "SUCC_CoA", "SUCC", "FUMARATE", "MALATE", "OXALO"]
# lists of reaction nodes
glyc_re = ["R1", "R2r", "R3", "R4", "R5r", "R6r", "R7r", "R8r", "R9", "RR9"]
tca_re = ["R23r", "R24", "R25", "R26r", "R27r", "R28r", "R29r", "R22"]
# ubiquitous small metabolites
small_sps = ["ATP", "ADP", "NAD", "NADH", "CO2", "NH3"]

not_spec_sps = [sps for sps in cmod.getSpeciesIds() if sps not in glyc_sps or sps not in tca_sps]
not_spec_re = [re for re in cmod.getReactionIds() if re not in glyc_re or re not in tca_re]

# create groups (not necessary if provided by sbml)
cmod.createGroup('not_specified')
path0 = cmod.getGroup('not_specified')
path0.setName('Not specified')
path0.addMember([cmod.getSpecies(s) for s in not_spec_sps] + [cmod.getReaction(r) for r in not_spec_re])

cmod.createGroup('glycolysis')
path1 = cmod.getGroup('glycolysis')
path1.setName('Glycolysis')
path1.addMember([cmod.getSpecies(s) for s in glyc_sps] + [cmod.getReaction(r) for r in glyc_re])

cmod.createGroup('tca_cycle')
path2 = cmod.getGroup('tca_cycle')
path2.setName('TCA Cycle')
path2.addMember([cmod.getSpecies(s) for s in tca_sps] + [cmod.getReaction(r) for r in tca_re])

# automatically create list of boundary species (those which connect the biol. pathway to rest of network)
#glyc_sps_bound = identifyBoundarySpecies(glyc_sps, glyc_re, cmod)
tca_sps_bound_tup = mvl.identifyBoundarySpecies(tca_sps, tca_re, cmod)
# put only the first item of the tuples (species) in a new list and only if it is not part of the glycolysis list
# we don't want to alias to much later on
tca_sps_bound = list(set([tup[0] for tup in tca_sps_bound_tup if tup[0] not in glyc_sps]))

# select biomass reaction for aliasing
#for r in cmod.getActiveObjectiveReactionIds():
    # add biomass reaction with substrates and products to alias lists (ub_sp_list, ub_re_list)?
    # -> commented because some species of tca cycle are also part of biomass (eg: AKG)
    # => aliasing them would break up the tca circle layout
    #ub_sp_list += [s for s in cmod.getReaction(r).getSpeciesIds() if s not in ub_sp_list]

ub_sp_list = tca_sps_bound
ub_re_list = []


# 3.) create alias node list based on a certain compartment <-> compartment info in sbml (strategy 3)

#comp_ids = cmod.getCompartmentIds() # get a list of all compartments
#comp_2 = cmod.getCompartment(comp_ids[1]) # get  compartment 2 (first entry in list of all compartments)
#comp_2_species = comp_2.containsSpecies() # list of species in compartment 2
#comp_2_reac = comp_2.containsReactions() # list of reactions in compartment 2
#ub_sp_list = comp_2_species
#ub_re_list = comp_2_reac


# 4.) manually provide alias node list <-> hard coded by user (for testing purposes)

# specific to toy_model_biofuel
#ub_sp_list = ['ATP', 'ADP', 'M8', 'M10', 'M11', 'M12', 'M13', 'biomass', 'S', 'fuel1', 'fuel2']
#ub_re_list = ['R09', 'R12', 'R13', 'R14', 'R15', 'R16', 'R17', 'R18']


# 5.) create alias node list of all boundary (exchange) reactions (for aliasNode() function, NOT Mimoza!)

ex_re_list = cmod.getExchangeReactionIds()
# ub_sp_list = []
# ub_re_list = []

# combine alias species and alias reaction lists (either for import_sbml or aliasNodes function)
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



#######################################################
# create Tulip graph and run Mimoza layout functions  #
#######################################################

# create new tulip graph
graph = tlp.newGraph()
graph.setName(model.getId())

# create graph properties for new tulip graph
mvl.create_props(graph)

# flag (true/false) a list of nodes (already saved as ub_list) so that they can be aliased later
mvl.markUbiquitousNodes(graph, ub_list)

# create nodes for all model species
id2n = mvl.speciesToNodes(graph, model, ub_list)

# run Mimoza functions for further graph manipulaption
# species to nodes, reactions to nodes, alias node cloning
root, c_id2info, c_id2outs, ub_sps = sbml2tlp.import_sbml(model, groups_sbml, ub_list, c_id2level, graph, id2n)

# actual layouting, creation of graph structure
fc, (n2lo, e2lo), hidden_c_ids, c_id_hidden_ubs, graph = tulip_graph2geojson.graph2geojson(c_id2info, c_id2outs, root)

# save raw Mimoza graph output as backup
tlp.saveGraph(graph, "{}_raw.tlp".format(model_id))



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

# TODO: add labels to edges
#viewLabel = graph.getStringProperty("viewLabel")
#for e in graph.getEdges():
#    viewLabel[e] = "test"
# -> edge labels are set in "viewLabel" but not visible!

# change node or edge colour
# define colours
COLOURS = {'red'  : tlp.Color(255, 0, 0),
           'green': tlp.Color(0, 255, 0),
           'blue' : tlp.Color(0, 0, 255),
           'grey' : tlp.Color(220, 220, 220),
           'black': tlp.Color(0, 0, 0)
           }

# set colour of all edges (complete graph)
mvl.setColours(graph, graph.getEdges(), COLOURS['grey'])

# change line types (Polyline, BezierCurve, CatmullRomCurve, CubicSplineCurve)
viewShape = graph.getIntegerProperty("viewShape")
for e in graph.getEdges():
    viewShape[e] = tlp.EdgeShape.BezierCurve



########################################################################
# apply different layouts on subsystems (compartments, pathways, etc.) #
########################################################################


# 1.) manually define subsystem node lists and apply layout
# specific to toy_model_biofuel
small_circle = ["M1", "M2", "M9", "R02", "R10", "R11"]
larger_circle = ["M3", "M4", "M5", "M6", "M7", "R04", "R05", "R06", "R07", "R08"]

if sbml == os.path.join(cDir, "TestModels", "toy_model_biofuel.xml"):
    # create toymodel subgraphs and apply centered layout
    small_circle_subgraph = mvl.createSubgraph(graph, "small_circle", small_circle, None)
    larger_circle_subgraph = mvl.createSubgraph(graph, "larger_circle", larger_circle, None)
    mvl.layoutSubgraph(small_circle_subgraph, "Circular (OGDF)", 100, 100, 0)
    mvl.layoutSubgraph(larger_circle_subgraph, "Circular (OGDF)", -50, 25, 0)
    # save created toymodel subgraphs
    tlp.saveGraph(small_circle_subgraph, "{}_SmallCircle.tlp".format(model_id))
    tlp.saveGraph(larger_circle_subgraph, "{}_LargerCircle.tlp".format(model_id))


# 2.) get subsystem node list from all species of a certain compartment and apply layout
#comp_ids = cmod.getCompartmentIds() # get a list of all compartments
#comp_2 = cmod.getCompartment(comp_ids[1]) # get  compartment 2 (second entry in list of all compartments)
#comp_2_species = comp_2.containsSpecies() # list of species in compartment 2

# layout all species from compartment 2 (e.g.: extracellular) seperately (i.e.: along the edges of the graph)
#comp_2_sg = mvl.createSubgraph(graph, "comp_2", comp_2_species, None)
#mvl.layoutSubgraph(comp_2_sg, "Hierarchical Graph", 0, 0, 0)

"""
# 3.) get subsystem node list from different biological pathways (e.g.: glycolysis, tca) and apply layout

# I.) initial subgraph generation
# either: pass list of metabolite nodes to function create subgraph
#glyc_sg = mvl.createSubgraph(graph, "glycolysis", glyc_sps, {"distance": 1L})
#tca_sg = mvl.createSubgraph(graph, "tca_cycle", tca_sps, {"distance": 1L})
# or: pass list of reaction nodes to function create subgraph
glyc_sg = mvl.createSubgraph(graph, "glycolysis", glyc_re, {"distance": 1L})
tca_sg = mvl.createSubgraph(graph, "tca_cycle", tca_re, {"distance": 1L})
# set colours of tca and glycolysis
mvl.setColours(graph, glyc_sg.getNodes(), COLOURS['red'])
mvl.setColours(graph, tca_sg.getNodes(), COLOURS['blue'])
# layout tca and glycolysis differently and set initial coordinates for tca and glycolysis
#mvl.layoutSubgraph(glyc_sg, "Sugiyama (OGDF)", {"node distance": 1, "layer distance": 10}, -45, 150, 0)
mvl.layoutSubgraph(glyc_sg, "Hierarchical Graph", {"orientation": "vertical", }, -45, 150, 0)
mvl.layoutSubgraph(tca_sg, "Circular (OGDF)", None, -45, -30, 0)

# II.) manipulate pathway positions relative to each other
# -> choose coordinates so that the two pathways don't overlap:
# -> known assumption: tca cycle is always below glycolysis
# -> move tca cycle!
viewLayout = graph.getLayoutProperty("viewLayout")
selection = graph.getBooleanProperty("viewSelection")

# get y-coord of lowest glycolysis node
bottom_gly_node = None
y_bottom_gly = 10000
for n in glyc_sg.getNodes():
    y = viewLayout[n][1]
    if y < y_bottom_gly:
        y_bottom_gly = y
        bottom_gly_node = n
# get y-coord of highest tca node
y_top_tca = -10000
for n in tca_sg.getNodes():
    y= viewLayout[n][1]
    if y > y_top_tca:
        y_top_tca = y

# get x value of node with lowest y-value (to center tca directly under glycolysis)
x_trans = viewLayout[bottom_gly_node][0]
# get difference for y translation (to avoind overlap between tca and glycolysis)
y_diff = y_top_tca - y_bottom_gly

# center tca at new position: x = x-value of lowest glyc node, y = old y - y_diff - buffer, z = 0
#buff = 60 # buffer space between glycolysis and tca
#tlp.LayoutProperty.center(viewLayout, tlp.Coord(x_trans, -30-y_diff-buff, 0), tca_sg)

# -> try moving the glycolysis instead => edges are still bent (are not moved!)
#buff = 60 # buffer space between glycolysis and tca
#tlp.LayoutProperty.center(viewLayout, tlp.Coord(x_trans, 150+y_diff+buff, 0), glyc_sg)

# ***********************************bounding box idea***********************************
# create a glycolysis bounding box from selection of all glycolysis nodes and edges
# reset selection
selection.setAllNodeValue(False)
selection.setAllEdgeValue(False)
for n in glyc_sg.getNodes():
    selection[n] = True
for e in glyc_sg.getEdges():
    selection[e] = True
# create BB
glyc_bb = tlp.computeBoundingBox(glyc_sg, selection)
# move BB
buff = 60 # buffer space between glycolysis and tca
tlp.BoundingBox.translate(glyc_bb, tlp.Coord(x_trans, 150+y_diff+buff, 0))
# -> Bounding Box does not move with graph and layout can't be confined to a bounding box!

# III.) generate combined graph of both pathways
# show both tca and glycolysis in same graph

# reset selection
selection.setAllNodeValue(False)
selection.setAllEdgeValue(False)
# select tca and glyc nodes
for n in glyc_sg.getNodes():
    selection[n] = True
for n in tca_sg.getNodes():
    selection[n] = True
# select tca and glyc edges
for e in glyc_sg.getEdges():
    selection[e] = True
for e in tca_sg.getEdges():
    selection[e] = True
# create subgraph with only tca and glycolysis from selection
tca_glyc_sg = graph.addSubGraph(selection, "tca and glyc")
# reset selection
selection.setAllNodeValue(False)
selection.setAllEdgeValue(False)

# TEST: call layout_components (from Mimoza)
#mvl.layout_components(tca_glyc_sg)
#-> also struggles with edge movement when called after own graph manipulation!

# ***********************************overlap removal idea***************************************
# call overlap removal on combined tca+glyc subgraph
tca_glyc_sg.applyLayoutAlgorithm('Fast Overlap Removal')
# -> only resolves overlapping nodes, edges still overlap

# ***********************************remove edge crossing idea***********************************
# remove edge crossings between the pathways:
#  tlp.computeLinesIntersection(line1, line2)
# ...

# fix position of a specific node (M1 (id 2) for small circle, M3 (id 8) for larger circle)
# specific to toy_model_biofuel
#viewLayout[tlp.node(2)] = tlp.Coord(50,75,0)
#viewLayout[tlp.node(8)] = tlp.Coord(-50,25,0)
#viewLayout = graph.getLayoutProperty("viewLayout")
#viewLayout[tlp.node(16)] = tlp.Coord(100, 100, 0)
"""

####################################################
# aliasing examples (using the function aliasNode) #
####################################################


# example 1: remove edge (edge 49) connecting biomass to TCA cycle at AKG (node 4)
#            get the edge that connects AKG (node 4) to biomass (R70, node 60)

# specific to carlsson_model
#aliasNode(graph, tlp.node(4), "R70", 22, 27)


# The code below works and only needs a Tulip graph
# take biomass reaction from a graph, alias it and all of its connected species and layout it individually
"""
# example 2: cut out the biomass function (individually alias all of its connected species)
#            and move it together with all of its now aliased species somewhere else
#            biomass function = R70 (node 60)

viewName = graph.getStringProperty("name")
viewLayout = graph.getLayoutProperty("viewLayout")

# specific to carlsson model
biomass_re_id = "R70"

# alias the biomass reaction
biomass_node_list = []
for n in graph.getNodes():
    # select biomass reaction node
    if viewName[n] == biomass_re_id:
        # add biomass reaction node to list
        biomass_node_list.append(n)
        # access biomass reaction node coords (to layout connected nodes around it)
        x, y, z = viewLayout[n]
        # iterate over all nodes which are connected to biomass reaction node
        for node in graph.getInOutNodes(n):
            # alias node and move it to a random position close to biomass reaction node (quick test hack)
            node_aliased = mvl.aliasNode(graph, node, biomass_re_id, x+random.random()*50, y+random.random()*50, None)
            # save newly created node in a list after aliasing
            biomass_node_list.append(node_aliased)


# move reaction node and all connected, aliased nodes 200 units to the right (x+200) and 200 units up (y+200)
re_node_coords = []
for node in biomass_node_list:
    x, y, z = viewLayout[node]
    viewLayout[node] = tlp.Coord(x+50, y+15, z) # shift coords here to move the biomass reaction independently
    if viewName[node] == "R70":
        # save coords of biomass reaction node
        x, y, z = viewLayout[node]
        re_node_coords.append(x)
        re_node_coords.append(y)
        re_node_coords.append(z)

# turn biomass node group into subgraph (but without reaction node, map it seperately!)
biomass_node_list = [n for n in biomass_node_list if viewName[n] != "R70"]
biomass_sg = mvl.createSubgraph(graph, "biomass reaction", biomass_node_list, None)

# increasing circle radius of layout "Circular" only possible by increasing size of nodes
viewSize = graph.getSizeProperty("viewSize")
for n in biomass_sg.getNodes():
    viewSize[n] = tlp.Size(3, 3, 0)

# apply layout to biomass subgraph
mvl.layoutSubgraph(biomass_sg, "Circular", None, re_node_coords[0], re_node_coords[1], 0)
"""

"""
# example 3: alias all exchange reactions (from ex_re_list, extracted from model with cbmpy if info available)
#            add them to a subgraph (new sg or the same as biomass) and layout them

# I.) construct first dictionary with key = ex reaction id and value = ex reaction object

# create a list of exchange reaction node objects (value list)
id_ = graph.getStringProperty("id")
ex_re_node_obj = [r for r in graph.getNodes() if id_[r] in ex_re_list]

# ex_re_list already contains the ids (key list)

# merge the key and value lists to a dictionary
ex_re_id_dict = {}
for i in range(len(ex_re_list)):
    ex_re_id_dict[ex_re_list[i]] = ex_re_node_obj[i]


# II.) construct second dictionary with key = ex reaction id and value = list of connected out nodes

# ex_re_list already contains the ids (key list)

# create a dict with ex reaction ids as keys and a list of connected in and out nodes as value
ex_re_con_dict = dict((id_[r], [n for n in graph.getInOutNodes(r) if id_[r] in ex_re_list]) for r in ex_re_node_obj)
print("\n dict of exchange reaction ids and lists of connected nodes:")
print(ex_re_con_dict)

# carlsson model hack: biomass reaction is classified as exchange reaction -> remove it from dict
del ex_re_con_dict['R70']

# III.) call aliasNode function for every node connected to an exchange reaction node
ex_re_alias_result_nodes = []
for r_id in ex_re_con_dict:
    for node in ex_re_con_dict[r_id]:
        n_new = mvl.aliasNode(graph, node, r_id, None, None)
        ex_re_alias_result_nodes.append(n_new)

print("\n list of nodes that were aliased (all exchange reaction nodes):")
print(ex_re_alias_result_nodes)

# IV.) add newly created, aliased nodes to a subgraph and layout it
ex_re_sg = mvl.createSubgraph(graph, "exchange reactions", ex_re_alias_result_nodes, None)
mvl.layoutSubgraph(ex_re_sg, "Tree Radial", None, 0, 0, 0)
"""

# The code below works and only needs a tulip graph
# take exchange reactions and biomass reaction from a graph and:
# 1.) get a subgraph with both and 2.) get individual subgraphs
"""
# example 4.1: basically the combination of examples 2 and 3
#              just create a subgraph with all the ex reaction nodes and the all biomass reaction nodes together
#              => test different layouts on that subgraph

# I.) exchange reaction nodes
# convert reaction node ids to reaction node objects
id_ = graph.getStringProperty("id")
ex_node_list = []
for r in graph.getNodes():
    # only consider reaction node whose id is in the exchage reaction list
    if id_[r] in ex_re_list:
        ex_node_list.append(r)
    # get species connected to ex reaction nodes
    for n in graph.getInOutNodes(r):
        if id_[r] in ex_re_list:
            ex_node_list.append(n)

# II.) all biomass reaction nodes
viewName = graph.getStringProperty("name")
biomass_node_list = []

# specific for carlsson model
biomass_name = "R70"

for n in graph.getNodes():
    # select biomass reaction node
    if viewName[n] == biomass_name:
        # add biomass reaction node to list
        biomass_node_list.append(n)
        # iterate over all nodes which are connected to biomass reaction node
        for node in graph.getInOutNodes(n):
            biomass_node_list.append(n)

# III.) create subgraph
biomass_and_ex_re_list = biomass_node_list + ex_node_list
biomass_and_ex_re_sg = mvl.createSubgraph(graph, "biomass and ex reactions", biomass_and_ex_re_list, None)
# => use this new subgraph to get to the exchange reactions instead of trying to extract them from the entire graph


# example 4.2 : alias all exchange reactions
#               add them as new subgraph to the subgraph biomass_and_ex_re_sg

id_ = biomass_and_ex_re_sg.getStringProperty("id")
clone_id = biomass_and_ex_re_sg.getStringProperty("clone_id")

# specific for carlsson model
biomass_id = "R70"

# I.)
# put all exchange reaction nodes in one list and all biomass related nodes in another
ex_re_nodes = []
biomass_nodes = []
for n in biomass_and_ex_re_sg.getNodes():
    if clone_id[n] != biomass_id and id_[n] != biomass_id:
        ex_re_nodes.append(n)
    else:
        biomass_nodes.append(n)

# II.) create seperate exchange reaction and biomass subgraphs
ex_re_only_sg = mvl.createSubgraph(biomass_and_ex_re_sg, "only exchange reactions", ex_re_nodes, None)
biomass_only_sg = mvl.createSubgraph(biomass_and_ex_re_sg, "only biomass reaction", biomass_nodes, None)

# III.) apply layout only on exchange reaction subgraph
mvl.layoutSubgraph(ex_re_only_sg, "Tree Radial", {"layer spacing": 10}, 0, 0, 0)
# Sugiyama: change node distance (default: 3) and layer distance (default: 3) to scale up the biomass reaction
mvl.layoutSubgraph(biomass_only_sg, "Sugiyama (OGDF)", {"node distance": 5, "layer distance": 5}, 0, 0, 0)

# TODO: fine tuning of the layout
# TODO: for tree raidal: straighten exchange reactions (in edge slope = out edge slope)

# VI.) save subgraph
tlp.saveGraph(biomass_and_ex_re_sg, "{}_biomass_ex_re.tlp".format(model_id))
"""


########################
# visualize FBA result #
########################
#"""
# color nodes and edges depending on flux value
# input format: dict = {"M1": -3, ...}

# specific for carlsson model
#flux_dict = {"R70": 2, "R93": 1, "R55": 0, "R42": -0.00001, "R7r": -1, "R10": -2}
flux_dict = fba_results

# save flux results as graph integer property "flux"
flux = graph.getDoubleProperty("flux")
id_ = graph.getStringProperty("id")

node_ids = [id_[r] for r in graph.getNodes()] # create a list of reaction node ids for entire network (keys)
node_objs = [r for r in graph.getNodes()] # create a list of reaction node objects for entire network (values)

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

# color total graph
mvl.setColours(graph, graph.getNodes(), COLOURS['grey'])
mvl.setColours(graph, graph.getEdges(), COLOURS['grey'])

# tulip provides a colour scale function!
# 1. initialize empty color scales
colorScale_pos = tlp.ColorScale([])
colorScale_neg = tlp.ColorScale([])

# 2. fill colorScales
colorScale_pos.setColorAtPos(0.0, tlp.Color(255, 255, 255))
colorScale_pos.setColorAtPos(1.0, COLOURS['green'])
colorScale_neg.setColorAtPos(0.0, tlp.Color(255, 255, 255))
colorScale_neg.setColorAtPos(1.0, COLOURS['red'])

# 3. get colors for flux values
# color in and out edges of flux nodes and reaction nodes
# TODO: use numpy to generate ranges for flux ranges numpy.linspace(start, end, number of intervalls)
viewColor = graph.getColorProperty("viewColor")
for n in network_nodes.values():
    if id_[n] in pos_flux:
        if flux[n] > 0 and flux[n] <= 0.5:
            viewColor[n] = colorScale_pos.getColorAtPos(0.25)
            for e in graph.getInOutEdges(n):
                viewColor[e] = colorScale_pos.getColorAtPos(0.25)

        if flux[n] > 0.5 and flux[n] <= 1:
            viewColor[n] = colorScale_pos.getColorAtPos(0.5)
            for e in graph.getInOutEdges(n):
                viewColor[e] = colorScale_pos.getColorAtPos(0.5)

        if flux[n] > 1 and flux[n] <= 3:
            viewColor[n] = colorScale_pos.getColorAtPos(0.75)
            for e in graph.getInOutEdges(n):
                viewColor[e] = colorScale_pos.getColorAtPos(0.75)

        if flux[n] > 3 :
            viewColor[n] = colorScale_pos.getColorAtPos(1)
            for e in graph.getInOutEdges(n):
                viewColor[e] = colorScale_pos.getColorAtPos(1)

    if id_[n] in neg_flux:
        if flux[n] < 0 and flux[n] <= -0.5:
            viewColor[n] = colorScale_neg.getColorAtPos(0.25)
            for e in graph.getInOutEdges(n):
                viewColor[e] = colorScale_neg.getColorAtPos(0.25)

        if flux[n] < -0.5 and flux[n] <= -1:
            viewColor[n] = colorScale_neg.getColorAtPos(0.5)
            for e in graph.getInOutEdges(n):
                viewColor[e] = colorScale_neg.getColorAtPos(0.5)

        if flux[n] < -1 and flux[n] <= -3:
            viewColor[n] = colorScale_neg.getColorAtPos(0.75)
            for e in graph.getInOutEdges(n):
                viewColor[e] = colorScale_neg.getColorAtPos(0.75)

        if flux[n] < -3 :
            viewColor[n] = colorScale_neg.getColorAtPos(1)
            for e in graph.getInOutEdges(n):
                viewColor[e] = colorScale_neg.getColorAtPos(1)

#"""


############################################################################
# save final graph, generate SBGN and SBML+Layout Files (Mimoza functions) #
############################################################################

# save total graph with all changes
tlp.saveGraph(graph, "{}_graph.tlp".format(model_id))

# export tulip graph as SVG file
#params = tlp.getDefaultPluginParameters('SVG Export', graph)
#tlp.exportGraph('SVG Export', graph, '{}.svg'.format(model_id), params
# NOT WORKING! -> so far: Exception: No plugin named  SVG Export.
# emailed tulip staff about it

# add layout to SBML file
groups_document = reader.readSBML(groups_sbml)
groups_model = groups_document.getModel()
gen_model = False
sbml_manager.save_as_layout_sbml(groups_model, gen_model, groups_sbml, n2lo, ub_sps)

# create SBGN file
groups_sbgn = os.path.join(cDir, '{}.sbgn'.format(model_id))
sbgn_helper.save_as_sbgn(n2lo, e2lo, groups_model, groups_sbgn)
