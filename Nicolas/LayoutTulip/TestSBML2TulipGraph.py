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
glyc_sps = ["GLU_6_P", "FRU_6_P", "FRU_BIS_P", "DHAP", "GA_3P", "PG", "PEP"]
glyc_re = ["R1", "R2r", "R3", "R4", "R5r", "R6r", "R7r", "R8r", "R9", "RR9"]
tca_sps = ["CITRATE", "ISOCIT", "AKG", "SUCC_CoA", "SUCC", "FUMARATE", "MALATE", "OXALO"]
tca_re = ["R21", "R23r", "R24", "R25", "R26r", "R27r", "R28r", "R29r", "R22"]
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



##################################################################
# different layouts on subsystems (compartments, pathways, etc.) #
##################################################################


# 1.) manually define subsystem node lists

# specific to toy_model_biofuel
small_circle = ["M1", "M2", "M9", "R02", "R10", "R11"]
larger_circle = ["M3", "M4", "M5", "M6", "M7", "R04", "R05", "R06", "R07", "R08"]

if sbml == os.path.join(cDir, "TestModels", "toy_model_biofuel.xml"):
    # create toymodel subgraphs and apply centered layout
    small_circle_subgraph = mvl.createSubgraph(graph, "small_circle", small_circle)
    larger_circle_subgraph = mvl.createSubgraph(graph, "larger_circle", larger_circle)
    mvl.layoutSubgraph(small_circle_subgraph, "Circular (OGDF)", 100, 100, 0)
    mvl.layoutSubgraph(larger_circle_subgraph, "Circular (OGDF)", -50, 25, 0)
    # save created toymodel subgraphs
    tlp.saveGraph(small_circle_subgraph, "{}_SmallCircle.tlp".format(model_id))
    tlp.saveGraph(larger_circle_subgraph, "{}_LargerCircle.tlp".format(model_id))


# 2.) get subsystem node list from all species of a certain compartment
#comp_ids = cmod.getCompartmentIds() # get a list of all compartments
#comp_2 = cmod.getCompartment(comp_ids[1]) # get  compartment 2 (second entry in list of all compartments)
#comp_2_species = comp_2.containsSpecies() # list of species in compartment 2

# layout all species from compartment 2 (e.g.: extracellular) seperately (i.e.: along the edges of the graph)
#comp_2_sg = mvl.createSubgraph(graph, "comp_2", comp_2_species)
#mvl.layoutSubgraph(comp_2_sg, "Hierarchical Graph", 0, 0, 0)


# 3.) get subsystem node list from different biological pathways (e.g.: glycolysis, tca and 'not specified')
glyc_sg = mvl.createSubgraph(graph, "glycolysis", glyc_sps)
tca_sg = mvl.createSubgraph(graph, "tca_cycle", tca_sps)
not_spec_sg = mvl.createSubgraph(graph, "not specified", not_spec_sps)

mvl.setColours(graph, glyc_sg.getNodes(), COLOURS['red'])
mvl.setColours(graph, tca_sg.getNodes(), COLOURS['blue'])

#mvl.layoutSubgraph(glyc_sg, "Mixed Model", -45, 90, 0)
#mvl.layoutSubgraph(tca_sg, "Circular", -45, -30, 0)
#mvl.layoutSubgraph(not_spec_sg, "Planarization Layout (OGDF)", None, None, 0)

# fix position of a specific node (M1 (id 2) for small circle, M3 (id 8) for larger circle)
# specific to toy_model_biofuel
#viewLayout[tlp.node(2)] = tlp.Coord(50,75,0)
#viewLayout[tlp.node(8)] = tlp.Coord(-50,25,0)
#viewLayout = graph.getLayoutProperty("viewLayout")
#viewLayout[tlp.node(16)] = tlp.Coord(100, 100, 0)


####################################################
# aliasing examples (using the function aliasNode) #
####################################################


# example 1: remove edge (edge 49) connecting biomass to TCA cycle at AKG (node 4)
#            get the edge that connects AKG (node 4) to biomass (R70, node 60)

# specific to carlsson_model
#aliasNode(graph, tlp.node(4), "R70", 22, 27)


# example 2: cut out the biomass function (individually alias all of its connected species)
#            and move it together with all of its now aliased species somewhere else
#            biomass function = R70 (node 60)

# alias the biomass reaction
# specific to carlsson_model (R70 is biomass reaction node)
viewName = graph.getStringProperty("name")
viewLayout = graph.getLayoutProperty("viewLayout")

# specific to carlsson model
biomass_re_id = "R70"
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


# move it and all connected. aliased nodes 200 units to the right (x+200) and 200 units up (y+200)
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
biomass_sg = mvl.createSubgraph(graph, "biomass reaction", biomass_node_list)

# increasing circle radius of layout "Circular" only possible by increasing size of nodes
viewSize = graph.getSizeProperty("viewSize")
for n in biomass_sg.getNodes():
    viewSize[n] = tlp.Size(3, 3, 0)

# apply layout to biomass subgraph
mvl.layoutSubgraph(biomass_sg, "Circular", re_node_coords[0], re_node_coords[1], 0)


# example 3: alias all exchange reactions (from ex_re_list, extracted from model with cbmpy if info available)
#            add them to a subgraph (new sg or the same as biomass) and layout them

# exchange reaction node aliasing attempts
"""
########################################################################################################################
# like example 2 but instead for only one reaction -> for a list of reactions => add one more loop
id_ = graph.getStringProperty("id")
viewName = graph.getStringProperty("name")
viewLayout = graph.getLayoutProperty("viewLayout")

# convert reaction node ids to reaction node objects
ex_node_list = []
for r in graph.getNodes():
    # only consider reaction node whose id is in the exchage reaction list
    if id_[r] in ex_re_list:
        ex_node_list.append(r)

ex_result_list = []
for re_node in ex_node_list:
    # iterate over all nodes which are connected to reaction node
    for node in graph.getInOutNodes(re_node):
        x, y, z = viewLayout[node]
        # alias node and move it to a random position close to reaction node
        node_aliased = mvl.aliasNode(graph, node, id_[re_node], x+random.random()*50, y+random.random()*50, None)
        # save newly created node in a list after aliasing
        ex_result_list.append(id_[node_aliased])

print(ex_result_list)
########################################################################################################################
id_ = graph.getStringProperty("id")

# convert reaction node ids to reaction node objects
ex_node_list = []
for r in graph.getNodes():
    # only consider reaction node whose id is in the exchage reaction list
    if id_[r] in ex_re_list:
        ex_node_list.append(r)

ex_result_list = []
for r in ex_node_list:
    for n in graph.getInOutNodes(r):
        # nodes of exchange reactions that have either only in or only out nodes already sit on the edges
        # -> no aliasing required, only add them to result list for subgraph creation
        if len(list(graph.getOutEdges(n))) == 0 or len(list(graph.getInEdges(n))) == 0:
            ex_node_list.append(n)
            # nodes of exchange reactions that have in and out nodes connect the exchange reactions to the rest of the
            # pathway => alias them
        elif len(list(graph.getOutEdges(n))) != 0 and len(list(graph.getInEdges(n))) != 0:
            # n = node to be aliased
            # the target strings come from a list of all nodes connected to the alias node (con_node)
            # except those nodes which are in ex_re_list (to not cut the exchange reactions themselves)
            for con_node in list(graph.getNodes()):
                if id_[con_node] not in ex_re_list:
                    n_new = mvl.aliasNode(graph, n, id_[con_node], 0, 0)
                    ex_node_list.append(n_new)
            else:
                pass
########################################################################################################################
"""


# example 4: basically the combination of examples 2 and 3
# just create a subgraph with all the ex reaction nodes and the all biomass reaction nodes
# => test different layouts on that subgraph


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
for n in graph.getNodes():
    # select biomass reaction node
    if viewName[n] == "R70":
        # add biomass reaction node to list
        biomass_node_list.append(n)
        # iterate over all nodes which are connected to biomass reaction node
        for node in graph.getInOutNodes(n):
            biomass_node_list.append(n)

# create subgraph
biomass_and_ex_re_list = biomass_node_list + ex_node_list
biomass_and_ex_re_sg = mvl.createSubgraph(graph, "biomass and ex reactions", biomass_and_ex_re_list)

# save subgraph seperately
tlp.saveGraph(biomass_and_ex_re_sg, "{}_biomass_ex_re.tlp".format(model_id))

# use this new subgraph to get to the exchange reactions instead of trying to extract them from the entire graph (TODO)

# example 5 (new version of example 3): alias all exchange reactions
#                                       add them to a subgraph (new sg or the same as biomass) and layout them



########################
# visualize FBA result #
########################
"""
# color nodes and edges depending on flux value (pos, neg, 0) -> do subgraphs for each (three colors)
# input format: dict = {"M1": -3, ...}

# specific for carlsson model
#flux_dict = {"R70": 2, "R93": 1, "R55": 0, "R42": -0.00001, "R7r": -1, "R10": -2}
flux_dict = fba_results

pos_flux = []
neg_flux = []
zero_flux = []

# extract reaction nodes from dict
for key in flux_dict:
    if flux_dict[key] > 0:
        pos_flux.append(key)
    elif flux_dict[key] < 0:
        neg_flux.append(key)
    else:
        zero_flux.append(key)

# create subgraphs for pos, neg and zero flux
pos_sg = createSubgraph(graph, "pos_flux", pos_flux)
neg_sg = createSubgraph(graph, "neg_flux", neg_flux)
zero_sg = createSubgraph(graph, "zero_flux", zero_flux)

# color total graph
setColours(graph, graph.getNodes(), COLOURS['grey'])
setColours(graph, graph.getEdges(), COLOURS['grey'])

# color in and out edges of flux nodes
viewColor = pos_sg.getColorProperty("viewColor")
for n in pos_sg.getNodes():
    for e in graph.getInOutEdges(n):
        viewColor[e] = COLOURS['green']

for n in neg_sg.getNodes():
    for e in graph.getInOutEdges(n):
        viewColor[e] = COLOURS['red']

#for n in zero_sg.getNodes():
#    for e in graph.getInOutEdges(n):
#        viewColor[e] = COLOURS['blue']

# color the flux subgraph nodes
setColours(graph, pos_sg.getNodes(), COLOURS['green'])
setColours(graph, neg_sg.getNodes(), COLOURS['red'])
#setColours(graph, zero_sg.getNodes(), COLOURS['blue'])

# reduce edge crossing
# in progress...
"""


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
