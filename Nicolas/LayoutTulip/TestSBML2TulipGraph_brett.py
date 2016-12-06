# Test transition from sbml file to structured tulip graph, sbml with layout and sbgn file
# necessary  dependencies (Wind64 Anaconda): pip install geojson tarjan libsbgnpy tulip-python

import os
# get current directory
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))

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



############################
# global control variables #
############################

# select a certain set of nodes and use mimoza aliasing + layouting
AUTO_ALIAS = False
PATHWAY_ALIAS = True

# custom aliasing and layouting of only a certain part of the graph (only exchange reactions, only obj reaction)
EX_RE_ALIAS = True
ALIAS_OBJ = True

LAYOUT_GLYC_TCA = True # needs PATHWAY_ALIAS true, so that necessary node lists are generated (for pathways)
LAYOUT_EX_RE_AND_BIOMASS = True # needs EX_RE_ALIAS true, so that necessary exchange reaction node list is generated

COLOUR_FBA = False



#############################################
# read SBML with libsbml and cbmpy, do FBA  #
#############################################

# read SBML file (with libsbml)
#sbml = os.path.join(cDir, "TestModels", "toy_model_biofuel.xml") # 4, 4, False
#sbml = os.path.join(cDir, "TestModels", "e_coli_core.xml") # 6, 4, True
sbml = os.path.join(cDir, "TestModels", "Ecoli_Carlson2003.xml") # 6, 4, True
#sbml = os.path.join(cDir, "TestModels", "Hpylori_iIT341.xml") # 6, 4, True
#sbml = os.path.join(cDir, "TestModels", "iJO1366.xml") # 6, 4, True
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

# get biomass reaction id
biomass_re_id = str(cmod.getActiveObjectiveReactionIds()[0])

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

# list of species nodes and reaction nodes to be aliased (filled by one or more of the subsequent statements)
ub_sp_list = []
ub_re_list = []

# 1.) manually provide alias node list <-> hard coded by user (for testing purposes)
# specific to toy_model_biofuel
if sbml == os.path.join(cDir, "TestModels", "toy_model_biofuel.xml"):
    ub_sp_list.extend(['ATP', 'ADP', 'M8', 'M10', 'M11', 'M12', 'M13', 'biomass', 'S', 'fuel1', 'fuel2'])
    ub_re_list.extend(['R09', 'R12', 'R13', 'R14', 'R15', 'R16', 'R17', 'R18'])

# 2.) create alias node list based on biological pathways <-> group info in sbml - if available (strategy 2)
# specific to carlsson model
glyc_sps = []
tca_sps = []
glyc_re = []
tca_re = []
if PATHWAY_ALIAS and sbml == os.path.join(cDir, "TestModels", "Ecoli_Carlson2003.xml"):
    # create lists of species ids and reaction ids for each pathway (not necessary if provided by sbml)
    # lists of metabolite nodes
    glyc_sps.extend(["GLU_6_P", "FRU_6_P", "FRU_BIS_P", "DHAP", "GA_3P", "PG", "PEP"])
    tca_sps.extend(["CITRATE", "ISOCIT", "AKG", "SUCC_CoA", "SUCC", "FUMARATE", "MALATE", "OXALO"])
    # lists of reaction nodes
    glyc_re.extend(["R1", "R2r", "R3", "R4", "R5r", "R6r", "R7r", "R8r", "R9", "RR9"])
    tca_re.extend(["R23r", "R24", "R25", "R26r", "R27r", "R28r", "R29r", "R22"])
    # ubiquitous small metabolites
    small_sps = ["ATP", "ADP", "NAD", "NADH", "CO2", "NH3"]
    # create list of all reaction and species nodes that are not part of tca or glyc
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
    ub_sp_list.extend(tca_sps_bound)
    ub_re_list.extend([])

# adding hand coded Pylori from BiGG2
if PATHWAY_ALIAS and sbml == os.path.join(cDir, "TestModels", "iIT341.xml"):
    # create lists of species ids and reaction ids for each pathway (not necessary if provided by sbml)
    gids = {}
    for g in cmod.groups:
        if g.getName() == 'Glycolysis':
            gids['Glycolysis'] = g.getId()
        elif g.getName() == 'TCA cycle':
            gids['TCA cycle'] = g.getId()

    # lists of reaction nodes
    glyc_re = cmod.getGroup(gids['Glycolysis']).getMemberIDs()
    tca_re = cmod.getGroup(gids['TCA cycle']).getMemberIDs()

    # lists of metabolite nodes
    glyc_sps = []
    tca_sps = []
    for r in glyc_re:
        glyc_sps.extend(cmod.getReaction(r).getSpeciesIds())
    for r in tca_re:
        tca_sps.extend(cmod.getReaction(r).getSpeciesIds())
    glyc_sps = list(set(glyc_sps))
    tca_sps = list(set(tca_sps))
    # ubiquitous small metabolites
    small_sps = [s.getId() for s in cmod.species if len(s.isReagentOf()) > 30]

    for ubs in small_sps:
        if ubs in glyc_sps:
            glyc_sps.remove(ubs)
        if ubs in tca_sps:
            tca_sps.remove(ubs)

    # create list of all reaction and species nodes that are not part of tca or glyc
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
    ub_sp_list.extend(tca_sps_bound)
    ub_re_list.extend([])



# 3.) create alias node list based on connectivity <-> user defined thresholds (strategy 1)
if AUTO_ALIAS:
    sp_list, re_list = mvl.identifyConnectedNodes(cmod, 6, 4, True)
    ub_sp_list.extend(sp_list)
    ub_re_list.extend(re_list)

# 4.) create alias node list based on a certain compartment <-> compartment info in sbml (strategy 3)
    #comp_ids = cmod.getCompartmentIds() # get a list of all compartments
    #comp_2 = cmod.getCompartment(comp_ids[1]) # get  compartment 2 (first entry in list of all compartments)
    #comp_2_species = comp_2.containsSpecies() # list of species in compartment 2
    #comp_2_reac = comp_2.containsReactions() # list of reactions in compartment 2
    #ub_sp_list.extend(comp_2_species)
    #ub_re_list.extend(comp_2_reac)

# 5.) create alias node list of all boundary (exchange) reactions (for aliasNode() function, NOT Mimoza!)
ex_re_list = []
if EX_RE_ALIAS:
    ex_re_list.extend(cmod.getExchangeReactionIds())

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

# create and map graph properties for new tulip graph
anc_term, anc_id, anc_name, term, c_id, id_, name, type, rev, stoich, transport, ub, viewLabel, \
viewLayout, viewShape, viewSize, viewColor, fake, clone_id, formula, rel_c_ids, edge_direct, flux, selection = mvl.create_props(graph)

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
    viewShape[e] = tlp.EdgeShape.BezierCurve



########################################################################
# apply different layouts on subsystems (compartments, pathways, etc.) #
########################################################################

# 1.) manually define subsystem node lists and apply layout
# specific to toy_model_biofuel
if sbml == os.path.join(cDir, "TestModels", "toy_model_biofuel.xml"):
    # create lists with pathways ids
    small_circle = ["M1", "M2", "M9", "R02", "R10", "R11"]
    larger_circle = ["M3", "M4", "M5", "M6", "M7", "R04", "R05", "R06", "R07", "R08"]
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

# 3.) get subsystem node list from different biological pathways (e.g.: glycolysis, tca) and apply layout
if LAYOUT_GLYC_TCA and len(glyc_re) > 1 and len(tca_re) > 1:
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
    # for edge properties
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
    # since only nodes (reaction and metabolite) are moved:
    # #delete edges from subgraph and redraw them (add) after nodes have been moved!
    # record edge information in dict before nodes are moved (multiple edges have the same name!)
    gly_edge_dict = {}
    for e in glyc_sg.getEdges():
        gly_edge_dict[str(e)] = [name[e], glyc_sg.source(e), glyc_sg.target(e), ub[e], viewColor[e], viewSize[e],
                             viewShape[e]]
    # move glycolysis => edges will stay behind but we delete them after moving and redraw them
    buff = 60 # buffer space between glycolysis and tca
    tlp.LayoutProperty.center(viewLayout, tlp.Coord(x_trans, 150+y_diff+buff, 0), glyc_sg)
    # delete edges from glyc subgraph (edges will still exist in the main graph)
    for e in glyc_sg.getEdges():
        glyc_sg.delEdge(e)
    # redraw the edges (add them to glyc subgraph with info stored in gly_edge_dict)
    for e_name, e_infos in gly_edge_dict.items():
        e_new = glyc_sg.addEdge(e_infos[1], e_infos[2])
        name[e_new] = e_infos[0]
        ub[e_new] = e_infos[3]
        viewColor[e_new] = e_infos[4]
        viewSize[e_new] = e_infos[5]
        viewShape[e_new] = e_infos[6]
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
    # remove edge crossings between the pathways:
    #  tlp.computeLinesIntersection(line1, line2)
    # ...

# fix position of a specific node (M1 (id 2) for small circle, M3 (id 8) for larger circle)
# specific to toy_model_biofuel
#viewLayout[tlp.node(2)] = tlp.Coord(50,75,0)
#viewLayout[tlp.node(8)] = tlp.Coord(-50,25,0)
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
if ALIAS_OBJ:
    # alias the biomass reaction
    mvl.aliasObjFunc(graph, biomass_re_id, "Sugiyama (OGDF)", [name, viewLayout, viewSize])

# example 3: alias all exchange reactions (from ex_re_list, extracted from model with cbmpy if info available)
#            add them to a subgraph (new sg or the same as biomass) and layout them

if EX_RE_ALIAS:
    # I.) construct first dictionary with key = ex reaction id and value = ex reaction object
    # create a list of exchange reaction node objects (value list)
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
    if sbml == os.path.join(cDir, "TestModels", "Ecoli_Carlson2003.xml"):
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
    mvl.setColours(ex_re_sg, ex_re_sg.getEdges(), COLOURS["grey"])
    mvl.layoutSubgraph(ex_re_sg, "Tree Radial", None, 0, 0, 0)

# take exchange reactions and biomass reaction from a graph and:
# 1.) get a subgraph with both and 2.) get individual subgraphs
if LAYOUT_EX_RE_AND_BIOMASS:
    #
    # example 4.1: basically the combination of examples 2 and 3
    #              just create a subgraph with all the ex reaction nodes and the all biomass reaction nodes together
    #              => test different layouts on that subgraph
    #
    # I.) exchange reaction nodes
    # convert reaction node ids to reaction node objects
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
    biomass_node_list = []
    for n in graph.getNodes():
        # select biomass reaction node
        if name[n] == biomass_re_id:
            # add biomass reaction node to list
            biomass_node_list.append(n)
            # iterate over all nodes which are connected to biomass reaction node
            for node in graph.getInOutNodes(n):
                biomass_node_list.append(node)
    # III.) create subgraph
    biomass_and_ex_re_list = biomass_node_list + ex_node_list
    biomass_and_ex_re_sg = graph.inducedSubGraph(biomass_and_ex_re_list)
    biomass_and_ex_re_sg.setName("biomass and ex reactions")
    #biomass_and_ex_re_sg = mvl.createSubgraph(graph, "biomass and ex reactions", biomass_and_ex_re_list, None)
    # => use this new subgraph to get to the exchange reactions instead of trying to extract them from the entire graph
    #
    # example 4.2 : alias all exchange reactions
    #               add them as new subgraph to the subgraph biomass_and_ex_re_sg
    # I.)
    # create subgraphs for biomass and ex reactions
    # ex_re_sg already exists as subgraph for only exchange reactions
    biomass_only_sg = graph.inducedSubGraph(biomass_node_list)
    biomass_only_sg.setName("only biomass reaction")
    #biomass_only_sg = mvl.createSubgraph(biomass_and_ex_re_sg, "only biomass reaction", biomass_nodes, None)
    # III.) apply layout only on exchange reaction subgraph
    mvl.layoutSubgraph(ex_re_sg, "Tree Radial", {"layer spacing": 10}, 0, 0, 0)
    # Sugiyama: change node distance (default: 3) and layer distance (default: 3) to scale up the biomass reaction
    mvl.layoutSubgraph(biomass_only_sg, "Sugiyama (OGDF)", {"node distance": 5, "layer distance": 5}, 0, 0, 0)
    # TODO: fine tuning of the layout
    # TODO: for tree raidal: straighten exchange reactions (in edge slope = out edge slope)
    # VI.) save subgraph
    tlp.saveGraph(biomass_and_ex_re_sg, "{}_biomass_ex_re.tlp".format(model_id))



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
