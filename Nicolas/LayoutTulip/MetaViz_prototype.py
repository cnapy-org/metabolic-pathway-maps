########################################################################################################################
#  MetaViz Prototype - Visualization Pipeline for Metabolic Networks                                                   #
#  Nicolas Huber, Dr. Brett Olivier, Prof. Bas Teusink - VU Amsterdam, 2016 - 2017                                     #
#  Generate a structured network map using only a SBML file (idealy with defined groups) as input. Resulting graph is  #
#  initially created as Tulip graph object and can be exported as SVG, JSON, GML, new SBML with layout plugin as well  #
#  as SBGN. Graph properties are also stored in a SQLite database.                                                     #
#  Necessary  dependencies for this pipeline (on Wind64 Anaconda):                                                     #
#  pip install cbmpy geojson tarjan libsbgnpy tulip-python                                                             #
#  SQLite module: sqlite 3 (https://www.tutorialspoint.com/sqlite/sqlite_installation.htm)                             #
########################################################################################################################


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
# MetaViz functions
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

# export generated graphs as .gml file
EXPORT_GML = True
# export generated graphs as .svg file
EXPORT_SVG = True
# export generated graphs as .json file
EXPORT_JSON = True

# colour the reaction nodes + in and outgoing edges based on flux values obtained from FBA
COLOR_FBA = False
COLOR_FBA_COLORBLIND = True



#################################################################
# read and initialize SBML file with libsbml and cbmpy, do FBA  #
#################################################################

## bgoli 180117: organising our tests
## bgoli 180117: other tests
#sbml = os.path.join(cDir, "TestModels", "others", "toy_model_biofuel.xml")
#sbml = os.path.join(cDir, "TestModels", "others", "e_coli_core.xml")
#sbml = os.path.join(cDir, "TestModels", "others", "Ecoli_Carlson2003.xml")
#sbml = os.path.join(cDir, "TestModels", "others", "Hpylori_iIT341.xml")
#sbml = os.path.join(cDir, "TestModels", "others", "Lactococcus_lactis_MG1363.xml")

## BiGG2 models (with pathway annotation)
#sbml = os.path.join(cDir, "TestModels", "bigg2", "iIT341.xml")
#sbml = os.path.join(cDir, "TestModels", "bigg2", "iJO1366.xml") # newest E. coli
#sbml = os.path.join(cDir, "TestModels", "bigg2", "iAF1260.xml") # E. coli about 2000+ reactions
sbml = os.path.join(cDir, "TestModels", "bigg2", "iJR904.xml") # E. coli about 1000 reactions
#sbml = os.path.join(cDir, "TestModels", "bigg2", "iND750.xml") # Yeast about 1000 reactions

## Willis models
#sbml = os.path.join(cDir, "TestModels", "willi", "M5_6.xml")
#sbml = os.path.join(cDir, "TestModels", "willi", "M5_6_bigggroups.xml")
#sbml = os.path.join(cDir, "TestModels", "willi", "M5_6_compgroups.xml")


# read SBML file and create libsbml model
reader = libsbml.SBMLReader()
doc = reader.readSBML(sbml)
model = doc.getModel()
model_id = model.getId()

# create CBMPy model
cmod = cbmpy.readSBML3FBC(sbml)
# change flux bounds values to not mess with the coloring later on
cmod.changeAllFluxBoundsWithValue(999999, cbmpy.INF)
cmod.changeAllFluxBoundsWithValue(-999999, cbmpy.NINF)
cmod.changeAllFluxBoundsWithValue(1000, cbmpy.INF)
cmod.changeAllFluxBoundsWithValue(-1000, cbmpy.NINF)
# remove non-reacting metabolites (since they would be visualized as disconnected nodes)
if len(cmod.deleteNonReactingSpecies()) != 0:
    cmod.deleteNonReactingSpecies(simulate=False)
# get stoichiometric matrix
stM = cmod.N
# get objective reaction id
objective_re_id = str(cmod.getActiveObjectiveReactionIds()[0])

# analyze model (FBA) and save results in a flux dictionary
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
# IF GROUPS DON'T EXIST -> provide pathway lists manually (example: e coli carlson model):
if PATHWAY_LISTS and sbml == os.path.join(cDir, "TestModels", "others", "Ecoli_Carlson2003.xml") and AUTO_ALIAS is \
        False:
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


# IF GROUPS EXIST -> fetch groups from sbml (this doesn't use ub_sps_list or ub_re_list)
# create dictionary of groups where each group name is a key and the values are again dicts with three keys: 'id',
# 'reactions' and 'metabolites'
path_dict = {}
if PATHWAY_LISTS and AUTO_ALIAS is False:

    # manually create id lists for a few selected pathways of model iIT341
    #glyc_sps, glyc_re = mvl.createPathwayLists(cmod, 'Glycolysis')
    #tca_sps, tca_re = mvl.createPathwayLists(cmod, 'TCA cycle')
    #fas_sps, fas_re = mvl.createPathwayLists(cmod, 'Fatty Acid Synthesis')
    #uc_sps, uc_re = mvl.createPathwayLists(cmod, 'Urea Cycle')
    #am_sps, am_re = mvl.createPathwayLists(cmod, 'Aminosugar Metabolism')

    # list of all reaction objects that are associated with a group
    all_groups_member_list = []
    for group in cmod.groups:
        for member in group.members:
            all_groups_member_list.append(member)
    # group all the ungrouped reactions and their metabolites in a new group "undefined" to have
    # them later (otherwise they won't show up in the output graphs since they are not in path_dict)
    cmod.createGroup('undefined')
    undef = cmod.getGroup('undefined')
    undef.setName('Undefined')
    for reaction in cmod.reactions:
        if reaction not in all_groups_member_list:
            undef.addMember(reaction)

    # loop over all groups in iIT341 and automatically create id lists for all groups
    # save all these lists in a new dictionary
    for group in cmod.groups:
        path_sps, path_re = mvl.createPathwayLists(cmod, group.name, small_threshold=30)
        path_dict.update({group.name: {"metabolites": path_sps, "reactions": path_re}})

    print("Pathway dictionary generated.")



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
# create initial nodes for all species
id2n = mvl.speciesToNodes(graph, cmod, ub_list, graph_props)

# TODO: initial tulip graph creation is slow because of Mimoza Black Box Code!
# TODO: reaction node creation still relies on copied code (problem: reaction node compartment assignment needs helper
# TODO: functions -> understand what the copied code does and rewrite!)
#START COPIED CODE######################################################################################################
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
#END COPIED CODE########################################################################################################

# clone/alias nodes selected above (in ub_list)
# -> chebi info was apparently only important for automatically extracting ubiquitous species ids! (we now provide them)
# => provide them by passing a list (ubs), converted to a set (cloneNodes needs unique ids)
root = sbml2tlp.cloneNodes(set(ub_list), graph)



##############################
#run Mimoza layout functions #
##############################

# actual layouting, creation of graph structure
# n2lo and e2lo are dictionaries containing node and edge layout info (id, color, coords, etc.) respectively
fc, (n2lo, e2lo), hidden_c_ids, c_id_hidden_ubs, graph = tulip_graph2geojson.graph2geojson(c_id2info, c_id2outs, root)

# save raw Mimoza graph output as backup
tlp.saveGraph(graph, "{}_raw.tlp".format(model_id))
print("Initial Tulip graph created.")



###############################################################
# map sbml element ids to initial graph nodes in a dictionary #
###############################################################

# Section Description:
# --------------------
# in the initial graph every model id is represented by exactly one node (if AUTO_ALIAS is off!)
# example: iIT341 initial graph (returned by Mimoza): 1039 nodes created
#          check: len(cmod.getSpeciesIds()+cmod.getReactionIds()) returns 1039 as well
#
# => create mapping between ids and their nodes in a dict that reflects the state of the raw graph
#    dict structure: id (key): value is another dict with layout (key): True/False and initial_node (key): node object
#    layout key serves as a flag to see whether an initial node has already been used for layouting later (True = has
#    been used)
#    -> if the downstream layouting function sees layout True => creates new node and copies properties
#                                                                (new node will not be added to this dict)

# use only if one model id = one graph node is true! (only if auto aliasing is off and  the current graph oject is
# in its raw, initial state)
graph_elements = {}
if AUTO_ALIAS == False:
    # all species and all reactions of sbml model = total set of model elements
    model_elements = cmod.species + cmod.reactions
    #dict with keys: model ids and values = dict with keys: pathway names and values = associated graph nodes
    for m_e in model_elements:
        graph_elements[m_e.getId()] = {}
        # use [0] to get first entry of that list so that the returned object is the nodes itself and not a list again
        graph_elements[m_e.getId()]["initial_node"] = [n for n in graph.getNodes() if graph_props["id"][n] ==
                                                       m_e.getId()][0]
        # edges do not exist in SBML model! (and created graph edges have no ids, just names)
        # -> edges are not included in dict!
        # track if graph element was used for layouting already (default value: False, gets changed to True once the
        # respective node was layouted in a subgraph)
        graph_elements[m_e.getId()]["layout"] = False

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

# Why is a mapping between ids and graph elements (nodes, edges) necessary?
# -------------------------------------------------------------------------
#
# There is a difference between a model species id and a graph node id!
# => It is not sufficient to select nodes based on model ids (at some point multiple nodes with the same id can exist)
# => handle node objects!
#
# Problem: if a node is shared between pathways, the last valid layouting call will determine its position meaning that
# previous positions of that node in other layouted subgraphs that also contain it are overwritten!
#
# example in Carlsson model: since the initial nodes with the ids "PEP" and "R1" are part of glycolysis AND exchange
# reactions, the assigned positions of the "PEP" and "R1" initial nodes after layouting glyc are overwritten by
# layouting exch re!
#
# -> give subgraphs their own versions of a shared node by tracking what graph elements have already been used for
#    layouting to determine what nodes and edges need to be recreated



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

# TODO: add labels to edges (possible? alternative?)
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
           'light grey': tlp.Color(235, 235, 235),
           'black': tlp.Color(0, 0, 0)
           }

# set colour of all edges (complete graph)
mvl.setColours(graph, graph.getEdges(), COLOURS['grey'])

# change line types (Polyline, BezierCurve, CatmullRomCurve, CubicSplineCurve)
#for e in graph.getEdges():
#    graph_props["viewShape"][e] = tlp.EdgeShape.BezierCurve



##############################################################
# layout pathways in subgraphs (Glycolysis, TCA Cycle, etc.) #
##############################################################

# initialize a dict of all existing graphs (key: name, value: graph object)
#  -> everytime a new, important graph is created it gets added to the graph registry
graph_registry = {}

# needs dict graph elements which is only created it AUTO_ALIAS is false since it needs 1:1 model id to graph elements
if LAYOUT_PATHWAYS is True and AUTO_ALIAS is False:

    # layout all pathways in graph (entire path_dict) to generate subgraphs
    for path_name, path_ids in path_dict.items():
        # fetch pathway ids
        path_list = list(path_ids["metabolites"]) + list(path_ids["reactions"])
        # choose appropriate layout algorithm
        path_layout = []
        if "ycle" in path_name:
            path_layout.append("Circular (OGDF)")
        else:
            path_layout.append("Sugiyama (OGDF)")
        # apply layout to pathway and track already layouted graph elements
        path_sg = mvl.layoutPathway(graph, graph_elements, graph_props, path_list, path_name, path_layout[0])
        # set name of subgraph
        path_sg.setName(path_name)

        # test if all reaction and metabolite ids of the group definition are part of the generated subgraph
        # unique set of ids of generated graph nodes
        generated_nodes = set([graph_props["id"][n] for n in path_sg.getNodes()])
        # unique set of input ids
        input_ids = set(path_list)
        # both need to have the same length since at least one node needs to exist for every unique input id
        # TODO. after having tested identity
        # TODO: -> do an intersection with the two sets (to catch duplications, etc.)
        if len(generated_nodes) == len(input_ids):
            print("{} subgraph complete.".format(path_name))

    # create and layout the objective reaction block (sg gets added to 'graph')
    obj_re_sg = mvl.layoutObjFunc(graph, objective_re_id, "Sugiyama (OGDF)", graph_props)
    obj_re_sg.setName("objective reaction")
    # after layouting add graph to graph registry
    graph_registry.update({"total graph": graph})


    # 1.) default dump (all pathways in one SVG file)
    # -> combine all pathways to create a new output graph of the whole network
    whole_output_graph = mvl.dumpToGraph(graph, graph_props, "output_whole", [sg.getName() for sg in graph.getSubGraphs()])
    graph_registry.update({whole_output_graph.getName(): whole_output_graph})

    # 2.) selected dump (selected pathways in one SVG file)
    # create a new output graph with a few selected pathways combined
    selected_output_graph = mvl.dumpToGraph(graph, graph_props, "output_selected", ["objective reaction", "Transport"])
    graph_registry.update({selected_output_graph.getName(): selected_output_graph})

    # 3.) sequential dump (all pathways in individual SVG files)
    individual_output_graphs = []
    for sg in graph.getSubGraphs():
        # append the returned output graph for each input subgraph to list
        sg_output_graph = mvl.dumpToGraph(graph, graph_props, sg.getName(), [sg.getName()])
        individual_output_graphs.append(sg_output_graph)
        graph_registry.update({sg_output_graph.getName(): sg_output_graph})

    # TODO: control positions of pathway subgraphs
    # TODO: get highest and lowest x and y coordinate values for each subgraph (max(list(x values)) and max(list(y
    # TODO: values)))
    # TODO: then. build a numpy.array() (a matrix) that serves as a template for canvas positions - with subgraph ids
    # TODO: as matrix entries:
    # TODO: example: if matrix entry at (1,1) is 1 -> put "Glycolysis" subgraph (id: 1) in top left corner etc.
    # TODO: moving subgraphs is possible with center Layout



########################
# visualize FBA result #
########################

# color nodes and edges depending on flux value
# input format: dict = {"M1": -3, ...}
if COLOR_FBA:
    for graph in graph_registry.values():
        mvl.colourFBA(graph, fba_results, COLOURS["green"], COLOURS["red"], COLOURS["light grey"])
if COLOR_FBA_COLORBLIND:
    for graph in graph_registry.values():
        mvl.colourFBA(graph, fba_results, tlp.Color.Cyan, tlp.Color.Magenta, COLOURS["light grey"])

# TODO: visualize existence of reactions of a smaller model on a big map (e coli vs mutated e coli eg.)



############################################################################
# save final graph, generate SBGN and SBML+Layout Files (Mimoza functions) #
############################################################################

# save total graph with all changes ('messy' layout)
tlp.saveGraph(graph_registry["total graph"], "{}_graph.tlp".format(model_id))
print("Base graph file: {}_graph.tlp".format(model_id))

if AUTO_ALIAS is False:
    for graph_name, graph_object in graph_registry.items():
        if graph_name == "output_whole":
            # save graph with whole network layout
            tlp.saveGraph(graph_object, "{0}_output_whole.tlp".format(model_id))
            print("Whole network visualization: {0}_output_whole.tlp".format(model_id))
        if graph_name == "output_selected":
            # save graph 'F' with selected pathway layouts
            tlp.saveGraph(graph_object, "{0}_output_selected.tlp".format(model_id))
            print("Visualization of selected pathways: {0}_output_combined.tlp".format(model_id))

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



#######################################
# export the graph (JSON, GML or SVG) #
#######################################

# create output folder
output_path = os.path.join(cDir, '{}_Output'.format(model_id))
try:
    os.makedirs(output_path)
except OSError:
    if not os.path.isdir(output_path):
        raise

# export every graph in graph registry
if EXPORT_GML:
    for graph in graph_registry.values():
        gml_path = os.path.join(output_path, "{0}_{1}.gml".format(model_id, graph.getName()))
        # replace forward dashes in graph name to not cause problems with path
        gml_path = gml_path.replace('/', '-')
        mvl.exportGML(graph, model_id, "gml_output_graph.tlp", gml_path)
        print("Graph exported as: {0}_{1}.gml".format(model_id, graph.getName()))

if EXPORT_SVG:
    # get default parameters of SVG export plugin
    params = tlp.getDefaultPluginParameters('SVG Export', graph)
    for graph in graph_registry.values():
        svg_path = os.path.join(output_path, "{0}_{1}.svg".format(model_id, graph.getName()))
        # replace forward dashes in graph name to not cause problems with path
        svg_path = svg_path.replace('/', '-')
        tlp.exportGraph('SVG Export', graph, svg_path, params)
        print("Graph exported as: {0}_{1}.svg".format(model_id, graph.getName()))

if EXPORT_JSON:
    # get default parameters of JSON export plugin
    params = tlp.getDefaultPluginParameters('JSON Export', graph)
    for graph in graph_registry.values():
        json_path = os.path.join(output_path, "{0}_{1}.json".format(model_id, graph.getName()))
        # replace forward dashes in graph name to not cause problems with path
        json_path = json_path.replace('/', '-')
        tlp.exportGraph('JSON Export', graph, json_path, params)
        print("Graph exported as: {0}_{1}.json".format(model_id, graph.getName()))



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

# select the graph object (e.g. total graph)
graph = graph_registry["total graph"]

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
    if prop_type == 'double': # in python: float (in SQL -> REAL)
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
    if prop_type == 'double': # in python: float (in SQL -> REAL)
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
# create empty connectivity table (each row is an edge with id, source node id and target node id as columns)
c.execute('''CREATE TABLE IF NOT EXISTS edge_connectivity(edge_id INTEGER PRIMARY KEY,
source_node_id TEXT,
target_node_id TEXT)
''')
# create node connectivity table aswell: for each node a list (= string) of all connected In and Out edges
c.execute('''CREATE TABLE IF NOT EXISTS node_connectivity(node_id INTEGER PRIMARY KEY,
in_out_edge_list TEXT)
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
# table graph_node
for l in graph_props_node_insertlist:
    c.execute('''INSERT OR IGNORE INTO graph_nodes VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,
    ?,?,?,?,?,?,?,?,?,?,?,?,?)''', l)
# table graph_edges
for l in graph_props_edge_insertlist:
    c.execute('''INSERT OR IGNORE INTO graph_edges VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,
    ?,?,?,?,?,?,?,?,?,?,?,?,?)''', l)
# table connectivity
# create dict with key = str(node) and value = integer id of the node
node_str_int_dict = {}
for node, i in zip(graph.getNodes(), range(len(list(graph.getNodes())))):
    node_str_int_dict[str(node)] = i
for edge, i in zip(graph.getEdges(), range(len(list(graph.getEdges())))):
    insert_row = [i, node_str_int_dict[str(graph.source(edge))], node_str_int_dict[str(graph.target(edge))]]
    c.execute('''INSERT OR IGNORE INTO edge_connectivity VALUES (?,?,?)''', insert_row)
for node, i in zip(graph.getNodes(), range(len(list(graph.getNodes())))):
    row = [i, str(list(graph.getInOutEdges(node)))]
    c.execute('''INSERT OR IGNORE INTO node_connectivity VALUES (?,?)''', row)
print("Graph data has been stored in SQL database.")



#####################################
# retrieve graphs from SQL database #
#####################################


# 4 tables in database:
# - graph_nodes: every row is a node, each column is a node property
# - graph_edges every row is an edge, each column is an edge property
# - node_connectivity: every row is a node, first column contains node ids and second column contains lists of
#                      connected edges
# - edge connectivity: every row is an edge, first column contains edge ids, second column contains source node ids
#                      and third column contains target node ids

# retrieve entire table 'graph_nodes'
c.execute('''SELECT * FROM graph_nodes''')
graph_nodes_rows = c.fetchall() # list of row tuples
# change list of row tuples to list of row lists
graph_nodes_rows = [list(tup) for tup in graph_nodes_rows]

# sqlite stores strings in unicode!
# -> 'test' is automatically converted to u'test'
# convert unicode strings to ascii strings
for ind, row_list in enumerate(graph_nodes_rows):
    graph_nodes_rows[ind] = [str(item) if isinstance(item, unicode) else item for item in row_list]
# decode queried data (string -> list)
for row_list in list(graph_nodes_rows):
    for ind, item in enumerate(row_list):
        # go over every item in a row list. for every string item:
        # if the string was generated from a tup of multiple numbers (x,x,x,x)
        # -> convert ascii string back to list
        # -> identify the strings that originate from tuples by only taking strings without letters
        if isinstance(item, str):
            item_list = []
            if ',' in item and any(e.isalpha() for e in item) == False:
                # convert item string to a list
                item_list.append([int(string_part) if isinstance(string_part, int) else float(string_part)
                     for string_part in item.split(',')])
            # replace encoded item with decoded item (if it was converted, otherwise keep old entry)
            if item_list:
                row_list[ind] = item_list
            else:
                row_list[ind] = item

# TODO: use retrieved information to rebuild the whole graph and save it as new tlp (see if layout stay constant)
# TODO: select just a pathway from the database, extract nodes and edges and generate the graph
# TODO: use biomass layout for testing purposes

# commit changes to database
db.commit()
# close connection to database
db.close()











