# Test transition from sbml file to structured tulip graph, sbml with layout and sbgn file
# necessary  dependencies (Wind64 Anaconda): pip install geojson tarjan libsbgnpy tulip-python

import os
# get current directory
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))

import libsbml
import cbmpy
from tulip import tlp

# Core Mimoza modules
import sbml2tlp
import sbml_manager
import tulip_graph2geojson
import sbgn_helper
import graph_properties


#################################################################
# read SBML, prepare Node Aliasing and extract Compartment Info #
#################################################################

# read SBML file (with libsbml)
#sbml = os.path.join(cDir, "TestModels", "toy_model_biofuel.xml") # 4, 4, False
#sbml = os.path.join(cDir, "TestModels", "e_coli_core.xml") # 6, 4, True
sbml = os.path.join(cDir, "TestModels", "Ecoli_Carlson2003.xml") # 6, 4, True
#sbml = os.path.join(cDir, "TestModels", "Hpylori_iIT341.xml") # 6, 4, True
#sbml = os.path.join(cDir, "TestModels", "iJO1366.xml") # 6, 4, True

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
    fba_results[r.getId()] = round(r.getValue(),6)

# create libsbml model
reader = libsbml.SBMLReader()
doc = reader.readSBML(sbml)
model = doc.getModel()
model_id = model.getId()

def autoSplitConnectedNodes(model, spec_split, reac_split, split_biomass=True):
    """
    Automatically tries to alias highly connected metabolites and nodes::

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

ub_sp_list, ub_re_list = autoSplitConnectedNodes(cmod, 6, 4, True)
# combine species and reactions to one list of ubiquitous nodes
ub_list = ub_sp_list + ub_re_list

# additionally alias node list from a certain compartment
#comp_ids = cmod.getCompartmentIds() # get a list of all compartments
#comp_2 = cmod.getCompartment(comp_ids[1]) # get  compartment 2 (first entry in list of all compartments)
#comp_2_species = comp_2.containsSpecies() # list of species in compartment 2
#comp_2_reac = comp_2.containsReactions() # list of reactions in compartment 2
#ub_sp_list = comp_2_species
#ub_re_list = comp_2_reac

# manually provide alias node list
#ub_sp_list = ['ATP', 'ADP', 'M8', 'M10', 'M11', 'M12', 'M13', 'biomass', 'S', 'fuel1', 'fuel2']
#ub_re_list = ['R09', 'R12', 'R13', 'R14', 'R15', 'R16', 'R17', 'R18']

# add ubiquitous species group to sbml file and save it as "..._with_groups.xml"
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

# extract compartment hierarchy, save to dictionary
# create c_id2level dictionary for function import_sbml
# c_id2level needs to be automatically extracted from input SBML file
# required form: c_id2level = {'Cell': (0, ''), 'Boundary': (1, 'Cell')}
comp_list = [''] + cmod.getCompartmentIds()
level_list = [(x, comp_list[x]) for x in range(len(comp_list))]
c_id2level = dict(zip(comp_list[1:], level_list))


#######################################################
# create Tulip Graph and run Mimoza layout functions  #
#######################################################

# create new tulip graph
graph = tlp.newGraph()
graph.setName(model.getId())

# create all graph properties for new tulip graph
def create_props(graph):
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

create_props(graph)

# Aliasing function (clone ubiquitous nodes)
def aliasNodes(g, ub_list):
    """
    Alias (clone) a list of nodes in a graph

    - *g* a tulip/mimoza graph
    - *ub_list* a list of nodes
    """
    ub = g.getBooleanProperty("ub")
    for n in g.getNodes():
        if n in ub_list:
            ub[n] = True

aliasNodes(graph, ub_list)

# run Mimoza functions for further graph manipulaption
root, c_id2info, c_id2outs, ub_sps = sbml2tlp.import_sbml(model, groups_sbml, ub_list, c_id2level, graph)
fc, (n2lo, e2lo), hidden_c_ids, c_id_hidden_ubs, graph =\
    tulip_graph2geojson.graph2geojson(c_id2info, c_id2outs, root)

# save raw Mimoza graph output as backup
tlp.saveGraph(graph, "{}_raw.tlp".format(model_id))


####################################################
# change basic graph layout and cosmetic properties#
####################################################

# add labels to nodes
def addLabels(g, node_type, label_type):
    """
    Add labels to nodes and/or edges of graph

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

# define label targets
NODE_TYPES = {'metabolites': 1,
              'reactions'  : 2
              }

addLabels(graph, NODE_TYPES['metabolites'], "id")
addLabels(graph, NODE_TYPES['reactions'], "id")

# change node or edge colour
def setColours(g, item_list, colour):
    """
    Apply colour to a list of edges

     - *g* a tulip/mimoza graph
     - *item_list* a list of edges (graph.getEdges) or nodes (graph.getNodes)
     - *colour* a predefined colour
    """
    # apply colour to graph
    viewColor = g.getColorProperty("viewColor")
    for e in item_list:
        viewColor[e] = colour

# define colours
COLOURS = {'red'  : tlp.Color(255, 0, 0),
           'green': tlp.Color(0, 255, 0),
           'blue' : tlp.Color(0, 0, 255),
           'grey' : tlp.Color(220, 220, 220)
           }
setColours(graph, graph.getEdges(), COLOURS['red'])

# change line types (Polyline, BezierCurve, CarmullRomCurve, CubicSplineCurve)
viewShape = graph.getIntegerProperty("viewShape")

for e in graph.getEdges():
    viewShape[e] = tlp.EdgeShape.BezierCurve

# different layouts on pathways (subsystems)
# manually define input pathways
small_circle = ["M1", "M2", "M9", "R02", "R10", "R11"]
larger_circle = ["M3", "M4", "M5", "M6", "M7", "R04", "R05", "R06", "R07", "R08"]

# get input pathway node list from a certain compartment
#comp_ids = cmod.getCompartmentIds() # get a list of all compartments
#comp_2 = cmod.getCompartment(comp_ids[1]) # get  compartment 1 (first entry in list of all compartments)
#comp_2_species = comp_2.containsSpecies() # list of species in compartment 1

# define function to create subgraph from graph and apply layout
def createSubgraph(graph, sg_name, n_list):
    """
    Creates a subgraph from a list of nodes

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
        if id_[n] in n_list:
            selection[n] = True
            sg.addNode(n)
    for e in graph.getEdges():
        if selection[graph.source(e)] and selection[graph.target(e)]:
            selection[e] = True
            sg.addEdge(graph.source(e), graph.target(e))
    # reset selection
    selection.setAllNodeValue(False)
    selection.setAllEdgeValue(False)

    return sg

# define function to apply layout to subgraph and center it
def layoutSubgraph(sg, sg_layout, x_center = None, y_center = None, z_center = 0):
    """
    Apply layout to subgraph and center it around a given point.

    - *sg* a tulip/mimoza subgraph
    - *sg_layout* subgraph layout (string), ex: "Circular (OGDF)"
    - *x_center* x coordinate of layout center which defines subgraph position (optional)
    - *y_center* y coordinate of layout center which defines subgraph position (optional)
    - *z_center* z coordinate of layout center which defines subgraph position, default 0
    """
    # pseudo-fixed layout method
    # list of available algorithms: tlp.getLayoutAlgorithmPluginsList()
    # apply algorithm to subgraph
    sg.applyLayoutAlgorithm(sg_layout)
    # center layout: center(viewLayout, newCenter, subgraph)
    # newCenter can be given as coordinates (tlp.Coord(x,y,z)) or
    viewLayout = graph.getLayoutProperty("viewLayout")
    if x_center is not None and y_center is not None:
        tlp.LayoutProperty.center(viewLayout, tlp.Coord(x_center, y_center, z_center), sg)

    return sg

# call functions only if toymodel is loaded
if sbml == os.path.join(cDir, "TestModels", "toy_model_biofuel.xml"):
    # create toymodel subgraphs and apply centered layout
    small_circle_subgraph = createSubgraph(graph, "small_circle", small_circle)
    layoutSubgraph(small_circle_subgraph, "Circular (OGDF)", 50, 75, 0)

    larger_circle_subgraph = createSubgraph(graph, "larger_circle", larger_circle)
    layoutSubgraph(larger_circle_subgraph, "Circular (OGDF)", -50, 25, 0)

    # save created toymodel subgraphs
    tlp.saveGraph(small_circle_subgraph, "{}_SmallCircle.tlp".format(model_id))
    tlp.saveGraph(larger_circle_subgraph, "{}_LargerCircle.tlp".format(model_id))

# fix position of a specific node (M1 (id 2) for small circle, M3 (id 8) for larger circle)
#viewLayout[tlp.node(2)] = tlp.Coord(50,75,0)
#viewLayout[tlp.node(8)] = tlp.Coord(-50,25,0)

# layout all species from compartment 2 (extracellular) seperately along the edges of the graph
#comp_2_sg = createSubgraph(graph, "comp_2", comp_2_species)
#layoutSubgraph(comp_2_sg, "Hierarchical Graph", 0, 0, 0)


#######################
# visualize FBA result#
#######################

# color nodes and edges depending on flux value (pos, neg, 0) -> do subgraphs for each (three colors)
# input format: dict = {"M1": -3, ...}
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

# save total graph with all changes
tlp.saveGraph(graph, "{}_graph.tlp".format(model_id))


#########################################################
# generate SBGN and SBML+Layout Files (Mimoza functions)#
#########################################################

# add layout to SBML file
groups_document = reader.readSBML(groups_sbml)
groups_model = groups_document.getModel()
gen_model = False
sbml_manager.save_as_layout_sbml(groups_model, gen_model, groups_sbml, n2lo, ub_sps)

# create SBGN file
groups_sbgn = os.path.join(cDir, '{}.sbgn'.format(model_id))
sbgn_helper.save_as_sbgn(n2lo, e2lo, groups_model, groups_sbgn)