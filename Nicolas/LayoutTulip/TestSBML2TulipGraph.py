# Test transition from sbml file to structured tulip graph, sbml with layout and sbgn file
# necessary  dependencies (Wind64 Anaconda): pip install geojson tarjan libsbgnpy tulip-python

import os
# get current directory
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))

import libsbml
import cbmpy
from tulip import tlp

# core Mimoza modules
import sbml2tlp
import sbml_manager
import tulip_graph2geojson
import sbgn_helper
import graph_properties
import resize

#############################################
# read SBML with libsbml and cbmpy, do FBA  #
#############################################

# read SBML file (with libsbml)
#sbml = os.path.join(cDir, "TestModels", "toy_model_biofuel.xml") # 4, 4, False
#sbml = os.path.join(cDir, "TestModels", "e_coli_core.xml") # 6, 4, True
sbml = os.path.join(cDir, "TestModels", "Ecoli_Carlson2003.xml") # 6, 4, True
#sbml = os.path.join(cDir, "TestModels", "Lactococcus_lactis_MG1363.xml") # 6, 4, True
#sbml = os.path.join(cDir, "TestModels", "Hpylori_iIT341.xml") # 6, 4, True
#sbml = os.path.join(cDir, "TestModels", "iJO1366.xml") # 6, 4, True

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

#ub_sp_list, ub_re_list = identifyConnectedNodes(cmod, 4, 4, False)

# 2.) create alias node list based on biological pathways <-> group info in sbml - if available (strategy 2)
# create lists of species ids and reaction ids for each pathway (not necessary if provided by sbml)

# specific to carlsson
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
def identifyBoundarySpecies(sps_list, re_list, cmod):
    """
    Identifies boundary species of given pathway by:
    1.) getting the nodes that have at least one reaction outside the pathway -> level0_nodes (still themselves inside
    the pathway - we don't want to alias them!)
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

#glyc_sps_bound = identifyBoundarySpecies(glyc_sps, glyc_re, cmod)
tca_sps_bound_tup = identifyBoundarySpecies(tca_sps, tca_re, cmod)
# put only the first item of the tuples (species) in a new list and only if it is not part of the glycolysis list
# we don't want to alias to much later on
tca_sps_bound = list(set([tup[0] for tup in tca_sps_bound_tup if tup[0] not in glyc_sps]))

# select biomass reaction
#for r in cmod.getActiveObjectiveReactionIds():
    # add biomass reaction with substrates and products to alias lists (ub_sp_list, ub_re_list)?
    # -> commented because some species of tca cycle are also part of biomass (eg: AKG)
    # => aliasing them would break up the tca circle layout
    #ub_sp_list += [s for s in cmod.getReaction(r).getSpeciesIds() if s not in ub_sp_list]

ub_sp_list = tca_sps_bound
#ub_sp_list = []
ub_re_list = []

# 3.) create alias node list based on a certain compartment <-> compartment info in sbml (strategy 3)
#comp_ids = cmod.getCompartmentIds() # get a list of all compartments
#comp_2 = cmod.getCompartment(comp_ids[1]) # get  compartment 2 (first entry in list of all compartments)
#comp_2_species = comp_2.containsSpecies() # list of species in compartment 2
#comp_2_reac = comp_2.containsReactions() # list of reactions in compartment 2
#ub_sp_list = comp_2_species
#ub_re_list = comp_2_reac

# 4.) manually provide alias node list <-> hard coded by user (for testing purposes)
#ub_sp_list = ['ATP', 'ADP', 'M8', 'M10', 'M11', 'M12', 'M13', 'biomass', 'S', 'fuel1', 'fuel2']
#ub_re_list = ['R09', 'R12', 'R13', 'R14', 'R15', 'R16', 'R17', 'R18']

# combine alias species and alias reaction lists (for aliasNodes function)
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
def create_props(graph):
    # create properties
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

    # return list with variables pointing to the properties
    #prop_list = [anc_term, anc_id, anc_name, term, c_id, id_, name, type_, rev, stoich, tr, ub, viewLayout,
    #             viewShape, viewSize, fake, clone_id, formula, rel_c_id]
    #return prop_list
create_props(graph)

# flag (true/false) a list of nodes (already saved as ub_list) so that they can be aliased later
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

markUbiquitousNodes(graph, ub_list)

# create nodes for all model species
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

id2n = speciesToNodes(graph, model, ub_list)

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

#function to alias any chosen node
def aliasNode(g, n_alias, target_string, x ,y, n_new=None):
    """
    Alias a specific node (reaction nodes and/or metabolite nodes)::

    - *g* a tulip/mimoza graph
    - *n_alias* a tulip node to alias (format: tlp.node(index))
    - *target_string* a name string of the target connected to n_alias
    - *x* x coordinate (int) of new node position
    - *y* y coordinate (int) of new node position
    - *n_new* newly created tulip node (cloned from n_alias)

    """
    viewName = g.getStringProperty("name")
    # iterate over all incoming and outgoing edges of chosen node
    for e in g.getInOutEdges(n_alias):
        # for each node take the edge that either begins at or ends et the chosen target
        if viewName[g.target(e)] == target_string or viewName[g.source(e)] == target_string:
            g.delEdge(e) # delete the edge
    # delete alias node if its not connected to any edge anymore
    if len(list(g.getInOutEdges(n_alias))) == 0:
        g.delNode(n_alias)

    # add new node with same properties as AKG
    n_new = g.addNode()
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
        clone_id = graph.getStringProperty("clone_id")
        clone_id[n_new] = target_string

    # add new edge for new node
    # select new node via tlp.node(#) internal id
    # graph.addEdge(source = new node, target = node that belongs to the name in target_string)
    for n in g.getNodes():
        if viewName[n] == target_string:
            n_edge = g.addEdge(n_new, n)
            # add edge name (same as name of old node)
            viewName[n_edge] = viewName[g.source(n_edge)]

    # adjust props to those of nodes connected to target
        # for n in graph.getInOutNodes(tlp.node(59)):
        #    print(viewName[n])
        # prop_obj[n_new] = prop_obj[tlp.node()]
        #
        # give node different size to old one
        # if prop_string == "viewSize":
        #    prop_obj = graph.getProperty(prop_string)
        #    prop_obj[n_new] = tlp.Size(2, 2, 0)
        #

    return n_new

# delete a set of nodes from the graph
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

# add labels to nodes
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

# define node types for labeling
NODE_TYPES = {'metabolites': 1,
              'reactions'  : 2
              }

addLabels(graph, NODE_TYPES['metabolites'], "id")
addLabels(graph, NODE_TYPES['reactions'], "id")

# change node or edge colour
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

# define colours
COLOURS = {'red'  : tlp.Color(255, 0, 0),
           'green': tlp.Color(0, 255, 0),
           'blue' : tlp.Color(0, 0, 255),
           'grey' : tlp.Color(220, 220, 220),
           'black': tlp.Color(0, 0, 0)
           }

# set colour of all edges (complete graph)
setColours(graph, graph.getEdges(), COLOURS['grey'])

# change line types (Polyline, BezierCurve, CatmullRomCurve, CubicSplineCurve)
viewShape = graph.getIntegerProperty("viewShape")
for e in graph.getEdges():
    viewShape[e] = tlp.EdgeShape.BezierCurve


##################################################################
# different layouts on subsystems (compartments, pathways, etc.) #
##################################################################

# define function to create subgraph from graph and apply layout
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
            sg.addEdge(graph.source(e), graph.target(e))
    # reset selection
    selection.setAllNodeValue(False)
    selection.setAllEdgeValue(False)

    return sg

# define function to apply layout to subgraph and center it
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
    if sg_layout is not None:
        sg.applyLayoutAlgorithm(sg_layout)
    # center layout: center(viewLayout, newCenter, subgraph)
    # newCenter can be given as coordinates (tlp.Coord(x,y,z)) or
    viewLayout = graph.getLayoutProperty("viewLayout")
    if x_center is not None and y_center is not None:
        tlp.LayoutProperty.center(viewLayout, tlp.Coord(x_center, y_center, z_center), sg)

    return sg

# 1.) manually define subsystem node lists
small_circle = ["M1", "M2", "M9", "R02", "R10", "R11"]
larger_circle = ["M3", "M4", "M5", "M6", "M7", "R04", "R05", "R06", "R07", "R08"]

# subsystems "small circle" and "larger circle" are specific to the toymodel
if sbml == os.path.join(cDir, "TestModels", "toy_model_biofuel.xml"):
    # create toymodel subgraphs and apply centered layout
    small_circle_subgraph = createSubgraph(graph, "small_circle", small_circle)
    larger_circle_subgraph = createSubgraph(graph, "larger_circle", larger_circle)
    layoutSubgraph(small_circle_subgraph, "Circular (OGDF)", 100, 100, 0)
    layoutSubgraph(larger_circle_subgraph, "Circular (OGDF)", -50, 25, 0)
    # save created toymodel subgraphs
    tlp.saveGraph(small_circle_subgraph, "{}_SmallCircle.tlp".format(model_id))
    tlp.saveGraph(larger_circle_subgraph, "{}_LargerCircle.tlp".format(model_id))

# 2.) get subsystem node list from all species of a certain compartment
#comp_ids = cmod.getCompartmentIds() # get a list of all compartments
#comp_2 = cmod.getCompartment(comp_ids[1]) # get  compartment 1 (first entry in list of all compartments)
#comp_2_species = comp_2.containsSpecies() # list of species in compartment 1

# layout all species from compartment 2 (extracellular) seperately along the edges of the graph
#comp_2_sg = createSubgraph(graph, "comp_2", comp_2_species)
#layoutSubgraph(comp_2_sg, "Hierarchical Graph", 0, 0, 0)

# 3.) get subsystem node list from different biological pathways (e.g.: glycolysis, tca and 'not specified')
glyc_sg = createSubgraph(graph, "glycolysis", glyc_sps)
tca_sg = createSubgraph(graph, "tca_cycle", tca_sps)
not_spec_sg = createSubgraph(graph, "not specified", not_spec_sps)

setColours(graph, glyc_sg.getNodes(), COLOURS['red'])
setColours(graph, tca_sg.getNodes(), COLOURS['blue'])

#layoutSubgraph(glyc_sg, "Mixed Model", -45, 90, 0)
#layoutSubgraph(tca_sg, "Circular", -45, -30, 0)
#layoutSubgraph(not_spec_sg, "Planarization Layout (OGDF)", None, None, 0)

# fix position of a specific node (M1 (id 2) for small circle, M3 (id 8) for larger circle)
#viewLayout[tlp.node(2)] = tlp.Coord(50,75,0)
#viewLayout[tlp.node(8)] = tlp.Coord(-50,25,0)
#viewLayout = graph.getLayoutProperty("viewLayout")
#viewLayout[tlp.node(16)] = tlp.Coord(100, 100, 0)


####################################################
# aliasing examples (using the function aliasNode) #
####################################################

# example 1: remove edge (edge 49) connecting biomass to TCA cycle at AKG (node 4)
#            get the edge that connects AKG (node 4) to biomass (R70, node 60)
#aliasNode(graph, tlp.node(4), "R70", 22, 27)

# example 2: cut out the biomass function (individually alias all of its connected species)
#            and move it together with all of its now aliased species somewhere else
#            biomass function = R70 (node 60)

import random

# alias the biomass reaction
viewName = graph.getStringProperty("name")
viewLayout = graph.getLayoutProperty("viewLayout")
biomass_node_list = []

for n in graph.getNodes():
    # select biomass reaction node
    if viewName[n] == "R70":
        # add biomass reaction node to list
        biomass_node_list.append(n)
        # access biomass reaction node coords (to layout connected nodes around it)
        x, y, z = viewLayout[n]
        # iterate over all nodes which are connected to biomass reaction node
        for node in graph.getInOutNodes(n):
            # alias node and move it to a random position close to biomass reaction node (quick test hack)
            node_aliased = aliasNode(graph, node, "R70", x+random.random()*50, y+random.random()*50, None)
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
biomass_sg = createSubgraph(graph, "biomass reaction", biomass_node_list)

print(re_node_coords)

# increase circle radius by increasing size of nodes?
viewSize = graph.getSizeProperty("viewSize")

for n in biomass_sg.getNodes():
    viewSize[n] = tlp.Size(3, 3, 0)

# apply layout to biomass subgraph
layoutSubgraph(biomass_sg, "Circular", re_node_coords[0], re_node_coords[1], 0)


########################
# visualize FBA result #
########################
"""
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
"""


############################################################################
# save final graph, generate SBGN and SBML+Layout Files (Mimoza functions) #
############################################################################

# save total graph with all changes
tlp.saveGraph(graph, "{}_graph.tlp".format(model_id))

# add layout to SBML file
groups_document = reader.readSBML(groups_sbml)
groups_model = groups_document.getModel()
gen_model = False
sbml_manager.save_as_layout_sbml(groups_model, gen_model, groups_sbml, n2lo, ub_sps)

# create SBGN file
groups_sbgn = os.path.join(cDir, '{}.sbgn'.format(model_id))
sbgn_helper.save_as_sbgn(n2lo, e2lo, groups_model, groups_sbgn)
