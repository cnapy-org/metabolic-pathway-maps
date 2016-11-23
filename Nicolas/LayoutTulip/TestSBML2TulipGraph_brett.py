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


######################################
# read SBML and run Mimoza for layout#
######################################

# read SBML file (with libsbml)
#sbml = os.path.join(cDir, "TestModels", "toy_model_biofuel.xml") # 4, 4, False
sbml = os.path.join(cDir, "TestModels", "Ecoli_Carlson2003.xml") # 6, 4, True
#sbml = os.path.join(cDir, "TestModels", "Hpylori_iIT341.xml") # 6, 4, True
#sbml = os.path.join(cDir, "TestModels", "iJO1366.xml") # 6, 4, True
#sbml = os.path.join(cDir, "TestModels", "M5_6.xml") # 4, 4, True

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

#ub_sp_list = ['ATP', 'ADP', 'M8', 'M10', 'M11', 'M12', 'M13', 'biomass', 'S', 'fuel1', 'fuel2']
#ub_re_list = ['R09', 'R12', 'R13', 'R14', 'R15', 'R16', 'R17', 'R18']

ub_sp_list, ub_re_list = autoSplitConnectedNodes(cmod, 4, 4, True)


print(ub_sp_list)
print(ub_re_list)

# create groups
cmod.createGroup('g_ubiquitous_sps')
g1 = cmod.getGroup('g_ubiquitous_sps')
g1.setName('ubiquitous species')
ub_list = ub_sp_list + ub_re_list
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

# create Tulip Graph (with Mimoza)
root, c_id2info, c_id2outs, ub_sps = sbml2tlp.import_sbml(model, groups_sbml, ub_list)
fc, (n2lo, e2lo), hidden_c_ids, c_id_hidden_ubs, graph =\
    tulip_graph2geojson.graph2geojson(c_id2info, c_id2outs, root)

# save total graph and different subgraphs
tlp.saveGraph(graph, "test_raw.tlp".format(model_id))

###########################################
# change graph layout with tulip functions#
###########################################

# backup original graph
backup_graph = graph.addCloneSubGraph("backup_graph")
tlp.saveGraph(backup_graph, "{}_backup_graph.tlp".format(model_id))

# add labels to nodes
viewLabel = graph.getStringProperty("viewLabel")
name = graph.getStringProperty("name")

for n in graph.getNodes():
    viewLabel[n] = name[n]

# change global line colour
def setEdgeColours(g, edge_list, colour):
    """
    Apply colour to a list of edges

     - *g* a tulip/mimoza graph
     - *edge_list* a list of edges
     - *colour* a predefined colout
    """
    # apply colour to graph
    viewColor = g.getColorProperty("viewColor")
    for e in edge_list:
        viewColor[e] = colour

# define colours
COLOURS = {'red'  : tlp.Color(255, 0, 0),
           'green': tlp.Color(0, 255, 0),
           'blue' : tlp.Color(0, 0, 255)
           }
setEdgeColours(graph, graph.getEdges(), COLOURS['red'])

# change line types (Polyline, BezierCurve, CarmullRomCurve, CubicSplineCurve)
viewShape = graph.getIntegerProperty("viewShape")

for e in graph.getEdges():
    viewShape[e] = tlp.EdgeShape.BezierCurve

# different layouts on subsystems
# create empty subgraphs
small_circle_subgraph = graph.addSubGraph("SmallCircle")
larger_circle_subgraph = graph.addSubGraph("LargerCircle")

# define properties
nodeType = graph.getIntegerProperty("type")
name = graph.getStringProperty("name")
selection = graph.getBooleanProperty("viewSelection")

# reset selection
selection.setAllNodeValue(False)
selection.setAllEdgeValue(False)

# input pathways
small_circle = ["M1", "M2", "M9", "R02", "R10", "R11"]
larger_circle = ["M3", "M4", "M5", "M6", "M7", "R04", "R05", "R06", "R07", "R08"]

# populate subgraphs
# small circle
for n in graph.getNodes():
    if name[n] in small_circle:
        selection[n] = True
        small_circle_subgraph.addNode(n)
for e in graph.getEdges():
    if selection[graph.source(e)] and selection[graph.target(e)]:
        selection[e] = True
        small_circle_subgraph.addEdge(graph.source(e), graph.target(e))

# reset selection
selection.setAllNodeValue(False)
selection.setAllEdgeValue(False)

# larger circle
for n in graph.getNodes():
    if name[n] in larger_circle:
        selection[n] = True
        larger_circle_subgraph.addNode(n)
for e in graph.getEdges():
    if selection[graph.source(e)] and selection[graph.target(e)]:
        selection[e] = True
        larger_circle_subgraph.addEdge(graph.source(e), graph.target(e))

# reset selection
selection.setAllNodeValue(False)
selection.setAllEdgeValue(False)

# pseudo-fixed layout method
# apply layout algorithms to subgraphs
# list of available algorithms: tlp.getLayoutAlgorithmPluginsList()
small_circle_subgraph.applyLayoutAlgorithm("Circular (OGDF)")
larger_circle_subgraph.applyLayoutAlgorithm("Circular (OGDF)")

# center layout: center(viewLayout, newCenter, subgraph)
# newCenter can be given as coordinates (tlp.Coord(x,y,z)) or
viewLayout = graph.getLayoutProperty("viewLayout")

tlp.LayoutProperty.center(viewLayout, tlp.Coord(50,75,0), small_circle_subgraph)
tlp.LayoutProperty.center(viewLayout, tlp.Coord(-50,25,0),larger_circle_subgraph)

# fix position of a specific node (M1 (id 2) for small circle, M3 (id 8) for larger circle)
#viewLayout[tlp.node(2)] = tlp.Coord(50,75,0)
#viewLayout[tlp.node(8)] = tlp.Coord(-50,25,0)

# reduce edge crossing
# in progress...

# save total graph and different subgraphs
tlp.saveGraph(graph, "{}_graph.tlp".format(model_id))
tlp.saveGraph(small_circle_subgraph, "{}_SmallCircle.tlp".format(model_id))
tlp.saveGraph(larger_circle_subgraph, "{}_LargerCircle.tlp".format(model_id))


######################################
# generate SBGN and SBML+Layout Files#
######################################

# add layout to SBML file
groups_document = reader.readSBML(groups_sbml)
groups_model = groups_document.getModel()
gen_model = False
sbml_manager.save_as_layout_sbml(groups_model, gen_model, groups_sbml, n2lo, ub_sps)

# create SBGN file
groups_sbgn = os.path.join(cDir, '{}.sbgn'.format(model_id))
sbgn_helper.save_as_sbgn(n2lo, e2lo, groups_model, groups_sbgn)



def applyCircularOGDF(g, params={}):
    """
    # params['minDistCircle'] = ...
    # params['minDistLevel'] = ...
    # params['minDistSibling'] = ...
    # params['minDistCC'] = ...
    # params['pageRatio'] = ...

    """

    sg_params = tlp.getDefaultPluginParameters('Circular (OGDF)', g)
    print(sg_params)
    for k in params:
        if k in sg_params:
            sg_params[k] = params[k]
    success = g.applyLayoutAlgorithm('Circular (OGDF)', sg_params)
    print('applyCircularOGDF: {}'.format(success))
