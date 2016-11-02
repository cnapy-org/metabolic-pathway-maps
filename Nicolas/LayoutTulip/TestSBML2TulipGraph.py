# Test transition from sbml file to structured tulip graph, sbml with layout and sbgn file

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

# read SBML file (with libsbml)
sbml = os.path.join(cDir, "TestModels", "toy_model_biofuel.xml")
reader = libsbml.SBMLReader()
doc = reader.readSBML(sbml)
model = doc.getModel()
model_id = model.getId()

# create SBML file with groups (with cbmpy)
cmod = cbmpy.readSBML3FBC(sbml)
cmod.createGroup('g_ubiquitous_sps')
g1 = cmod.getGroup('g_ubiquitous_sps')
g1.setName('ubiquitous species')
# add ubiquitous species to group.
# bgoli: can be any metabolite?
ub_sps = ['ATP', 'ADP']
ub_sps = [cmod.getSpecies(i) for i in ub_sps]
g1.addMember(ub_sps)
# add RDF annotation to group:
# MIRIAM annotation and SBO term (important for drawing ub species seperately)
# bgoli: requires any ChEBI annotation to work, remember to remove in output!
g1.setSBOterm('SBO:0000248')
g1.addMIRIAMannotation('isDescribedBy', 'ChEBI', 'CHEBI:000000') #dummy annotation
g1.miriam.isDescribedBy[0] = "http://identifiers.org/ubiquitous" #change dummy to fit mimoza

cbmpy.writeSBML3FBC(cmod, '{}_with_groups.xml'.format(model_id), add_groups=True)
groups_sbml = os.path.join(cDir, "{}_with_groups.xml".format(model_id))

# create Tulip Graph (with Mimoza)
root, c_id2info, c_id2outs, ub_sps = sbml2tlp.import_sbml(model, groups_sbml)
fc, (n2lo, e2lo), hidden_c_ids, c_id_hidden_ubs, graph =\
    tulip_graph2geojson.graph2geojson(c_id2info, c_id2outs, root)

# change graph layout with tulip functions
# add labels to nodes
viewLabel = graph.getStringProperty("viewLabel")
name = graph.getStringProperty("name")

for n in graph.getNodes():
    viewLabel[n] = name[n]

# change line colour (dashes seem to be not possible)
# add example of different line "types"

# global line colour

# define colours
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

COLOURS = {'red' : tlp.Color(255, 0, 0),
           'green': tlp.Color(0, 255, 0),
           'blue' : tlp.Color(0, 0, 255)
           }

setEdgeColours(graph, graph.getEdges(), COLOURS['red'])

# different layouts on subsystems
# Brett: what possibilities there are for manipulation subgraphs
test = graph.addSubGraph("test")

nodeType = graph.getIntegerProperty("type")
name = graph.getStringProperty("name")
selection = graph.getBooleanProperty("viewSelection")

selection.setAllNodeValue(False)
selection.setAllEdgeValue(False)

selected_nodes = ["M1", "M2", "M9", "R02", "R10", "R11"]

# Brett: look into cloning, nodes, edges, subgraphs, etc.

for n in graph.getNodes():
    if name[n] in selected_nodes:
        selection[n] = True
        test.addNode(n)

for e in graph.getEdges():
    if selection[graph.source(e)] and selection[graph.target(e)]:
        selection[e] = True
        test.addEdge(graph.source(e), graph.target(e))

# list of available layouts: tlp.getLayoutAlgorithmPluginsList()
#test.applyLayoutAlgorithm("Circular (OGDF)")

# reduce edge crossing
# in progress...

# save graph in its final state
tlp.saveGraph(graph, "{}_graph.tlp".format(model_id))
tlp.saveGraph(test, "{}_test.tlp".format(model_id))

# add layout to SBML file
groups_document = reader.readSBML(groups_sbml)
groups_model = groups_document.getModel()
gen_model = False
sbml_manager.save_as_layout_sbml(groups_model, gen_model, groups_sbml, n2lo, ub_sps)

# create SBGN file
groups_sbgn = os.path.join(cDir, '{}.sbgn'.format(model_id))
sbgn_helper.save_as_sbgn(n2lo, e2lo, groups_model, groups_sbgn)
