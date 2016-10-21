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
# add ubiquitous species to group
ub_sps = ['ATP', 'ADP']
ub_sps = [cmod.getSpecies(i) for i in ub_sps]
g1.addMember(ub_sps)
# add RDF annotation to group:
# MIRIAM annotation and SBO term (important for drawing ub species seperately)
g1.setSBOterm('SBO:0000248')
g1.addMIRIAMannotation('isDescribedBy', 'ChEBI', 'CHEBI:123456') #dummy annotation
g1.miriam.isDescribedBy[0] = "http://identifiers.org/ubiquitous" #change dummy to fit mimoza
cbmpy.writeSBML3FBC(cmod, '{}_with_groups.xml'.format(model_id), add_groups=True)
groups_sbml = os.path.join(cDir, "{}_with_groups.xml".format(model_id))

# generalization does not seem to be mandatory

# create Tulip Graph (with Mimoza)
root, c_id2info, c_id2outs, chebi, ub_sps = sbml2tlp.import_sbml(model, groups_sbml)
fc, (n2lo, e2lo), hidden_c_ids, c_id_hidden_ubs, graph =\
    tulip_graph2geojson.graph2geojson(c_id2info, c_id2outs, root)

# add labels to nodes
viewLabel = graph.getStringProperty("viewLabel")
name = graph.getStringProperty("name")

for n in graph.getNodes():
    viewLabel[n] = name[n]

# TEST GRAPH MANIPULATION
# change graph layout with tulip functions




tlp.saveGraph(graph, "{}_graph.tlp".format(model_id))

# add layout to SBML file
groups_document = reader.readSBML(groups_sbml)
groups_model = groups_document.getModel()
gen_model = False
sbml_manager.save_as_layout_sbml(groups_model, gen_model, groups_sbml, n2lo, ub_sps)

# create SBGN file
groups_sbgn = os.path.join(cDir, '{}.sbgn'.format(model_id))
sbgn_helper.save_as_sbgn(n2lo, e2lo, groups_model, groups_sbgn)
