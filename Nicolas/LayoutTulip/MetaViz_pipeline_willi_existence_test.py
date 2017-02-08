import os
import MetaViz

##################
# 0.) User Input #
##################

# 1. global control variables
# activate/deactivate parts of the pipeline when starting the script from a console (including the one in PyCharm)
global_ops = {"AUTO_ALIAS": True,  # apply aliasing to every node of the graph with a certain connectivity
                                   # and option to alias biomass (result in maingraph)
                                   # if auto alias is active the automatically generated layout will still be
                                   # overwritten if any other layouting functions are active aswell
              "AUTO_ALIAS_OPTIONS": {"spec_split": 4, "reac_split": 2},  # threshold values for species and reac
                                                                          # splitting

              "PATHWAY_LISTS": True,  # layout annotated groups of sbml file as seperate pathways
              "LAYOUT_PATHWAYS": False,  # needs PATHWAY_LISTS true, so that necessary node lists are generated
              "ADD_COFACTORS": ["Respiratory chain"],  # add nodes for small sps/cofactors after they have been
              # deleted in path_dict

              "TRANSPORT_RADIAL": False,  # apply Tree Radial algorithm to transport reactions in output_whole
              "SELECTED_PATHWAYS": ["objective reaction",  # pathways to be saved in output_selected
                                    "Transport"],
              "GRID_OUTPUT": ["Glycolysis",  # align pathways on a grid, control position
                              "TCA cycle",
                              "Respiratory chain",
                              "Aminosugar Metabolism"],

              "EXPORT_TLP": True,  # save generated graphs as .tlp files
              "EXPORT_GML": True,  # export generated graphs as .gml files
              "EXPORT_SVG": True,  # export generated graphs as .svg files
              "EXPORT_JSON": True,  # export generated graphs as .json files
              "EXPORT_SBML": True,  # export generated graphs as .sbml files (in layout plugin)
              "EXPORT_SBGN": True,  # export generated graphs as .sbgn files

              "COLOR_FBA": True,  # color the reaction nodes and in and outgoing edges
              "COLOR_FBA_COLORBLIND": False,  # based on flux values obtained from FBA
              "COLOR_PATHWAYS": False  # pathways specific coloring
              }
# 2. read and initialize SBML file with libsbml and cbmpy, do FBA
# get current directory
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
## bgoli 180117: organising our tests
## bgoli 180117: other tests
# sbml = os.path.join(cDir, "TestModels", "others", "toy_model_biofuel.xml")
# sbml = os.path.join(cDir, "TestModels", "others", "e_coli_core.xml")
# sbml = os.path.join(cDir, "TestModels", "others", "Ecoli_Carlson2003.xml")
# sbml = os.path.join(cDir, "TestModels", "others", "Hpylori_iIT341.xml")
# sbml = os.path.join(cDir, "TestModels", "others", "Lactococcus_lactis_MG1363.xml")
## BiGG2 models (with pathway annotation)
#sbml = os.path.join(cDir, "TestModels", "bigg2", "iIT341.xml")
# sbml = os.path.join(cDir, "TestModels", "bigg2", "iJO1366.xml") # newest E. coli
# sbml = os.path.join(cDir, "TestModels", "bigg2", "iAF1260.xml") # E. coli about 2000+ reactions
# sbml = os.path.join(cDir, "TestModels", "bigg2", "iJR904.xml") # E. coli about 1000 reactions
# sbml = os.path.join(cDir, "TestModels", "bigg2", "iND750.xml") # Yeast about 1000 reactions
## Willis models
# sbml = os.path.join(cDir, "TestModels", "willi", "M5_6.xml")
#sbml = os.path.join(cDir, "TestModels", "willi", "M5_6_bigggroups.xml")
#sbml = os.path.join(cDir, "TestModels", "willi", "M5_6_compgroups.xml")
## Willis existence test
# whole E.coli
sbml = os.path.join(cDir, "TestModels", "willi", "(bigg1-eco)-(Ecoli_iAF1260.glc.l3).seqplus.xml")
# subset E.coli
sbml2 = os.path.join(cDir, "TestModels", "willi", "322b_EcoliiAF1260.xml")
# whole Plantarum
# sbml = os.path.join(cDir, "TestModels", "willi", "(lab-lplawcfs1)-(Lactobacillus_plantarum_WCFS1).seqplus.xml")
# subset Plantarum
# sbml = os.path.join(cDir, "TestModels", "willi", "322b_lplantarum.xml")
## Stijns models
# sbml = os.path.join(cDir, "TestModels", "stijn", "Y7.xml") #M.getNumReactions: 3498, M.getNumSpecies: 2220
# Error message when trying to run Metaviz with Y7xml:
"""
Traceback (most recent call last):

File "C:/Users/Nicolas/surfdrive/LabR1_Teusink_Amsterdam/Project/FindStructuredLayout/MetaViz.py", line 1174, in <module>
graph, graph_props)
File "C:/Users/Nicolas/surfdrive/LabR1_Teusink_Amsterdam/Project/FindStructuredLayout/MetaViz.py", line 100, in layoutModel
root)
File "C:\Users\Nicolas\surfdrive\LabR1_Teusink_Amsterdam\Project\FindStructuredLayout\tulip_graph2geojson.py", line 352,      in graph2geojson
process_compartments(c_id2info, meta_graph, n2xy)
File "C:\Users\Nicolas\surfdrive\LabR1_Teusink_Amsterdam\Project\FindStructuredLayout\tulip_graph2geojson.py", line 457,      in process_compartments
layout(meta_graph)
File "C:\Users\Nicolas\surfdrive\LabR1_Teusink_Amsterdam\Project\FindStructuredLayout\ubiquitous_layout.py", line 382, in     layout
layout_components(gr)
File "C:\Users\Nicolas\surfdrive\LabR1_Teusink_Amsterdam\Project\FindStructuredLayout\layout_utils.py", line 115, in          layout_components
iter((followers_rev_or_gen, followers_rev, followers_irrev, followers_irrev_gen_only)))
File "C:\Users\Nicolas\surfdrive\LabR1_Teusink_Amsterdam\Project\FindStructuredLayout\layout_utils.py", line 102, in          process_cc
if dfs(cc[0], gr, set(), None, cycle_number_threshold, follower_method) \

File "C:\Users\Nicolas\surfdrive\LabR1_Teusink_Amsterdam\Project\FindStructuredLayout\layout_utils.py", line 82, in <lambda>
followers_rev_or_gen = lambda n, gr: (gr.opposite(e, n) for e in gr.getInOutEdges(n) if n == gr.source(e)
RuntimeError: maximum recursion depth exceeded while calling a Python object

Process finished with exit code 1
"""
# TODO: maximum recursion depth (999?) reached -> change recursive function dfs in layout_utils.py to iterative
# TODO: function?
# Y5 with fewer reactions works so its most likely an issue with big models!
# sbml = os.path.join(cDir, "TestModels", "stijn", "Y5.xml") # M.getNumReactions: 2110, M.getNumSpecies: 1655


########################
# I.) MetaViz Pipeline #
########################

# 1. create Tulip graph and generate graph objects
graph, graph_props = MetaViz.createInitialGraph("initialGraph")
print("Initial (empty) Tulip graph created. \n")

# 2. take sbml model file and initial graph and run all layouting function calls to return necessary data structures
cmod, fba_results, path_dict, COLORS, ub_sps, n2lo, e2lo, groups_sbml, obj_re_nodes, graph_registry, cDir, \
output_path = MetaViz.layoutModel(global_ops, sbml, graph, graph_props)
print("Layout generated. Subgraph(s) added to initial graph. Output graph(s) created from subgraph(s). \n")


## existence test
# read sub model from sbml2
# TODO: write a function that can take in a list of sbml files and output a dict with keys=reac_ids and values=diff
# TODO: numbers for diff conditions
# TODO: Multiploting function
cmod2 = MetaViz.cbmpy.readSBML3FBC(sbml2)
existence_dict = {}
# get reaction id list from sub model cmod2
for re in cmod2.getReactionIds():
    # test if reaction is also in main model cmod
    # if True: 1
    if re in cmod.getReactionIds():
        existence_dict[re] = 1
    # if False: 0
    else:
        existence_dict[re] = 0


# 3. apply coloring algorithms to graphs
# pass existence_dict instead of fba_results to coloring function
MetaViz.colorGraphs(COLORS, global_ops, existence_dict, path_dict, graph_registry, output_path, cmod.getId())
print("Graph(s) colored. \n")

# 4. export graph SVG, GML, JSON, SBML+layout Plugin and/or SBGN file (controlled by global_ops)
MetaViz.exportGraph(output_path, graph_registry, global_ops, cDir, cmod.getId(), groups_sbml, n2lo, e2lo, ub_sps)


###################################
# II.) SQLite Database Operations #
###################################

# reason for database: to have a tool-independent way of getting network map info -> feed into d3, networkx,
# cytoscape, etc...!


# 1. create SQLite db to store a graph object
db_path = MetaViz.buildSQLiteDatabase(graph_registry['output_whole'][0], output_path)

# DB INFO ##########################################################################################################
# 5 tables in database:
# - graph_nodes: every row is a node, each column is a node property
# - graph_edges: every row is an edge, each column is an edge property
# - node_connectivity: every row is a node, first column contains node ids and second column contains lists of
#                      connected edges
# - edge connectivity: every row is an edge, first column contains edge ids, second column contains source node ids
#                      and third column contains target node ids
# - edge color: every row is an edge, first column contains edge id, subsequent columns contain color sets (
#                      tuples of (flux value, pickled tulip color prop))
####################################################################################################################

# 2. retrieve graph data from SQL database
db_graph_node_list, node_col_headers = MetaViz.retrieveDataFromDatabaseTable(db_path, 'graph_nodes')
db_graph_edge_list, edge_col_headers = MetaViz.retrieveDataFromDatabaseTable(db_path, 'graph_edges')
db_graph_edgeconn_list, edgeconn_col_headers = MetaViz.retrieveDataFromDatabaseTable(db_path, 'edge_connectivity')

# 3. SQL --> back to Tulip (to see if everything is preserved)
db_graph = MetaViz.buildTulipGraphFromSQL(db_graph_node_list, db_graph_edge_list, db_graph_edgeconn_list,
                                  node_col_headers, edge_col_headers, "graph_fromDB", output_path)


# TODO: SQL --> d3py/networkx/someother gaph program
# TODO: extract only nodes with coords, and ids from graph_nodes table and connectivity between the nodes
# TODO: use networkx, and/or D3 (d3py in python) to rebuild a part of the graph only from SQL database output
# TODO: do the networkx/d3 thing on a subset of SQL output (only nodes of aspartate pathway, etc.)
