import os
import cbmviz

# User Input
global_ops = cbmviz.cfg_default
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
sbml = os.path.join(cDir, "TestModels", "bigg2", "iIT341.xml") # H.pylori

# MetaViz Pipeline
graph, graph_props = cbmviz.createInitialGraph("initialGraph")
cmod, fba_results, path_dict, COLORS, ub_sps, n2lo, e2lo, groups_sbml, obj_re_nodes, graph_registry, cDir, \
output_path = cbmviz.layoutModel(global_ops, sbml, graph, graph_props)



# test graph to graph_registry converter function
# -----------------------------------------------
# delete old graph_registry object
del(graph_registry)
# import a .tlp file and load graph object
test_graph_obj = cbmviz.loadGraphObjFromTulipFile("C:\Users\Nicolas\surfdrive\LabR1_Teusink_Amsterdam\Project\CBMViz_TestScripts\iIT341_Output\iIT341_output_whole.tlp")
# store the graph object in a new graph_registry dict
# (same structure -> keys: graph names, values: [graph_obj, {graph_prop dict}])
test_graph_registry = cbmviz.createGraphDictFromGraphObj(test_graph_obj)



# use this new test graph_registry dict to do the coloring and export new results
cbmviz.colorGraphs(COLORS, global_ops, fba_results, path_dict, test_graph_registry, output_path, cmod.getId())
cbmviz.exportGraph(output_path, test_graph_registry, global_ops, cDir, cmod.getId(), groups_sbml, n2lo, e2lo, ub_sps)


