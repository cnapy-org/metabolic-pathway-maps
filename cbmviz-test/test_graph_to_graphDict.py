import os
import cbmviz

# User Input
global_ops = cbmviz.cfg_default
# turn on pathway coloring
global_ops.update({"COLOR_FBA": False})
global_ops.update({"COLOR_PATHWAYS": True})

cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
sbml = os.path.join(cDir, "TestModels", "bigg2", "iIT341.xml") # H.pylori
file_name = "iIT341"
output_path = os.path.join(cDir, "{}_Output".format(file_name))


# test graph to graph_registry converter function
# -----------------------------------------------
# import a .tlp file and load graph object
test_path = os.path.join(cDir, "OutputExamples", "pathwayLayout_FBAcolored", "iIT341_Output", "iIT341_output_whole.tlp")
test_graph_obj = cbmviz.loadGraphObjFromFile(test_path)

# store the graph object in a new graph_registry dict
# (same structure -> keys: graph names, values: [graph_obj, {graph_prop dict}])
test_graph_registry = cbmviz.createGraphDictFromGraphObj(test_graph_obj)

# use this new test graph_registry dict to do the coloring and export new results
cmod, stM, obj_re_id = cbmviz.createCBMPYModel(sbml)
path_dict, small_sps_dict = cbmviz.automPathDict(cmod, global_ops)
color_input_dict = cbmviz.createPathColoringDict(path_dict)
cbmviz.colorGraphs(global_ops, color_input_dict, test_graph_registry, output_path, file_name)
cbmviz.exportGraph(output_path, test_graph_registry, global_ops, cDir, file_name)


