import os
import cbmviz

# User Input
global_ops = cbmviz.cfg_default
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
sbml = os.path.join(cDir, "TestModels", "bigg2", "iIT341.xml") # H.pylori
file_name = "iIT341"
output_path = os.path.join(cDir, "{}_Output".format(file_name))

# load a graph object from an already existing layout (e.g. a JSON file)
existing_layout = os.path.join('iIT341_Output', 'TCAcycle.json')
input_graph = cbmviz.loadGraphObjFromFile(existing_layout)
# generate a dictionary from the graph object
input_graph_dict = cbmviz.createGraphDictFromGraphObj(input_graph)

# test: color edges based on group membership instead of flux values
global_ops.update({"COLOR_FBA": False})
global_ops.update({"COLOR_PATHWAYS": True})

# read SBML file to get info to build path dict which will be used for pathway coloring
cmod, stM, obj_re_id  = cbmviz.createCBMPYModel(sbml)
path_dict, small_sps_dict = cbmviz.automPathDict(cmod, global_ops)
fba_results = cbmviz.createPathColoringDict(path_dict)

# run coloring and exporting functions
cbmviz.colorGraphs(global_ops, fba_results, input_graph_dict, output_path, "iIT341")
cbmviz.exportGraph(output_path, input_graph_dict, global_ops, cDir, "iIT341")
