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



# test combine subgraphs function
# -------------------------------
# 1.) iIT341 "TCA cycle" (directly from pipeline)
test_pylori_graph = graph_registry["output_whole"][0]
test_pylori_graph.cbmviz_file_name = "iIT341"
test_pylori_list = [test_pylori_graph, ["TCA cycle", "Glycolysis"]]
# 2.) iJO1366 "Citric Acid Cycle" (from tlp file)
test_ecoli_path = os.path.join(cDir, "OutputExamples", "pathwayLayout_FBAcolored", "iJO1366_Output",
                               "iJO1366_output_whole.tlp")
test_ecoli_graph = cbmviz.loadGraphObjFromFile(test_ecoli_path)
test_ecoli_graph.cbmviz_file_name = "iJO1366"
test_ecoli_list = [test_ecoli_graph, ["Citric Acid Cycle", "Glycolysis/Gluconeogenesis"]]
# 3.) IND750 "Citric Acid Cycle" (from tlp file)
test_yeast_path = os.path.join(cDir, "OutputExamples", "pathwayLayout_FBAcolored", "iND750_Output",
                               "iND750_output_whole.tlp")
test_yeast_graph = cbmviz.loadGraphObjFromFile(test_yeast_path)
test_yeast_graph.cbmviz_file_name = "iND750"
test_yeast_list = [test_yeast_graph, ["Citric Acid Cycle", "Glycolysis/Gluconeogenesis"]]
# append the different dicts to the input list
test_input_list = list()
test_input_list.append(test_pylori_list)
test_input_list.append(test_ecoli_list)
test_input_list.append(test_yeast_list)
# generate the combined graph
output_combine_graph = cbmviz.combineSubgraphsFromGraphs(test_input_list)
print("output_combine_graph", output_combine_graph)
# save combined graph in new output folder
cbmviz.saveGraphObjToTulipFile(output_combine_graph, "test_combine_output.tlp")


# TODO: put the graphs in one (or more?) graph dicts and run coloring+exporting functions
#cbmviz.colorGraphs(global_ops, fba_results, graph_registry, output_path, cmod.getId())
#cbmviz.exportGraph(output_path, graph_registry, global_ops, cDir, cmod.getId(), groups_sbml, n2lo, e2lo, ub_sps)


