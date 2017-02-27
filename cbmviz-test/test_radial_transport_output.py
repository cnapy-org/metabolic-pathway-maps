import os
import cbmviz

global_ops = cbmviz.cfg_default


global_ops.update({"METABOLITE_NODE_SIZE": 15})
global_ops.update({"COFACTOR_NODE_SIZE": 12.5})
global_ops.update({"REACTION_NODE_SIZE": 10})

# test radial transport output
global_ops.update({"TRANSPORT_RADIAL": True})
global_ops.update({"ADD_COFACTORS": ["TCA cycle"]})

sbml = os.path.join("TestModels", "bigg2", "iIT341.xml")
graph, graph_props = cbmviz.createInitialGraph("iIT341")
cmod, fba_results, graph_registry, cDir, output_path = cbmviz.layoutModel(global_ops, sbml, graph, graph_props)
cbmviz.colorGraphs(global_ops, fba_results, graph_registry, output_path, "iIT341")
cbmviz.exportGraph(output_path, graph_registry, global_ops, cDir, "iIT341")