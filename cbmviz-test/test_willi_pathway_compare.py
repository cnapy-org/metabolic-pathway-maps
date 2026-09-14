import os
import cbmviz

# User Input
global_ops = cbmviz.cfg_default
global_ops.update({"COLOR_FBA": False})
global_ops.update({"COLOR_GAPS": True})
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
sbml = os.path.join(cDir, "TestModels", "bigg2", "iIT341.xml")

# MetaViz Pipeline
graph, graph_props = cbmviz.createInitialGraph("initialGraph")
cmod, fba_results, graph_registry, cDir, output_path = cbmviz.layoutModel(global_ops, sbml, graph, graph_props)


wmod = cmod.clone()
# new_sub_model is a dictionary with reaction ids (keys) and flux values (values)
# flux values get either deleted or overwritten by the value 1000 (to have big contrast for binary coloring)
new_sub_model = wmod.getReactionValues()
cntr = 0
for x in list(new_sub_model.keys()):
    if cntr <= 200:
        new_sub_model.pop(x)
    else:
        new_sub_model[x] = 1000.0
    cntr += 1


cbmviz.colorGraphs(global_ops, new_sub_model, graph_registry, output_path, cmod.getId())
cbmviz.exportGraph(output_path, graph_registry, global_ops, cDir, cmod.getId())