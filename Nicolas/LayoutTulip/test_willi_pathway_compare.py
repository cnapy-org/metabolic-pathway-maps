#from MetaViz import *
import os
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
import MetaViz

global_ops = {"AUTO_ALIAS": False,  # apply aliasing to every node of the graph with a certain connectivity
              # and option to alias biomass (result in maingraph)
              # if auto alias is active the automatically generated layout will still be
              # overwritten if any other layouting functions are active aswell

              "PATHWAY_LISTS": True,  # layout annotated groups of sbml file as seperate pathways
              "LAYOUT_PATHWAYS": True, # needs PATHWAY_LISTS true, so that necessary node lists are generated

              "TRANSPORT_RADIAL": True,  # apply Tree Radial algorithm to transport reactions in output_whole
              "SELECTED_PATHWAYS": ["objective reaction", "Transport"],  # pathways to be saved in output_selected
              #"SELECTED_PATHWAYS": False,

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

sbml = os.path.join(cDir, "TestModels", "bigg2", "iIT341.xml")
graph, graph_props = MetaViz.createInitialGraph("initialGraph")
cmod, fba_results, path_dict, COLORS, ub_sps, n2lo, e2lo, groups_sbml, obj_re_nodes, graph_registry = \
    MetaViz.layoutModel(global_ops, sbml, graph, graph_props)

wmod = cmod.clone()

new_sub_model = wmod.getReactionValues()
cntr = 0
for x in new_sub_model.keys():
    if cntr <= 200:
        new_sub_model.pop(x)
    else:
        new_sub_model[x] = 1000.0
    cntr += 1

MetaViz.colorGraphs(global_ops, new_sub_model, path_dict, graph_registry)

new_sub_model2 = wmod.getReactionValues()
cntr = 0
for x in new_sub_model2.keys():
    if cntr > 200:
        new_sub_model2.pop(x)
    else:
        new_sub_model2[x] = -1000.0
    cntr += 1


MetaViz.colorGraphs(global_ops, new_sub_model2, path_dict, graph_registry)
print("Graph(s) colored. \n")
MetaViz.exportGraph(graph_registry, global_ops, cDir, cmod.getId(), groups_sbml, n2lo, e2lo, ub_sps)
