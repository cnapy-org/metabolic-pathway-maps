#from MetaViz import *
import os
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
import MetaViz

global_ops = {"AUTO_ALIAS": False,  # apply aliasing to every node of the graph with a certain connectivity
              # and option to alias biomass (result in maingraph)
              # if auto alias is active the automatically generated layout will still be
              # overwritten if any other layouting functions are active aswell
              "AUTO_ALIAS_OPTIONS": {"spec_split": 6, "reac_split": 4},  # threshold values for species and reac
                                                                          # splitting

              "PATHWAY_LISTS": True,  # layout annotated groups of sbml file as seperate pathways
              "SMALL_THRSH": 30,  # the threshold value of connections a species must have to be a small sps
              "LAYOUT_PATHWAYS": True,  # needs PATHWAY_LISTS true, so that necessary node lists are generated
              "LABEL_TYPE": "id",  # can either be "id" or "name"
              "ADD_COFACTORS": ["Respiratory chain",    # add nodes for small sps/cofactors after they have been
                                ],  # deleted in path_dict

              "TRANSPORT_RADIAL": True,  # apply Tree Radial algorithm to transport reactions in output_whole
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

              "COLOR_FBA": False,  # color the reaction nodes and in and outgoing edges
              "COLOR_FBA_COLORBLIND": False,  # based on flux values obtained from FBA
              "COLOR_PATHWAYS": False,  # pathways specific coloring
              "COLOR_GAPS": True  # highlight gaps (in red), reactions that exist in both (all) models = black
              }
sbml = os.path.join(cDir, "TestModels", "bigg2", "iIT341.xml")

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


wmod = cmod.clone()

new_sub_model = wmod.getReactionValues()
cntr = 0
for x in new_sub_model.keys():
    if cntr <= 200:
        new_sub_model.pop(x)
    else:
        new_sub_model[x] = 1000.0
    cntr += 1

# 3. apply coloring algorithms to graphs
MetaViz.colorGraphs(COLORS, global_ops, new_sub_model, path_dict, graph_registry, output_path, cmod.getId())
print("Graph(s) colored. \n")
# 4. export graph SVG, GML, JSON, SBML+layout Plugin and/or SBGN file (controlled by global_ops)
MetaViz.exportGraph(output_path, graph_registry, global_ops, cDir, cmod.getId(), groups_sbml, n2lo, e2lo, ub_sps)