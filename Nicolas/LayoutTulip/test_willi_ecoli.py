#from MetaViz import *
import os, time
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
import MetaViz

global_ops = {"AUTO_ALIAS": False,  # apply aliasing to every node of the graph with a certain connectivity
              # and option to alias biomass (result in maingraph)
              # if auto alias is active the automatically generated layout will still be
              # overwritten if any other layouting functions are active aswell
              "AUTO_ALIAS_OPTIONS": {"spec_split": 6, "reac_split": 4},  # threshold values for species and reac
              "SMALL_THRSH": 30,  # the threshold value of connections a species must have to be a small sps
              # splitting

              "PATHWAY_LISTS": True,  # layout annotated groups of sbml file as seperate pathways
              "LAYOUT_PATHWAYS": True,  # needs PATHWAY_LISTS true, so that necessary node lists are generated
              "LABEL_TYPE": "id",  # can either be "id" or "name"
              "ADD_COFACTORS": ["Glycolysis-Gluconeogenesis", 'Citric Acid Cycle'],  # add nodes for small sps/cofactors after they have been
              # deleted in path_dict

              "TRANSPORT_RADIAL": False,  # apply Tree Radial algorithm to transport reactions in output_whole
              "SELECTED_PATHWAYS": ["Glycolysis-Gluconeogenesis",  "Citric Acid Cycle", "Transport, Outer Membrane", 'Transport, Outer Membrane Porin'],
              "GRID_OUTPUT": ["Glycolysis-Gluconeogenesis",  # align pathways on a grid, control position
                              "Citric Acid Cycle",
                              "Alternate Carbon Metabolism",
                              "Glyoxylate Metabolism"],

              "EXPORT_TLP": True,  # save generated graphs as .tlp files
              "EXPORT_GML": False,  # export generated graphs as .gml files
              "EXPORT_SVG": True,  # export generated graphs as .svg files
              "EXPORT_JSON": False,  # export generated graphs as .json files
              "EXPORT_SBML": False,  # export generated graphs as .sbml files (in layout plugin)
              "EXPORT_SBGN": False,  # export generated graphs as .sbgn files

              "COLOR_FBA": True,  # color the reaction nodes and in and outgoing edges
              "COLOR_FBA_COLORBLIND": False,  # based on flux values obtained from FBA
              "COLOR_PATHWAYS": False,  # pathways specific coloring
              "COLOR_GAPS": True  # highlight gaps (in red), reactions that exist in both (all) models = black
              }

#sbml = os.path.join(cDir, "TestModels", "bigg2", "iAF1260.xml") # E. coli about 2000+ reactions
sbml = os.path.join(cDir, "TestModels", "willi", "iAF1260_max_ac.xml") # E. coli about 2000+ reactions

TIME0 = time.time()
graph, graph_props = MetaViz.createInitialGraph("initialGraph")
TIME1 = time.time()
cmod, fba_results, path_dict, COLORS, ub_sps, n2lo, e2lo, groups_sbml, obj_re_nodes, graph_registry, cDir, output_path =\
                                                                MetaViz.layoutModel(global_ops, sbml, graph, graph_props)
TIME2 = time.time()
MetaViz.colorGraphs(COLORS, global_ops, fba_results, path_dict, graph_registry, output_path, cmod.getId())
TIME3 = time.time()
MetaViz.exportGraph(output_path, graph_registry, global_ops, cDir, cmod.getId(), groups_sbml, n2lo, e2lo, ub_sps)
TIME4 = time.time()
#db_path = MetaViz.buildSQLiteDatabase(graph_registry['output_whole'][0], output_path)
#db_graph_node_list, node_col_headers = MetaViz.retrieveDataFromDatabaseTable(db_path, 'graph_nodes')
#db_graph_edge_list, edge_col_headers = MetaViz.retrieveDataFromDatabaseTable(db_path, 'graph_edges')
#db_graph_edgeconn_list, edgeconn_col_headers = MetaViz.retrieveDataFromDatabaseTable(db_path, 'edge_connectivity')
#db_graph = MetaViz.buildTulipGraphFromSQL(db_graph_node_list, db_graph_edge_list, db_graph_edgeconn_list,
                                          #node_col_headers, edge_col_headers, "graph_fromDB", output_path)


print('\nTotal time:    {}'.format(round(TIME4-TIME0, 2)))
print('- create time: {}'.format(round(TIME1-TIME0, 2)))
print('- layout time: {}'.format(round(TIME2-TIME1, 2)))
print('- colour time: {}'.format(round(TIME3-TIME2, 2)))
print('- export time: {}'.format(round(TIME4-TIME3, 2)))


#python -m cProfile [-o output_file] [-s sort_order] myscript.py
#import pstats
#p = pstats.Stats('restats')
#p.strip_dirs().sort_stats(-1).print_stats()

"""
In [17]: P = p.strip_dirs().sort_stats('cumulative').print_stats(50)
Thu Feb 09 16:13:11 2017    xlayout.prf

         116165942 function calls (115595154 primitive calls) in 91.128 seconds

   Ordered by: cumulative time
   List reduced from 9544 to 50 due to restriction <50>

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.058    0.058   91.140   91.140 test_willi_ecoli.py:2(<module>)
        1    0.003    0.003   86.404   86.404 MetaViz.py:40(layoutModel)
        1    0.002    0.002   45.786   45.786 tulip_graph2geojson.py:345(graph2geojson)
        1    0.005    0.005   33.791   33.791 tulip_graph2geojson.py:266(meta_graph2features)
        4    0.000    0.000   31.957    7.989 tulip_graph2geojson.py:258(export_elements)
       10    0.493    0.049   31.504    3.150 tulip_graph2geojson.py:35(export_edges)
    27754    1.254    0.000   30.882    0.001 tlp2geojson.py:65(e2feature)
      391   17.624    0.045   17.624    0.045 {built-in method applyLayoutAlgorithm}
        1   11.410   11.410   12.644   12.644 MetaViz.py:600(mapModelToGraph)
24462/11112    0.219    0.000   12.480    0.001 function.py:362(__new__)
    69359    0.387    0.000   12.266    0.000 evalf.py:1329(evalf)
    60365    0.115    0.000   12.228    0.000 expr.py:217(__float__)
24461/11111    0.198    0.000   11.504    0.001 function.py:194(__new__)
334565/69359    1.018    0.000   11.487    0.000 evalf.py:1274(evalf)
        1    0.001    0.001   11.360   11.360 tulip_graph2geojson.py:441(process_compartments)
        6    0.005    0.001   11.309    1.885 ubiquitous_layout.py:379(layout)
        6    0.004    0.001   10.971    1.829 layout_utils.py:81(layout_components)
    11111    0.113    0.000   10.970    0.001 trigonometric.py:2619(eval)
   261138    0.561    0.000   10.757    0.000 assumptions.py:239(getit)
    56326    0.418    0.000   10.462    0.000 evalf.py:475(evalf_add)
        1    0.525    0.525   10.186   10.186 MVLayout.py:311(reactionsToNodes)
        5    0.000    0.000    9.883    1.977 layout_utils.py:56(layout_force)
    55508    0.206    0.000    9.293    0.000 tlp2geojson.py:38(get_border_coord)
        1    0.019    0.019    8.274    8.274 MetaViz.py:701(layoutAllPathwaysFromGroups)
       40    1.122    0.028    8.255    0.206 MVLayout.py:677(layoutPathway)
   490051    4.174    0.000    8.086    0.000 facts.py:503(deduce_all_facts)
13335/9020    0.299    0.000    7.694    0.001 trigonometric.py:2086(eval)
       41    0.001    0.000    7.047    0.172 MVLayout.py:642(layoutSubgraph)
   241685    0.200    0.000    7.018    0.000 assumptions.py:223(copy)
   242056    0.492    0.000    6.820    0.000 assumptions.py:212(__init__)
    27693    4.401    0.000    6.035    0.000 CBModel.py:1693(getSpecies)
        5    0.000    0.000    5.607    1.121 layout_utils.py:71(remove_overlaps)
   112712    0.635    0.000    5.059    0.000 evalf.py:392(add_terms)
367686/279323    0.345    0.000    4.925    0.000 decorators.py:70(__sympifyit_wrapper)
    35646    0.146    0.000    4.813    0.000 expr.py:1881(extract_multiplicatively)
        6    0.076    0.013    4.704    0.784 factoring.py:84(comp_to_meta_node)
    35878    0.584    0.000    4.650    0.000 evalf.py:517(evalf_mul)
"""
