import os
import cbmviz


# User Input
global_ops = cbmviz.cfg_default
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
sbml = os.path.join(cDir, "TestModels", "bigg2", "iIT341.xml") # H.pylori

# run CBMViz core pipeline to create a layout
graph, graph_props = cbmviz.createInitialGraph("initialGraph")
cmod, fba_results, graph_registry, cDir, output_path = cbmviz.layoutModel(global_ops, sbml, graph, graph_props)
cbmviz.colorGraphs(global_ops, fba_results, graph_registry, output_path, cmod.getId())
cbmviz.exportGraph(output_path, graph_registry, global_ops, cDir, cmod.getId())

# test: SQL operations
# create SQLite db to store a graph object
db_path = cbmviz.buildSQLiteDatabase(graph_registry['output_whole'][0], output_path)
# retrieve graph data from SQL database
db_graph_node_list, node_col_headers = cbmviz.retrieveDataFromDatabaseTable(db_path, 'graph_nodes')
db_graph_edge_list, edge_col_headers = cbmviz.retrieveDataFromDatabaseTable(db_path, 'graph_edges')
db_graph_edgeconn_list, edgeconn_col_headers = cbmviz.retrieveDataFromDatabaseTable(db_path, 'edge_connectivity')
# SQL --> back to CBMViz
db_graph = cbmviz.buildTulipGraphFromSQL(db_graph_node_list, db_graph_edge_list, db_graph_edgeconn_list,
                                  node_col_headers, edge_col_headers, "graph_fromDB", output_path)
# export database graph
db_graph_dict = cbmviz.createGraphDictFromGraphObj(db_graph)
cbmviz.exportGraph(output_path, db_graph_dict, global_ops, cDir, cmod.getId())
