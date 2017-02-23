import os
import cbmviz
import networkx as nx

# test SQL -> NetworkX -> export as DOT file (graphviz format) -> load to get graph object back!
# user input
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
file_name = "iIT341"
output_path = os.path.join(cDir, '{}_Output'.format(file_name))
db_path = os.path.join(cDir, "OutputExamples", "pathwayLayout_FBAcolored", "iIT341_Output", "iIT341_graph.db")

# access SQL db
db_graph_node_list, node_col_headers = cbmviz.retrieveDataFromDatabaseTable(db_path, 'graph_nodes')
db_graph_edgeconn_list, edgeconn_col_headers = cbmviz.retrieveDataFromDatabaseTable(db_path, 'edge_connectivity')
# select only Glycolysis nodes
db_glyc_node_list = list()
for sublist in db_graph_node_list:
    if sublist[node_col_headers.index("pathway")] == "Glycolysis":
        db_glyc_node_list.append(sublist[node_col_headers.index("node_id")])

#use networkx to rebuild only Glycolysis
G = nx.Graph()
# build edges in networkx graph G between the nodes by looking up the connectivity per node id
for sublist in db_graph_edgeconn_list:
    for node_id in db_glyc_node_list:
        if str(node_id) in sublist[edgeconn_col_headers.index("source_node_id")]:
            G.add_edge(node_id, sublist[edgeconn_col_headers.index("target_node_id")])
        elif str(node_id) in sublist[edgeconn_col_headers.index("target_node_id")]:
            G.add_edge(sublist[edgeconn_col_headers.index("source_node_id")], node_id)

# get plugin to write dot files (copied from networkx documentation)
try:
    import pygraphviz
    from networkx.drawing.nx_agraph import write_dot
    print("using package pygraphviz")
except ImportError:
    try:
        import pydotplus
        from networkx.drawing.nx_pydot import write_dot
        print("using package pydotplus")
    except ImportError:
        print()
        print("Both pygraphviz and pydotplus were not found ")
        print("see http://networkx.github.io/documentation"
              "/latest/reference/drawing.html for info")
        print()
        raise
# write networkx graph as dot file
write_dot(G, "glycolysis.dot")
# read dot file with cbmviz
test_glyc_dot_file_path = os.path.join(cDir, "glycolysis.dot")
test_glyc_dot_graph = cbmviz.loadGraphObjFromFile(test_glyc_dot_file_path)
# build graph_registry dict from new dot test graph
test_glyc_graph_dict = cbmviz.createGraphDictFromGraphObj(test_glyc_dot_graph)
# export test dot graph as new tlp file
cbmviz.exportGraph(output_path, test_glyc_graph_dict, cbmviz.cfg_default, cDir, file_name)