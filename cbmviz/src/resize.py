from tulip import *
from .graph_properties import *

# default node sizes
SPECIES_SIZE = 2.5
UBIQUITOUS_SPECIES_SIZE = 2
REACTION_SIZE = 1

# default edge sizes
UBIQUITOUS_EDGE_SIZE = 1
TRANSPORT_EDGE_SIZE = 1
EDGE_SIZE = 1.5
# GENERALIZED_EDGE_SIZE = 2
# COMP_EDGE_SIZE = 2.5


def get_n_length(graph, n):
    return graph.getRoot()[VIEW_META_GRAPH][n].numberOfNodes() if graph.isMetaNode(n) else 1


def _get_e_size(root, e, global_ops):
    if root[UBIQUITOUS][e]:
        return UBIQUITOUS_EDGE_SIZE
    if not root.isMetaEdge(e) or not root.isMetaNode(root.source(e)) or not root.isMetaNode(root.target(e)):
        return TRANSPORT_EDGE_SIZE if root[TRANSPORT][root.source(e)] or root[TRANSPORT][root.target(e)] else \
            EDGE_SIZE
    return sum(_get_e_size(root, edg, global_ops) for edg in root[VIEW_META_GRAPH][e])


def get_mn_size(n, root):
    bb = tlp.computeBoundingBox(root[VIEW_META_GRAPH][n])
    return bb.width(), bb.height()


def get_n_size(graph, n, global_ops):
    root = graph.getRoot()
    n_type = root[TYPE][n]
    # take global_ops input for node size
    SPECIES_SIZE = global_ops["METABOLITE_NODE_SIZE"]  # default 2.5
    UBIQUITOUS_SPECIES_SIZE = 2
    REACTION_SIZE = global_ops["REACTION_NODE_SIZE"]  # default 1
    if TYPE_REACTION == n_type:
        if root[FAKE][n]:
            w = h = max(get_mn_size(n, root))
        else:
            w = h = REACTION_SIZE * get_n_length(graph, n)
    elif TYPE_COMPARTMENT == n_type:
        w, h = get_mn_size(n, root)
    elif root[UBIQUITOUS][n]:
        w = h = UBIQUITOUS_SPECIES_SIZE
    else:
        w = h = SPECIES_SIZE * get_n_length(graph, n)
    return tlp.Size(w, h)


def get_e_size(graph, e, global_ops):
    sz = _get_e_size(graph.getRoot(), e, global_ops)
    return tlp.Size(sz, sz)
