import logging
from tulip import tlp
import cbmpy

from .mod_sbml.annotation.gene_ontology.go_annotator import annotate_compartments, get_go_id
from .sbml_helper import parse_group_sbml, GrPlError
from .sbml_manager import check_names, check_compartments
from .mod_sbml.sbml.compartment.compartment_positioner import comp2level
from .mod_sbml.onto import parse_simple
from .mod_sbml.onto.term import Term
from .mod_sbml.annotation.chebi.chebi_serializer import get_chebi
from .mod_sbml.annotation.gene_ontology.go_serializer import get_go
from .mod_sbml.sbml.sbml_manager import get_gene_association, get_formulas
from .mod_sbml.annotation.chebi.chebi_annotator import get_chebi_id
from .node_cloner import clone_node
from .graph_properties import *
from .resize import get_n_size

from . import MVLayout as mvl

__author__ = 'anna'

SKIP_UBS = False


def species2nodes(graph, input_model, ub_sps):
    id2n = {}
    for s in input_model.getListOfSpecies():
        _id = s.getId()
        ub = (_id in ub_sps)

        # todo
        if SKIP_UBS and ub:
            continue

        n = graph.addNode()
        comp = input_model.getCompartment(s.getCompartment())
        graph[COMPARTMENT_ID][n] = comp.getId()
        graph[ID][n] = _id
        id2n[_id] = n

        chebi_id = get_chebi_id(s)
        if chebi_id:
            graph[TERM][n] = get_chebi_id(s)
        name = s.getName()
        if name:
            # remove compartment from the name,
            # e.g. H2O [peroxisome] --> H2O
            if comp.getName():
                name = name.replace('[%s]' % comp.getName(), '').strip()
                name = name.replace('[%s]' % comp.getId(), '').strip()
            graph[NAME][n] = name
        graph[TYPE][n] = TYPE_SPECIES
        graph[UBIQUITOUS][n] = ub
        graph[VIEW_SHAPE][n] = SPECIES_SHAPE
        graph[VIEW_SIZE][n] = get_n_size(graph, n)
        formulas = get_formulas(s)
        if formulas:
            formula = next(iter(formulas))
            if formula and formula != '.':
                graph[FORMULA][n] = formula
    return id2n


def get_quotient_maps(chebi, sbml_file):
    try:
        r_id2g_id, s_id2gr_id, ub_sps = parse_group_sbml(sbml_file, chebi)
        return r_id2g_id, s_id2gr_id, ub_sps
    except GrPlError as e:
        logging.error("Could not parse the groups SBML %s: %s" % (sbml_file, e.msg))
        return None, None, None


def compute_c_id2info(c_id2level, input_model):
    c_id2info = {}
    for comp in input_model.getListOfCompartments():
        c_id = comp.getId()
        c_name = comp.getName()
        if not c_name:
            c_name = c_id
        go_id = get_go_id(comp)
        c_id2info[c_id] = (c_name, go_id if go_id else '', c_id2level[c_id])

    c_id2outs = {}
    for c_id in c_id2info.keys():
        _, _, (_, out_c_id) = c_id2info[c_id]
        outs = []
        while out_c_id:
            outs.append(out_c_id)
            _, _, (_, out_c_id) = c_id2info[out_c_id]
        c_id2outs[c_id] = outs

    return c_id2info, c_id2outs

# COMPLETELY CHANGED FUNCTION
# former function name: import_sbml
# -> all other functionality has been moved to TestSBML2TulipGraph.py
def cloneNodes(ubs, graph):
   duplicate_nodes(graph)
   return graph


def duplicate_nodes(graph):
    root = graph.getRoot()
    duplicate_ub_nodes(graph, root)
    # duplicate_tr_nodes(graph, root)


def duplicate_tr_nodes(graph, root):
    for n in (n for n in graph.getNodes() if TYPE_SPECIES == root[TYPE][n] and not root.isMetaNode(n) and
            not root[UBIQUITOUS][n] and next((r for r in root.getInOutNodes(n) if root[TRANSPORT][r]), False)):
        clone_node(graph, n, get_neighbour_nodes=lambda m, gr: (t for t in gr.getInOutNodes(m) if gr[TRANSPORT][t]))


def duplicate_ub_nodes(graph, root):
    for n in (n for n in graph.getNodes() if root[UBIQUITOUS][n]):
        clone_node(graph, n)


def mark_ancestors(graph, r_eq2clu, s2clu, c_id2info):
    root = graph.getRoot()
    id_ = root.getStringProperty(ID)
    anc_id = root.getStringProperty(ANCESTOR_ID)
    anc_name = root.getStringProperty(ANCESTOR_NAME)
    anc_ch_id = root.getStringProperty(ANCESTOR_TERM)
    type_ = root.getIntegerProperty(TYPE)
    for n in graph.getNodes():
        gr_id, gr_name, term = None, None, None
        if TYPE_REACTION == type_[n]:
            if id_[n] in r_eq2clu:
                gr_id, gr_name, el_num = r_eq2clu[id_[n]]
                gr_name = "%s (%d)" % (gr_name, el_num)
        elif id_[n] in s2clu:
            gr_id, term, el_num = s2clu[id_[n]]
            if term:
                if isinstance(term, Term):
                    gr_name = term.get_name()
                else:
                    gr_name = term
                    term = None
                c_name, _, _ = c_id2info[root[COMPARTMENT_ID][n]]
                gr_name = "%s (%d) [%s]" % (gr_name, el_num, c_name)
        if gr_name:
            anc_name[n] = gr_name
        if gr_id:
            anc_id[n] = gr_id
        if term:
            anc_ch_id[n] = term.get_id()

