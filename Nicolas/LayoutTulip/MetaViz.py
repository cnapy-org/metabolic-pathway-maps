########################################################################################################################
#  MetaViz Prototype - Visualization Pipeline for Metabolic Networks                                                   #
#  Nicolas Huber, Dr. Brett Olivier, Prof. Bas Teusink - VU Amsterdam, 2016 - 2017                                     #
#  Generate a structured network map using only a SBML file (idealy with defined groups) as input. Resulting graph is  #
#  initially created as Tulip graph object and can be exported as SVG, JSON, GML, new SBML with layout plugin as well  #
#  as SBGN. Graph properties are also stored in a SQLite database.                                                     #
#  Necessary  dependencies for this pipeline (on Wind64 Anaconda):                                                     #
#  pip install cbmpy geojson tarjan libsbgnpy tulip-python tulipgui-python                                             #
#  SQLite module: sqlite 3 (https://www.tutorialspoint.com/sqlite/sqlite_installation.htm)                             #
########################################################################################################################


import os
# tulip imports
import sys
sys.path.append('C:\\Program Files\\Tulip-4.10.0\\bin')
from tulipgui import * # for SVG export
from tulip import tlp
# modules for SBML model interaction
import libsbml
import cbmpy
# SQL database API
import sqlite3
# type conversion for db storage
import pickle
# necessary Mimoza modules for graph visualization
# TODO: put all the Mimoza in a mimoza_functions folder with __ini__.py
import sbml2tlp
import tulip_graph2geojson
# sbml and sbgn export
import sbml_manager
import sbgn_helper
# MetaViz functions
import MVLayout as mvl


########################################################################################################################
# wrapping functions to execute entire blocks of logically connected function calls
# Part 1: Layout
def layoutModel(global_ops, sbml, graph, graph_props):
    """
    MetaViz Part 1: Layout
    ----------------------

    Take an sbml file, a Tulip graph and its prop dict and return:
    - the CBMPY model object (cmod)
    - a dict with fba results - reaction ids and fluxes (fba_results)
    - a dict with colors and their RGB values (COLORS)
    - a list of ubiquitous species identified by autoAlias (ub_sps)
    - two dict with node and edge layout info of initial graph (n2lo and e2lo)
    - a path to the modified sbml file - with group for ub_sps (groups_sbml)
    - a list of nodes associated with the objective reaction (obj_re_nodes)
    - a dict with all created graph objects (graph_registry)
    - the current directory (cDir)
    - the path for all output files (output_path)::

    - *sbml* a path to an sbml model file
    """
    # get current directory
    cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))
    # create CBMPY model (separate data structure to SBML file)
    cmod, stM, objective_re_id = createCBMPYModel(sbml)
    # create output folder
    output_path = os.path.join(cDir, '{}_Output'.format(cmod.getId()))
    try:
        os.makedirs(output_path)
    except OSError:
        if not os.path.isdir(output_path):
            raise
    # do FBA and save results
    fba_results = runFBA(cmod)
    #################################################
    # select nodes for aliasing (different methods) #
    #################################################
    # 1.) create alias node list based on connectivity <-> user defined thresholds (aliasing by Mimoza)
    ub_list, ub_sps, groups_sbml = autoAlias(cmod, global_ops, output_path,
                                             global_ops["AUTO_ALIAS_OPTIONS"]["spec_split"],
                                             global_ops["AUTO_ALIAS_OPTIONS"]["reac_split"])
    # 2.) create alias node list based on biological pathways <-> group info in sbml - if available (aliasing done by own
    # functions)
    # IF GROUPS DON'T EXIST -> provide pathway lists manually:
    # use function manualPathDict
    # IF GROUPS EXIST -> fetch groups from sbml (this doesn't use ub_sps_list or ub_re_list)
    # create dictionary of groups where each group name is a key and the values are again dicts with three keys: 'id',
    # 'reactions' and 'metabolites'
    if global_ops["PATHWAY_LISTS"] is True and global_ops["AUTO_ALIAS"] is False:
        path_dict, small_sps_dict = automPathDict(cmod, global_ops, global_ops["SMALL_THRSH"])
    else:
        # if model has no groups defined -> return empty dicts
        path_dict = dict()
        small_sps_dict = dict()

    # flag (true/false) a list of nodes (already saved as ub_list) so that they can be aliased later
    mvl.markUbiquitousNodes(graph, ub_list)
    # create initial nodes for all species
    id2n = mvl.speciesToNodes(graph, cmod, ub_list, graph_props)
    # save compartment hierarchy to dictionary "c_id2level" for function import_sbml
    # required syntax: c_id2level = {'Cell': (0, ''), 'Boundary': (1, 'Cell')}
    c_id2level = compDictForMimoza(cmod)
    c_id2info, c_id2outs = compute_c_id2info(c_id2level, cmod)
    # create nodes for all reactions
    mvl.reactionsToNodes(get_r_comp, graph, id2n, cmod, c_id2level, c_id2outs)
    # clone/alias nodes selected above (in ub_list)
    # -> chebi info was apparently only important for automatically extracting ubiquitous species ids! (we now provide them)
    # => provide them by passing a list (ubs), converted to a set (cloneNodes needs unique ids)
    root = sbml2tlp.cloneNodes(set(ub_list), graph)
    ##############################
    # run Mimoza layout functions #
    ##############################
    # actual layouting, creation of graph structure
    # n2lo and e2lo are dictionaries containing node and edge layout info (id, color, coords, etc.) respectively
    fc, (n2lo, e2lo), hidden_c_ids, c_id_hidden_ubs, graph = tulip_graph2geojson.graph2geojson(c_id2info, c_id2outs,
                                                                                               root)
    # initialize a dict of all existing graphs (key: name, value: a list [graph object, graph props])
    #  -> everytime a new, important graph is created it gets added to the graph registry
    graph_registry = {}
    # add initial graph to graph registry (before layouting)
    graph.setName('initial graph')
    graph_registry['initial graph'] = [graph, graph_props]
    # graph_registry.update({"initial graph": [graph, graph_props]})
    # save raw Mimoza graph output as backup
    initial_path = os.path.join(output_path, "{}_initial.tlp".format(cmod.getId()))
    tlp.saveGraph(graph, initial_path)
    print("Initial Tulip graph saved as: {}_initial.tlp. \n".format(cmod.getId()))
    ###############################################################
    # map sbml element ids to initial graph nodes in a dictionary #
    ###############################################################
    # Section Description:
    # --------------------
    # in the initial graph every model id is represented by exactly one node (if AUTO_ALIAS is off!)
    # example: iIT341 initial graph (returned by Mimoza): 1039 nodes created
    #          check: len(cmod.getSpeciesIds()+cmod.getReactionIds()) returns 1039 as well
    #
    # => create mapping between ids and their nodes in a dict that reflects the state of the raw graph
    #    dict structure: id (key): value is another dict with layout (key): True/False and initial_node (key): node object
    #    layout key serves as a flag to see whether an initial node has already been used for layouting later (True = has
    #    been used)
    #    -> if the downstream layouting function sees layout True => creates new node and copies properties
    #                                                                (new node will not be added to this dict)
    # use only if one model id = one graph node is true! (only if auto aliasing is off and  the current graph oject is
    # in its raw, initial state)
    # subsequent functions either get an empty version or the actual graph_elements dict
    graph_elements = {}
    if global_ops["AUTO_ALIAS"] is False:
        graph_elements = mapModelToGraph(cmod, graph, graph_props)
    obj_re_nodes = objReNodeList(graph, graph_props, objective_re_id)
    # Why is a mapping between ids and graph elements (nodes, edges) necessary?
    # -------------------------------------------------------------------------
    #
    # There is a difference between a model species id and a graph node id!
    # => It is not sufficient to select nodes based on model ids (at some point multiple nodes with the same id can exist)
    # => handle node objects!
    #
    # Problem: if a node is shared between pathways, the last valid layouting call will determine its position meaning that
    # previous positions of that node in other layouted subgraphs that also contain it are overwritten!
    #
    # example in Carlsson model: since the initial nodes with the ids "PEP" and "R1" are part of glycolysis AND exchange
    # reactions, the assigned positions of the "PEP" and "R1" initial nodes after layouting glyc are overwritten by
    # layouting exch re!
    #
    # -> give subgraphs their own versions of a shared node by tracking what graph elements have already been used for
    #    layouting to determine what nodes and edges need to be recreated
    #################################
    # change basic graph properties #
    #################################
    # TODO: what happens to the output if labels are names (which are longer than ids)
    # TODO: in global_ops: control type of label
    # TODO: also differentiate between metabolite labels and reaction labels
    labelGraphNodes(graph, global_ops["LABEL_TYPE"])
    # TODO: add labels to edges (possible? alternative?)
    # for e in graph.getEdges():
    #    viewLabel[e] = "test"
    # -> edge labels are set in "viewLabel" but not visible!
    COLORS = colorAllGraphElements(graph, "egdes", "grey")
    # change line types (Polyline, BezierCurve, CatmullRomCurve, CubicSplineCurve)
    # for e in graph.getEdges():
    #    graph_props["viewShape"][e] = tlp.EdgeShape.BezierCurve
    # TODO: change font globally from 18 -> 8 and 24 -> 12 or something
    # for n in graph.getNodes():
    #    graph_props["font_size"] = 8
    # TODO: also write a function that applies "fit lable to node size" thing from tulipgui
    # TODO: see if calling that before layouting changes the layout
    ##############################################################
    # layout pathways in subgraphs (Glycolysis, TCA Cycle, etc.) #
    ##############################################################
    # needs dict graph elements which is only created it global_ops["AUTO_ALIAS"] is false since it needs 1:1 model id
    # to graph elements
    if global_ops["LAYOUT_PATHWAYS"] is True and global_ops["AUTO_ALIAS"] is False:
        layoutAllPathwaysFromGroups(cmod, global_ops, path_dict, graph, graph_elements, graph_props, small_sps_dict)
        # create and layout the objective reaction block (sg gets added to 'graph')
        mvl.layoutObjFunc(graph, objective_re_id, "Sugiyama (OGDF)", graph_props)
        pathwayDump(graph, global_ops, graph_registry, graph_props, global_ops["SELECTED_PATHWAYS"], global_ops["GRID_OUTPUT"])

    return cmod, fba_results, path_dict, COLORS, ub_sps, n2lo, e2lo, groups_sbml, obj_re_nodes, graph_registry, cDir,\
           output_path
# Part 2: Coloring
def colorGraphs(COLORS, global_ops, fba_results, path_dict, graph_registry, output_path, model_id):
    """
    MetaViz Part 2: Coloring
    ------------------------

    Take a dict of graphs (graph_registry) and return:
    - ::

    - *COLORS* a dict with predefined colors (RGB values)
    - *global_ops* a dict with global control variables
    - *fba_results* a dict with flux values for all reactions
    - *path_dict* a dict with pathway ids and associated metabolites and reactions
    - *graph_registry* a dict of Tulip graphs and their properties
    - *output_path* an os.path for all output files
    - *model_id* the id of the sbml model

    """
    #######################################
    # COLORING: visualize FBA result, ... #
    #######################################
    # color nodes and edges depending on flux value
    # input format: dict = {"M1": -3, ...}
    # TODO: create new dict: key -> round(flux value, 1), value -> color of the associated node with that flux value
    if global_ops["COLOR_FBA"]:
        for graph_data in graph_registry.values():
            graph_obj = graph_data[0]
            mvl.colourFBA(graph_obj, fba_results, COLORS["light grey"])
        print("Applied flux coloring to all graphs. \n")
    if global_ops["COLOR_GAPS"]:
        # highlight gaps in red, reactions with 0 value (exist in both) -> black
        for graph_data in graph_registry.values():
            graph_obj = graph_data[0]
            mvl.colourFBA(graph_obj, fba_results, COLORS["red"], binary_scale=[tlp.Color.Black, tlp.Color.Black])
        print("Applied gap coloring to all graphs. \n")
    if global_ops["COLOR_FBA_COLORBLIND"]:
        for graph_data in graph_registry.values():
            graph_obj = graph_data[0]
            mvl.colourFBA(graph_obj, fba_results, tlp.Color.Cyan, tlp.Color.Magenta, COLORS["light grey"])
        print("Applied colorblind flux coloring to all graphs. \n")
    if global_ops["COLOR_PATHWAYS"]:
        for graph_data in graph_registry.values():
            graph_obj = graph_data[0]
            # use path_dict for coloring, create path_coloring_dict where keys = reaction ids and values = pathway
            # specific numbers
            path_coloring_dict = dict()
            i = 1
            for path_name, path_data in path_dict.items():
                for re_id in path_data["reactions"]:
                    path_coloring_dict[re_id] = i
                i += 1
            mvl.colourFBA(graph_obj, path_coloring_dict, COLORS["light grey"], discrete=True) # have specific colors!
        print("Applied pathway coloring to all graphs. \n")
    # TODO: have whole map, all pathways that belong to groups are colored (differently), rest grey
    # TODO: return dict: key -> pathway name, value -> assigned color of the patwhay
    # after ALL the other coloring:
    # color metabolites light blue, all reaction nodes light red and all cofactors light yellow (AFTER all other
    # coloring functions!)
    for g_name, g_data in graph_registry.items():
        g_obj = g_data[0]
        g_props = g_data[1]
        colorMetaboliteReactionNodes(g_obj, g_props, tlp.Color(103, 169, 207), # metabolite color
                                                     tlp.Color(252, 146, 114), # reaction color
                                                     tlp.Color(255, 255, 153, 255)) # cofactor color
    print("Applied reaction and metabolite node specific coloring to all graphs. \n")
    # TODO: visualize existence of reactions of a smaller model on a big map (e coli vs mutated e coli eg.)
    # save initial graph with all changes ('messy' layout)
    graph_path = os.path.join(output_path, "{}_graph.tlp".format(model_id))
    tlp.saveGraph(graph_registry['initial graph'][0], graph_path)
    print("Initial graph after coloring saved as: {}_graph.tlp \n".format(model_id))
# Part 3: Exporting
def exportGraph(output_path, graph_registry, global_ops, cDir, file_name, groups_sbml, n2lo, e2lo, ub_sps):
    """
    Save graphs as .tlp, .sbml and/or .sbgn files and use Tulip export plugins to export graphs as .gml,
    .svg and/or .json files::

    - *output_path* an os.path string for all output files
    - *graph_registry* a dict of all generated Tulip graph objects
    - *global_ops* a dict of the current boolean status of all control variables (turn different export plugins on/off)
    - *cDir* current working directory
    - *file_name* the name of all output files (different file types specified by adding graph names to file names
    and the file extensions)

    """
    # save graphs as .tlp files
    if global_ops['EXPORT_TLP'] and global_ops["AUTO_ALIAS"] is False:
        for graph_name, graph_data in graph_registry.items():
            graph_object = graph_data[0]
            if graph_name == "output_whole":
                # save graph with whole network layout
                whole_path = os.path.join(output_path, "{0}_{1}.tlp".format(file_name, graph_name))
                tlp.saveGraph(graph_object, whole_path)
                print("Whole network visualization saved at: {}".format(whole_path))
            if graph_name == "output_selected":
                # save graph with selected pathway layouts
                selected_path = os.path.join(output_path, "{0}_{1}.tlp".format(file_name, graph_name))
                tlp.saveGraph(graph_object, selected_path)
                print("Visualization of selected pathways saved at: {}".format(selected_path))
            if graph_name == "transport_radial_obj_re":
                # save graph with radial transport reactions and obj reaction in the middle
                rad_trans_obj_re_path = os.path.join(output_path, "{0}_{1}.tlp".format(file_name, graph_name))
                tlp.saveGraph(graph_object, rad_trans_obj_re_path)
                print("Visualization of selected pathways saved at: {}".format(rad_trans_obj_re_path))
            if graph_name == "output_grid":
                # save graph with grid pathway layouts
                grid_path = os.path.join(output_path, "{0}_{1}.gml".format(file_name, graph_name))
                tlp.saveGraph(graph_object, grid_path)
                print("Visualization of aligned pathways saved at: {}".format(grid_path))
    # export every graph in graph registry
    if global_ops["EXPORT_GML"]:
        for graph_data in graph_registry.values():
            graph = graph_data[0]
            gml_path = os.path.join(output_path, "{0}_{1}.gml".format(file_name, graph.getName()))
            # replace forward dashes in graph name to not cause problems with path
            gml_path = gml_path.replace('/', '-')
            mvl.exportGML(graph, file_name, gml_path)
            print("Graph exported as: {0}_{1}.gml".format(file_name, graph.getName()))
    if global_ops["EXPORT_SVG"]:
        for graph_data in graph_registry.values():
            graph = graph_data[0]
            # get default parameters of SVG export plugin
            params = tlp.getDefaultPluginParameters('SVG Export', graph)
            # define output path
            svg_path = os.path.join(output_path, "{0}_{1}.svg".format(file_name, graph.getName()))
            # replace forward dashes in graph name to not cause problems with path
            svg_path = svg_path.replace('/', '-')
            tlp.exportGraph('SVG Export', graph, svg_path, params)
            print("Graph exported as: {0}_{1}.svg".format(file_name, graph.getName()))
    if global_ops["EXPORT_JSON"]:
        for graph_data in graph_registry.values():
            graph = graph_data[0]
            # get default parameters of JSON export plugin
            params = tlp.getDefaultPluginParameters('JSON Export', graph)
            # define output path
            json_path = os.path.join(output_path, "{0}_{1}.json".format(file_name, graph.getName()))
            # replace forward dashes in graph name to not cause problems with path
            json_path = json_path.replace('/', '-')
            tlp.exportGraph('JSON Export', graph, json_path, params)
            print("Graph exported as: {0}_{1}.json".format(file_name, graph.getName()))
    # generate SBGN and SBML+Layout Files (Mimoza functions)
    # save graph layout in layout plugin of an SBML file
    if global_ops['EXPORT_SBML']:
        saveLayoutSBML(groups_sbml, ub_sps, n2lo)
    # save graph layout as SBGN file
    if global_ops['EXPORT_SBGN']:
        saveLayoutSBGN(groups_sbml, file_name, output_path, n2lo, e2lo)
# Part 4: Database
def buildSQLiteDatabase(selected_graph, output_path):
    """
    Create an SQLite database to store a graph object (properties of nodes, properties if edges, node connectivity and
    edge connectivity)::

    - *selected_graph* a Tulip graph object
    - *output_path* an os.path for all output files

    """
    # extract data from selected graph
    graph_props_ordered_n_list, graph_props_ordered_e_list, graph_prop_ordered_types = extractGraphPropertyData(
        selected_graph)
    # encode extracted node data
    graph_props_node_tuplelist = encodeGraphPropertyData(selected_graph, graph_props_ordered_n_list,
                                                         graph_prop_ordered_types, nodes=True)
    graph_props_edge_tuplelist = encodeGraphPropertyData(selected_graph, graph_props_ordered_e_list,
                                                         graph_prop_ordered_types, edges=True)
    # build SQLite database from extracted and encoded data (to SQLite compatible types)
    db_path = createSQLiteDB(selected_graph, selected_graph.getName(), graph_props_node_tuplelist,
                             graph_props_edge_tuplelist,
                             output_path)

    return db_path
########################################################################################################################


# individual higher level functions that are called in wrapping functions
def createCBMPYModel(sbml_path, change_flux_bounds=True, ):
    """
    Use CBMPY to read an SBML file and construct a CBMPY model. Compared to the model object created by libSBML,
    a CBMPY model is a new data structure, completely separate from the SBML. This makes it possible to efficiently
    perform refactoring for instance of high integer value flux bounds or non reacting species. Return the CBMPY
    model, the stoichiometric matrix and the id of the objective reaction::

    - *sbml_path* a path to an SBML file
    - *change_flux_bounds* [default=True] refactor all occurrences of very high flux values to inf and -inf
    """
    # create CBMPy model
    cmod = cbmpy.readSBML3FBC(sbml_path)
    # change flux bounds values to not mess with the coloring later on
    if change_flux_bounds:
        cmod.changeAllFluxBoundsWithValue(999999, cbmpy.INF)
        cmod.changeAllFluxBoundsWithValue(-999999, cbmpy.NINF)
        cmod.changeAllFluxBoundsWithValue(1000, cbmpy.INF)
        cmod.changeAllFluxBoundsWithValue(-1000, cbmpy.NINF)
    # remove non-reacting metabolites (since they would be visualized as disconnected nodes)
    if len(cmod.deleteNonReactingSpecies()) != 0:
        cmod.deleteNonReactingSpecies(simulate=False)
    # get stoichiometric matrix
    stM = cmod.N
    # get objective reaction id
    objective_re_id = str(cmod.getActiveObjectiveReactionIds()[0])

    return cmod, stM, objective_re_id

def runFBA(cmod):
    """
    Do an FBA as implemented in CBMPY and save the results to a dictionary::

    - *cmod* a CBMPY model

    """
    # analyze model (FBA) and save results in a flux dictionary
    cbmpy.doFBA(cmod)
    fba_results = {}
    for r in cmod.reactions:
        fba_results[r.getId()] = round(r.getValue(), 6)

    return fba_results

def autoAlias(cmod, global_ops, output_path, spec_split, reac_split, split_biomass=True):
    """
    Alias a number of metabolites/reactions which are selected based on their connectivity. If the amount of connected
    elements is higher than a threshold value, add the respective species to ub_list for Mimoza function node_cloner::

    - *cmod* a CBMPY model
    - *global_ops* a dict with boolean control variables
    - *output_path* an os.path for all output files
    - *spec_split* a threshold value for species splitting (integer)
    - *reac_split* a threshold value for reaction splitting (integer)
    - *split_biomass* [default=True] if true, splits biomass reaction

    """
    # create alias node list based on connectivity <-> user defined thresholds (aliasing by Mimoza)
    ub_sp_list = []
    ub_re_list = []
    if global_ops["AUTO_ALIAS"]:
        sp_list, re_list = mvl.identifyConnectedNodes(cmod, spec_split, reac_split, split_biomass)
        ub_sp_list.extend(sp_list)
        ub_re_list.extend(re_list)
    # combine alias species ids and alias reaction ids lists (for Mimoza node cloning functions)
    ub_list = ub_sp_list + ub_re_list
    # add ubiquitous (alias) species group to sbml file and save it as "..._with_groups.xml"
    # create group for ubiquitous species
    cmod.createGroup('g_ubiquitous_sps')
    g1 = cmod.getGroup('g_ubiquitous_sps')
    g1.setName('ubiquitous species')
    ub_sps = [cmod.getSpecies(i) for i in ub_sp_list]
    ub_react = [cmod.getReaction(i) for i in ub_re_list]
    # add reactions and species to group
    g1.addMember(ub_sps+ub_react)
    # add RDF annotation to group:
    # MIRIAM annotation and SBO term (important for drawing ub species seperately)
    g1.setSBOterm('SBO:0000248')
    g1.addMIRIAMannotation('isDescribedBy', 'ChEBI', 'CHEBI:000000') #dummy annotation
    g1.miriam.isDescribedBy[0] = "http://identifiers.org/ubiquitous" #change dummy to fit mimoza
    cbmpy.writeSBML3FBC(cmod, '{}_with_groups.xml'.format(cmod.getId()), output_path, add_groups=True)
    groups_sbml = os.path.join(output_path, '{}_with_groups.xml'.format(cmod.getId()))

    return ub_list, ub_sps, groups_sbml

def manualPathDict(path_name, sps_id_list, re_id_list):
    """
    Create a dictionary and add one entry containing id lists of reactions and metabolites associated to a manually
    defined pathway in the form:

        path_dict = {
                    "Glycolysis":
                                {
                                "metabolites": ["GLU_6_P", "FRU_6_P", "FRU_BIS_P", "DHAP", "GA_3P", "PG", "PEP"],
                                "reactions": ["R2r", "R3", "R4", "R5r", "R6r", "R7r", "R8r", "R9", "RR9"]
                                }
                    }

    - *path_name* a name of the pathway (string)
    - *sps_id_list* a list of string ids of species of the pathway
    - *re_id_list* a list of string ids of reactions of the pathway

    """

    path_dict = dict()
    path_dict[path_name] = {}
    path_dict[path_name].update({"metabolites": sps_id_list, "reactions": re_id_list})

    return path_dict

def automPathDict(cmod, global_ops, small_threshold=30):
    """
    Automatically create a pathway dictionary for all annotated groups in an CBMPY model. Add all reactions which
    are not yet part of a group to a new group ("undefined")::

    - *cmod* a CBMPY model
    - *global_ops* a dict with boolean control variables
    - *small_treshold* [default=30] an integer threshold level of connections associated to a metabolite for small
    species identification

    """
    path_dict = {}
    if global_ops["PATHWAY_LISTS"]:
        # list of all reaction objects that are associated with a group
        all_groups_member_list = []
        for group in cmod.groups:
            for member in group.members:
                all_groups_member_list.append(member)
        # group all the ungrouped reactions and their metabolites in a new group "undefined" to have
        # them later (otherwise they won't show up in the output graphs since they are not in path_dict)
        cmod.createGroup('undefined')
        undef = cmod.getGroup('undefined')
        undef.setName('Undefined')
        for reaction in cmod.reactions:
            if reaction not in all_groups_member_list:
                undef.addMember(reaction)

        # loop over all groups in iIT341 and automatically create id lists for all groups
        # save all these lists in a new dictionary

        # example:
                  # glyc_sps, glyc_re = mvl.createPathwayLists(cmod, 'Glycolysis')
                  # tca_sps, tca_re = mvl.createPathwayLists(cmod, 'TCA cycle')
                  # fas_sps, fas_re = mvl.createPathwayLists(cmod, 'Fatty Acid Synthesis')
                  # uc_sps, uc_re = mvl.createPathwayLists(cmod, 'Urea Cycle')
                  # am_sps, am_re = mvl.createPathwayLists(cmod, 'Aminosugar Metabolism')
                  # etc. ...

        for group in cmod.groups:
            path_sps, path_re, small_sps_dict = mvl.createPathwayLists(cmod, group.name, small_threshold)
            path_dict.update({group.name: {"metabolites": path_sps, "reactions": path_re}})

        print("Pathway dictionary generated. \n")


        return path_dict, small_sps_dict

def createInitialGraph(graph_name):
    """
    Create a new, empty Tulip graph, instantiate its properties and return the graph object as well as a graph
    properties dictionary.

    - *graph_name* the name of the new graph (string)
    """
    # create new (empty) tulip graph
    graph = tlp.newGraph()
    graph.setName(graph_name)
    # create and store graph properties for new tulip graph elements (nodes and egdes)
    props_tup = mvl.create_props(graph)
    # save graph props in dictionary
    graph_props = {}
    for prop in props_tup:
        graph_props[prop.getName()] = prop

    return graph, graph_props

def compDictForMimoza(cmod):
    """
    Helper function which creates a specific dictionary required to run Mimoza function 'compute_c_id2info'::

    - *cmod* a CBMPY model object

    """
    comp_list = [''] + cmod.getCompartmentIds()
    level_list = [(x, comp_list[x]) for x in range(len(comp_list))]
    c_id2level = dict(zip(comp_list[1:], level_list))

    return c_id2level

# TODO: initial tulip graph creation is slow because of Mimoza Black Box Code!
# TODO: reaction node creation still relies on copied code (problem: reaction node compartment assignment needs helper
# TODO: functions -> understand what the copied code does and rewrite!)
#START COPIED CODE######################################################################################################
# function copied from Mimoza
def compute_c_id2info(c_id2level, input_model):
    c_id2info = {}
    for comp in input_model.compartments:
        c_id = comp.getId()
        c_name = comp.getName()
        if not c_name:
            c_name = c_id
        c_id2info[c_id] = (c_name, '', c_id2level[c_id])
    c_id2outs = {}
    for c_id in c_id2info.keys():
        _, _, (_, out_c_id) = c_id2info[c_id]
        outs = []
        while out_c_id:
            outs.append(out_c_id)
            _, _, (_, out_c_id) = c_id2info[out_c_id]
        c_id2outs[c_id] = outs
    return c_id2info, c_id2outs

# helper function for reaction compartment assignment (copied from Mimoza)
def get_r_comp(all_comp_ids, c_id2level, c_id2outs):
    if len(all_comp_ids) == 1:
        return all_comp_ids.pop()
    get_level = lambda c_id: c_id2level[c_id][0]
    outer_most = min(all_comp_ids, key=get_level)
    inner_most = max(all_comp_ids, key=get_level)
    outer_level, inner_level = get_level(outer_most), get_level(inner_most)
    if outer_level == inner_level or (outer_most not in c_id2outs[inner_most]):
        candidates = set(c_id2outs[inner_most]) & set(c_id2outs[outer_most])
        if candidates:
            return max(candidates, key=get_level)
        else:
            return outer_most
    if inner_level - outer_level > 1:
        return max(c_id2outs[inner_most], key=get_level)
    return outer_most
#END COPIED CODE########################################################################################################

def mapModelToGraph(cmod, graph, graph_props):
    """
    Map model objects (metabolites and reactions) to generated graph objects (nodes) and save in a dictionary. Use
    resulting dictionary to track whether or not an initially generated node has already been used for layouting (in
    which case a new, identical node needs to be generated by other functions)::

    - *cmod* a CBMPY model
    - *graph* a TULIP graph
    - *graph_props* a dict containing properties of the TULIP graph

    """
    graph_elements = {}
    # all species and all reactions of sbml model = total set of model elements
    model_elements = cmod.species + cmod.reactions
    # dict with keys: model ids and values = dict with keys: pathway names and values = associated graph nodes
    for m_e in model_elements:
        graph_elements[m_e.getId()] = {}
        # use [0] to get first entry of that list so that the returned object is the nodes itself and not a list again
        graph_elements[m_e.getId()]["initial_node"] = [n for n in graph.getNodes() if graph_props["id"][n] ==
                                                       m_e.getId()][0]
        # edges do not exist in SBML model! (and created graph edges have no ids, just names)
        # -> edges are not included in dict!
        # track if graph element was used for layouting already (default value: False, gets changed to True once the
        # respective node was layouted in a subgraph)
        graph_elements[m_e.getId()]["layout"] = False

    return graph_elements

def objReNodeList(graph, graph_props, objective_re_id):
    """
    Create a list with all nodes that correspond to ids of species associated with the objective reaction::

    - *graph* a TULIP graph
    - *graph_props* a dict with graph properties
    - *objective_re_id* the id of the objective reaction (string)

    """
    # create objective reaction nodes list:
    obj_re_nodes = []
    for n in graph.getNodes():
        # select obj reaction node
        if graph_props["id"][n] == objective_re_id:
            # add obj reaction node to the list
            obj_re_nodes.append(n)
        # iterate over all nodes which are connected to obj reaction node and add them to the list
        # the new PEP node will be added instead of the old(raw) node because the connection to the raw node was cut and a
            # new connection to the new PEP node was created
            for node in graph.getInOutNodes(n):
                obj_re_nodes.append(node)

    return obj_re_nodes

def labelGraphNodes(graph, label_type):
    """
    Add Labels to graph nodes (id or name)::

    - *graph* a Tulip graph
    - *label_type* labels are based either on id or name of elements

    """
    # add labels to nodes
    # define node types for labeling
    NODE_TYPES = {'metabolites': 1,
                  'reactions'  : 2
                  }
    if label_type == "id":
        mvl.addLabels(graph, NODE_TYPES['metabolites'], "id")
        mvl.addLabels(graph, NODE_TYPES['reactions'], "id")
    elif label_type == "name":
        mvl.addLabels(graph, NODE_TYPES['metabolites'], "name")
        mvl.addLabels(graph, NODE_TYPES['reactions'], "name")

def colorAllGraphElements(graph, graph_element, color):
    """
    Color graph elements (nodes and/or edges). Returns a dict with all possible colors so that other functions can
    also use them::

    - *graph* a TULIP graph
    - *graph_element* either "nodes" or "edges"
    - *color* a color that is defined in dict COLORS

    """
    # change node or edge colour
    # define colours
    COLORS = {'red'  : tlp.Color(255, 0, 0),
               'green': tlp.Color(0, 255, 0),
               'blue' : tlp.Color(0, 0, 255),
               'white': tlp.Color(255, 255, 255),
               'grey' : tlp.Color(220, 220, 220),
               'light grey': tlp.Color(235, 235, 235),
               'black': tlp.Color(0, 0, 0)
               }

    # set colour of all edges (complete graph)
    if graph_element == "nodes":
        mvl.setColours(graph, graph.getNodes(), COLORS[color])
    elif graph_element == "edges":
        mvl.setColours(graph, graph.getEdges(), COLORS[color])

    return COLORS

def layoutAllPathwaysFromGroups(cmod, global_ops, path_dict, graph, graph_elements, graph_props, small_sps_dict):
    """
    Generate subgraphs for all pathways contained in path_dict and create subgraph layouts::

    - *cmod* a CBMPY model
    - *global_ops* a dict with global control variables
    - *path_dict* a dict with pathways extracted from groups
    - *graph* a Tulip graph
    - *graph_elements* a dict of graph nodes mapped to model ids
    - *graph_props* a dict of graph properties
    - *small_sps_dict* a dict of small species of each reaction

    """
    # layout all pathways in graph (entire path_dict) to generate subgraphs
    for path_name, path_ids in path_dict.items():
        # fetch pathway ids
        path_list = list(path_ids["metabolites"]) + list(path_ids["reactions"])
        # choose appropriate layout algorithm
        path_layout = []
        if "ycle" in path_name:
            path_layout.append("Circular (OGDF)")
        else:
            path_layout.append("Sugiyama (OGDF)")
        # apply layout to pathway and track already layouted graph elements
        if path_name in global_ops["ADD_COFACTORS"]:
            path_sg = mvl.layoutPathway(cmod, graph, graph_elements, graph_props, path_list, path_name,
                                                path_layout[0], small_sps_dict, add_cofactors=True)
        else:
            path_sg = mvl.layoutPathway(cmod, graph, graph_elements, graph_props, path_list, path_name,
                                                path_layout[0], small_sps_dict, add_cofactors=False)
        # set name of subgraph
        path_sg.setName(path_name)

        # test if all reaction and metabolite ids of the group definition are part of the generated subgraph
        # unique set of ids of generated graph nodes
        generated_nodes = set([graph_props["id"][n] for n in path_sg.getNodes()])
        # unique set of input ids
        input_ids = set(path_list)
        # both need to have the same length since at least one node needs to exist for every unique input id
        # TODO. after having tested identity
        # TODO: -> do an intersection with the two sets (to catch duplications, etc.)
        if len(generated_nodes) != len(input_ids):
            print("{} subgraph incomplete! \n".format(path_name))

def pathwayDump(graph, global_ops, graph_registry, graph_props, selected_pathways, grid_pathways):
    """
    Define new output graphs and add them to graph registry::

    - *graph* a parent Tulip graph to take subgraphs from
    - *global_ops* a dict containing global control variables
    - *graph_registry* a dict containing all generated graph objects (keys: graph names)
    - *graph_props* a dict containing the properties of the parent graph
    - *selected_pathways* a list of selected pathways, positioning auto (list of strings)
    - *grid_pathways* a list of pathways for grid layout, positioning controlled (list of strings)

    """
    # 1.) default dump (all pathways in one SVG file)
    # -> combine all pathways to create a new output graph of the whole network
    if global_ops["TRANSPORT_RADIAL"]:
        whole_output_graph, whole_output_props = mvl.dumpToGraph(graph, graph_props, "output_whole",
                                             [sg.getName() for sg in graph.getSubGraphs()])
        graph_registry.update({whole_output_graph.getName(): [whole_output_graph,whole_output_props]})
    else:
        whole_output_graph, whole_output_props = mvl.dumpToGraph(graph, graph_props, "output_whole",
                                                                 [sg.getName() for sg in graph.getSubGraphs()],
                                                                 transport_radial=False)
        graph_registry.update({whole_output_graph.getName(): [whole_output_graph, whole_output_props]})

    # 2.) selected dump (selected pathways in one SVG file)
    # create a new output graph with a few selected pathways combined
    if selected_pathways:
        selected_output_graph, selected_output_props = mvl.dumpToGraph(graph, graph_props, "output_selected",
                                                                       selected_pathways)
        graph_registry.update({selected_output_graph.getName(): [selected_output_graph,selected_output_props]})

    # 3.) sequential dump (all pathways in individual SVG files)
    individual_output_graphs = []
    for sg in graph.getSubGraphs():
        # append the returned output graph for each input subgraph to list
        sg_output_graph, sg_output_props = mvl.dumpToGraph(graph, graph_props, sg.getName(), [sg.getName()])
        # TODO: add pathway label to sg_out_graph
        individual_output_graphs.append(sg_output_graph)
        graph_registry.update({sg_output_graph.getName(): [sg_output_graph, sg_output_props]})

    # 4.) radial transport reactions with obj reaction in the middle
    radial_trans_obj_re_graph, radial_trans_obj_re_props = mvl.dumpToGraph(graph, graph_props, "transport_radial_obj_re",
                                                                           ["objective reaction", "Transport"],                                                                                          center_obj_trans_re=True)
    graph_registry.update({radial_trans_obj_re_graph.getName(): [radial_trans_obj_re_graph, radial_trans_obj_re_props]})

    # 5.) grid dump (select pathways AND control their location)
    if grid_pathways:
        grid_output_graph, grid_output_props = mvl.dumpToGraph(graph, graph_props, "output_grid",
                                                               grid_pathways, grid_positioning=True)
        graph_registry.update({grid_output_graph.getName(): [grid_output_graph, grid_output_props]})

    # TODO: do - place 4 subgraphs in a new graph and control their location!
    # TODO:    - all small pathways (and labels), less than 10 reactions, layed out on a grid
    # TODO: i need to be able to control the location -> write a function to create grid layout (somehow)

def colorMetaboliteReactionNodes(graph, graph_props, met_color, reac_color, cofac_color):
    """
    Give metabolite, reaction and cofactor nodes distinct colors (call after all other coloring functions to not
    interfere with them!)

    - *graph* a Tulip graph
    - *graph_props* a dict of graph properties
    - *met_color* a Tulip color object (wraped around RGB value) for metabolite nodes
    - *reac_color* a Tulip color object (wraped around RGB value) for reaction nodes
    - *cofac_color* a Tulip color object (wraped around RGB value) for cofactor nodes

    """
    # change node color of metabolite and reaction nodes
    for n in graph.getNodes():
        if graph_props["type"][n] == 1:  # metabolite node
            graph_props["viewColor"][n] = met_color
        if graph_props["type"][n] == 2:  # reaction node
            graph_props["viewColor"][n] = reac_color
        if graph_props["type"][n] == 3:  # cofactor node
            graph_props["viewColor"][n] = cofac_color

def saveLayoutSBML(groups_sbml, ub_sps, n2lo):
    """
    Save graph layout in the layout plugin of an SBML file::

    - *groups_sbml* an SBML file with groups
    - *ub_sps* a list of ubiquitous species objects generated by autoAlias
    - *n2lo* a dict with node layout information from Mimoza function tulip_graph2geojson

    """
    # add layout to SBML file
    SBML_reader = libsbml.SBMLReader()
    groups_document = SBML_reader.readSBML(groups_sbml)
    groups_model = groups_document.getModel()
    gen_model = False
    sbml_manager.save_as_layout_sbml(groups_model, gen_model, groups_sbml, n2lo, ub_sps)

def saveLayoutSBGN(groups_sbml, file_name, output_path, n2lo, e2lo):
    """
    Save graph layout as SBGN file in current working directory cDir::

    - *groups_sbml* an SBML file with groups
    - *file_name* the name of the SBGN file (string)
    - *output_path* an os.path for all output files
    - *n2lo* a dict with node layout information
    - *e2lo* a dict with edge layout information

    """
    # create SBGN file
    # TODO: check SBGN: is it valid?, what information is stored?
    # TODO: -> Error on: http://sysbioapps.dyndns.org/Layout
    """
    Error

    Object reference not set to an instance of an object.

   at SBMLExtension.Formats.LibSbgnML.ReadLayoutsFrom(XmlDocument doc)
   at SBMLExtension.Formats.LayoutFormats.ReadLayouts(XmlDocument oDoc)
   at SBMLExtension.Util.readLayout(XmlDocument oDoc, Boolean overWriteWithSBGNDefault)
   at SBMLExtension.Util.readLayout(String sXML, Boolean overWriteWithSBGNDefault)
   at Welcome.Controllers.LayoutController.UpdateImage(String file, Boolean sbgn)

    If this problem persists, please notify the author.
    """
    groups_sbgn = os.path.join(output_path, '{}.sbgn'.format(file_name))
    # TODO: since Mimoza sbgn export function gets two dicts (n2lo and e2lo) with graph info instead of a graph object
    # TODO: can we link this function to the other output graphs?
    # TODO: maybe write code that extracts dicts like n2lo and e2lo from any graph object!
    # TODO: whats the structure of these dicts?
    SBML_reader = libsbml.SBMLReader()
    groups_document = SBML_reader.readSBML(groups_sbml)
    groups_model = groups_document.getModel()
    sbgn_helper.save_as_sbgn(n2lo, e2lo, groups_model, groups_sbgn)

def extractGraphPropertyData(graph):
    """
    Extract all node and edge data associated with all properties of a Tulip graph and save it in lists::

    - *graph* a Tulip graph

    """
    # create a dictionary of dicts for the graph properties to link types of properties to their data:
    # dict setup: keys = property names with values = dict with 1. key 'type' => type of property
    #                                                       and 2. key 'data' => list of data
    #                                                                              (same order as columns in graph object)
    #
    # for SQLite: convert boolean values to integers (True = 1, False = 0) since sqlite has no boolean data type!
    graph_props_node_dict = {}
    graph_props_edge_dict = {}
    graph_props_ordered_n_list = []
    graph_props_ordered_e_list = []
    graph_prop_ordered_types = []
    for prop in graph.getObjectProperties():
        # add types of props to a list
        graph_prop_ordered_types.append(prop.getTypename())
        # add property data of nodes to list
        prop_node_sublist = []
        for n in graph.getNodes():
            # convert True/False to 1/0
            if prop.getTypename() == "bool":
                prop_node_sublist.append(prop[n] * 1)
            # convert all instances of _tulip.Vec3f, _tulip.Size and _tulip.Color objects to tuples
            elif prop.getTypename() == "layout" or prop.getTypename() == "color" or prop.getTypename() == "size":
                prop_node_sublist.append(tuple(prop[n]))
            else:
                prop_node_sublist.append(prop[n])
        # add graph node prop info to ordered list
        graph_props_ordered_n_list.append(prop_node_sublist)
        # fill dict with graph node prop info
        graph_props_node_dict[prop.getName()] = {}
        graph_props_node_dict[prop.getName()]["type"] = prop.getTypename()
        graph_props_node_dict[prop.getName()]["data"] = prop_node_sublist
        # add property data of edges to list
        prop_edge_sublist = []
        for e in graph.getEdges():
            if prop.getTypename() == "bool":
                prop_edge_sublist.append(prop[e] * 1) # convert True/False to 1/0
            # convert all instances of _tulip.Vec3f, _tulip.Size and _tulip.Color objects to tuples
            elif prop.getTypename() == "layout" or prop.getTypename() == "color" or prop.getTypename() == "size":
                prop_edge_sublist.append(tuple(prop[e]))
            else:
                prop_edge_sublist.append(prop[e])
        # add graph edge prop info to ordered list
        graph_props_ordered_e_list.append(prop_edge_sublist)
        # fill dict with graph edge prop info
        graph_props_edge_dict[prop.getName()] = {}
        graph_props_edge_dict[prop.getName()]["type"] = prop.getTypename()
        graph_props_edge_dict[prop.getName()]["data"] = prop_edge_sublist

    return graph_props_ordered_n_list, graph_props_ordered_e_list, graph_prop_ordered_types

def encodeGraphPropertyData(graph, graph_props_ordered_data, graph_prop_ordered_types, nodes=None, edges=None):
    """
    Encode node/edge data of all graph properties to SQLite compatible data types (INTEGER, REAL or TEXT). Specify
    type of graph elements with 'nodes'/'edges' keywords since lists of nodes and egdes can have different sizes::

    - *graph* a Tulip graph
    - *graph_props_ordered_data* a list of node/edge data associated with all graph props
    - *graph_prop_ordered_types* a list of graph property types
    - *nodes* [default=None] include parameter to encode node properties
    - *edges* [default=None] include parameter to encode edge properties

    """
    # encode graph property data to SQLite compatible types
    # go through the list of graph property lists:
    # since a sublist contains only one graph property for all nodes/edges -> the first element of the sublist is
    # exemplary for the whole list
    graph_props_tuplelist = list() # after encoding is finished: save encoded prop lists in tuples to prepare for
    # INSERT
    # add node ids to beginning of tuple list to serve as primary key in table
    if nodes:
        graph_props_tuplelist.append([i for i in range(len(list(graph.getNodes())))])
    elif edges:
        graph_props_tuplelist.append([i for i in range(len(list(graph.getEdges())))])
    else:
        print("No graph element type specified! \n")
    #  data types that require encoding
    # every item of properties with these types is a list-like thing itself!
    # 'size' contains elements of type: _tulip.Size
    # 'color' contains elements of type: _tulip.Color
    # 'layout' and 'vector<string>' contain elements of type: _tulip.Vec3f
    tulip_prop_types = ('size', 'color', 'layout', 'vector<string>')
    # replace items of sublists that need to be encoded
    for prop_type, prop_data_list in zip(graph_prop_ordered_types, graph_props_ordered_data):
        # allowed data types, no encoding needed
        if prop_type == 'string':
            pass
        if prop_type == 'int':
            pass
        if prop_type == 'double': # in python: float (in SQL -> REAL)
            pass
        if prop_type == 'bool': # was converted to INTEGER
            pass
        # layout contains: tuples of _tulip.Vec3f objects!
        # ex: ([50.0,60.0],)
        if prop_type == 'layout':
            for i, item in enumerate(prop_data_list):
                # convert _tulip.Vec3f object inside tuple to a list inside tuple
                # ex: from ([50.0,60.0],) -> to [[50.0, 60.0]]
                list_item = list()
                for e in item:
                    if isinstance(e, float):
                        list_item.append(e)
                    else:
                        # if e is not a float is a fucking tulip vector thing -> use .get() to return tuple
                        # since e now contains an entire tuple => its the entire list item
                        # already, no need for appending -> just cast it as a list and save it as the whole list item
                        list_item = list((e.get()))
                # serialize list with pickle (returns string)
                str_item = pickle.dumps(list_item)
                prop_data_list[i] = str_item
        # convert data of a type that requires encoding to string
        if prop_type in tulip_prop_types and prop_type != 'layout':
            for i, item in enumerate(prop_data_list):
                try:
                    prop_data_list[i] = pickle.dumps(item, protocol=0)
                except Exception as ex:
                    print(ex)
                    print(prop_type, item, type(item))
        # exlucde props with type 'graph' from tuplelist (viewMetaGraph)
        included = True
        if prop_type == 'graph':
            included = False
        # final step: save included, encoded sublist as tuple
        if included is True:
            graph_props_tuplelist.append(tuple(prop_data_list))

    return graph_props_tuplelist

def createSQLiteDB(graph, db_name, graph_props_node_tuplelist, graph_props_edge_tuplelist, output_path):
    """
    Create an SQLite database storing the encoded graph props data of nodes, edges as well as their connectivity::

    - *graph* a Tulip graph
    - *db_name* the name of the database
    - *graph_props_node_tuplelist* a list of tuples containing all graph prop data of one node each
    - *graph_props_edge_tuplelist* a list of tuples containing all graph prop data of one edge each
    - *output_path* an os.path for all output files

    """
    # create an SQL database with sqlite3
    # connect to database
    db_path = os.path.join(output_path, '{0}_graph.db'.format(db_name))
    db = sqlite3.connect(db_path)
    # get a cursor object
    c = db.cursor()
    # pass SQL statements into cursor
    # create empty graph node table (skips if table already exists)
    c.execute('''CREATE TABLE IF NOT EXISTS graph_nodes(node_id INTEGER PRIMARY KEY,
    anc_id TEXT,
    anc_name TEXT,
    anc_term TEXT,
    c_id TEXT,
    clone_id TEXT,
    edge_directionality TEXT,
    f TEXT,
    fake INTEGER,
    flux REAL,
    id TEXT,
    name TEXT,
    pathway TEXT,
    rel_c_ids TEXT,
    rev INTEGER,
    stoich REAL,
    term TEXT,
    tr INTEGER,
    type INTEGER,
    ub INTEGER,
    viewBorderColor TEXT,
    viewBorderWidth REAL,
    viewColor TEXT,
    viewFont TEXT,
    viewFontAwesomeIcon TEXT,
    viewFontSize INTEGER,
    viewLabel TEXT,
    viewLabelBorderColor TEXT,
    viewLabelBorderWidth REAL,
    viewLabelColor TEXT,
    viewLabelPosition INTEGER,
    viewLayout TEXT,
    viewMetric REAL,
    viewRotation REAL,
    viewSelection INTEGER,
    viewShape INTEGER,
    viewSize INTEGER,
    viewSrcAnchorShape INTEGER,
    viewSrcAnchorSize TEXT,
    viewTexture TEXT,
    viewTgtAnchorShape INTEGER,
    viewTgtAnchorSize TEXT)
    ''')
    # create empty graph edge table (skips if table already exists)
    c.execute('''CREATE TABLE IF NOT EXISTS graph_edges(edge_id INTEGER PRIMARY KEY,
    anc_id TEXT,
    anc_name TEXT,
    anc_term TEXT,
    c_id TEXT,
    clone_id TEXT,
    edge_directionality TEXT,
    f TEXT,
    fake INTEGER,
    flux REAL,
    id TEXT,
    name TEXT,
    pathway TEXT,
    rel_c_ids TEXT,
    rev INTEGER,
    stoich REAL,
    term TEXT,
    tr INTEGER,
    type INTEGER,
    ub INTEGER,
    viewBorderColor TEXT,
    viewBorderWidth REAL,
    viewColor TEXT,
    viewFont TEXT,
    viewFontAwesomeIcon TEXT,
    viewFontSize INTEGER,
    viewLabel TEXT,
    viewLabelBorderColor TEXT,
    viewLabelBorderWidth REAL,
    viewLabelColor TEXT,
    viewLabelPosition INTEGER,
    viewLayout TEXT,
    viewMetric REAL,
    viewRotation REAL,
    viewSelection INTEGER,
    viewShape INTEGER,
    viewSize INTEGER,
    viewSrcAnchorShape INTEGER,
    viewSrcAnchorSize TEXT,
    viewTexture TEXT,
    viewTgtAnchorShape INTEGER,
    viewTgtAnchorSize TEXT)
    ''')
    # create empty connectivity table (each row is an edge with id, source node id and target node id as columns)
    c.execute('''CREATE TABLE IF NOT EXISTS edge_connectivity(edge_id INTEGER PRIMARY KEY,
    source_node_id TEXT,
    target_node_id TEXT)
    ''')
    # create node connectivity table: for each node a list (= string) of all connected In and Out edges
    c.execute('''CREATE TABLE IF NOT EXISTS node_connectivity(node_id INTEGER PRIMARY KEY,
    in_out_edge_list TEXT)
    ''')
    # create edge color table: for each edge multiple columns with color set tuples->(flux, edge color)
    c.execute('''CREATE TABLE IF NOT EXISTS edge_color(edge_id INTEGER PRIMARY KEY,
    color_set1 TEXT,
    color_set2 TEXT,
    color_set3 TEXT)
    ''')
    # reorder tuple list for INSERT INTO statement:
    # instead of [(all prop1), (all prop2), ...] I want [(prop1[1], prop2[1], ...), (prop1[2], prop2[2], ...), ...]
    graph_props_node_insertlist = []
    graph_props_edge_insertlist = []
    for i in range(len(list(graph.getNodes()))):
        graph_props_node_insertlist.append([tup[i] for tup in graph_props_node_tuplelist])
    for i in range(len(list(graph.getEdges()))):
        graph_props_edge_insertlist.append([tup[i] for tup in graph_props_edge_tuplelist])
    # fill tables with encoded values
    # table graph_node
    for l in graph_props_node_insertlist:
        c.execute('''INSERT OR IGNORE INTO graph_nodes VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,
        ?,?,?,?,?,?,?,?,?,?,?,?,?)''', l)
    # table graph_edges
    for l in graph_props_edge_insertlist:
        c.execute('''INSERT OR IGNORE INTO graph_edges VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,
        ?,?,?,?,?,?,?,?,?,?,?,?,?)''', l)
    # table edge_connectivity
    # create dict with key = str(node) and value = integer id of the node
    node_str_int_dict = {}
    for node, i in zip(graph.getNodes(), range(len(list(graph.getNodes())))):
        node_str_int_dict[str(node)] = i
    for edge, i in zip(graph.getEdges(), range(len(list(graph.getEdges())))):
        insert_row = [i, node_str_int_dict[str(graph.source(edge))], node_str_int_dict[str(graph.target(edge))]]
        c.execute('''INSERT OR IGNORE INTO edge_connectivity VALUES (?,?,?)''', insert_row)
    # table node_connectivity
    for node, i in zip(graph.getNodes(), range(len(list(graph.getNodes())))):
        row = [i, str(list(graph.getInOutEdges(node)))]
        c.execute('''INSERT OR IGNORE INTO node_connectivity VALUES (?,?)''', row)
    # table edge_color
    # create an ordered list of prop names of the graph
    graph_prop_ordered_names = []
    for prop in graph.getObjectProperties():
        graph_prop_ordered_names.append(prop.getName())
    # create dict with key = prop name and value = tuple of prop data for all edges
    prop_name_data_dict = {}
    for prop_name, prop_tup in zip(graph_prop_ordered_names, graph_props_edge_tuplelist[1:]):
        prop_name_data_dict[prop_name] = prop_tup
    # create color set tuples (flux value, edge color)
    # TODO: figure out a way of adding additional color sets if the whole pipeline only applies one coloring at a time
    # TODO: for now: just add tuples with empty strings as placeholder as color sets 2 and 3
    for i in range(len(list(graph.getEdges()))):
        db_row = [i, (prop_name_data_dict["flux"][i], prop_name_data_dict["viewColor"][i]),
                     (prop_name_data_dict["flux"][i], ""),
                     (prop_name_data_dict["flux"][i], "")]
        # pickle the color set tuples
        for i,item in enumerate(db_row):
            if isinstance(item, tuple):
                str_item = pickle.dumps(item)
                db_row[i] = str_item
        # add the contents of the row list [integer, pickled tuples] to sqlite db
        c.execute('''INSERT OR IGNORE INTO edge_color VALUES (?,?,?,?)''', db_row)
    print("Graph data has been stored in SQL database. \n")
    # commit changes to database
    db.commit()
    # close connection to database
    db.close()

    return db_path

def retrieveGraphElementData(db_path, table_name):
    """
    Retrieve and decode node and egde data of a stored tulip graph as lists of lists (sublists are rows).
    Also return a list of column headers::

    - *db_path* an os.path of the database
    - *table_name* a string table name contained in the current database

    """
    # connect to database
    db = sqlite3.connect(db_path)
    # get a cursor object
    c = db.cursor()
    # retrieve entire table
    c.execute('''SELECT * FROM {}'''.format(table_name))
    graph_rows_list = c.fetchall() # list of row tuples
    # change list of row tuples to list of row lists
    graph_rows_list = [list(tup) for tup in graph_rows_list]
    # sqlite stores strings in unicode!
    # -> 'test' is automatically converted to u'test'
    # convert unicode strings to ascii strings
    for ind, row_list in enumerate(graph_rows_list):
        graph_rows_list[ind] = [str(item) if isinstance(item, unicode) else item for item in row_list]
    # decode queried data (string -> list)
    for row_list in list(graph_rows_list):
        for ind, item in enumerate(row_list):
            # go over every item in a row list. for every string item:
            # if the string was generated from a tup of multiple numbers (x,x,x,x)
            # -> convert ascii string back to list
            # -> identify the strings that originate from tuples by only taking strings without letters
            if isinstance(item, str):
                item_list = []
                if ',' in item and any(e.isalpha() for e in item) == False:
                    # convert item string to a list
                    item_list.append([int(string_part) if isinstance(string_part, int) else float(string_part)
                         for string_part in item.split(',')])
                # replace encoded item with decoded item (if it was converted, otherwise keep old entry)
                if item_list:
                    row_list[ind] = item_list
                else:
                    row_list[ind] = item
    # retrieve a list of column headers
    col_headers = [description[0] for description in c.description]
    # commit changes to database
    db.commit()
    # close connection to database
    db.close()

    return graph_rows_list, col_headers

def retrieveDataFromDatabaseTable(db_path, table_name):
    """
    Retrieve and decode data from a Tulip graph stored in an SQLite database. Each row of an SQLite db table is saved in
    a list and all row lists are appended and returned as items of a nested list.Also return a list of column headers::

    - *db_path* an os.path of the database
    - *table_name* name of a table inside the db (string)

    """
    # retrieve encoded data
    db_graph_data_list, col_headers = retrieveGraphElementData(db_path, table_name)
    # decode data:
    # unpickle the content of a list returned from the db
    # returned lists contain row lists -> ordered list of all properties for one node
    for row_list in db_graph_data_list:
        for i, prop_value in enumerate(row_list):
            try:
                list_obj = pickle.loads(prop_value)
                row_list[i] = list_obj
            except:
                pass

    return db_graph_data_list, col_headers

def buildTulipGraphFromSQL(db_graph_node_list, db_graph_edge_list, db_graph_edgeconn_list, node_col_headers,
                           edge_col_headers, graph_name, output_path):
    """
    Generate a new Tulip graph based on data retrieved from an SQLite database. Save and return generated graph::

    - *db_graph_node_list* a list of row sublists containing all node property data
    - *db_graph_edge_list* a list of row sublists containing all edge property data
    - *db_graph_edgeconn_list* a list of row sublists containing all edges and their source and target node ids
    - *node_col_headers* a list of column headers of the graph_nodes table
    - *edge_col_headers* a list of column headers of the graph_edges table
    - *graph_name* a name for the new Tulip graph (string)
    - *output_path* an os.path for all output files

    """
    # create new graph to reconstruct network from data retrieved from db
    db_graph = tlp.newGraph()
    db_graph_props = mvl.create_props(db_graph)
    db_graph.setName(graph_name)
    # add new, empty nodes to graph
    db_graph.addNodes(len(db_graph_node_list))
    # convert data in db_graph_node_list back to tulip types
    for row_sublist in db_graph_node_list:
        # remove ids (first entry) from sublists of db_graph_node list to have equal amount of props again
        row_sublist.pop(0)
        # get a list of column headers (which will have same length and order as row list!)
        # if I want an entry of row_sublist: row_sublist[col_header_list.index("viewColor")] = viewColor data
        # go through each sublist of db retreived list and change entries of row list back to tulip specific types if
        #  necessary
        # since first column gets removed, the index needs to be always shifted back by 1
        for i, prop_item in enumerate(row_sublist):
            if prop_item == row_sublist[node_col_headers.index("viewColor")-1] or \
               prop_item == row_sublist[node_col_headers.index("viewBorderColor")-1] or \
               prop_item == row_sublist[node_col_headers.index("viewLabelColor")-1] or \
               prop_item == row_sublist[node_col_headers.index("viewLabelBorderColor")-1]:
                # it needs to be converted to tulip.Color
                row_sublist[i] = tlp.Color(prop_item[0], prop_item[1], prop_item[2], prop_item[3])
            elif prop_item == row_sublist[node_col_headers.index("viewSize")-1] or \
                 prop_item == row_sublist[node_col_headers.index("viewSrcAnchorSize")-1] or \
                 prop_item == row_sublist[node_col_headers.index("viewTgtAnchorSize")-1]:
                # it needs to be converted to tulip.Size
                row_sublist[i] = tlp.Size(prop_item[0], prop_item[1], prop_item[2])
            elif prop_item == row_sublist[node_col_headers.index("viewLayout")-1]:
                # it needs to be converted to tulip.Layout
                row_sublist[i] = tlp.Vec3f(prop_item[0], prop_item[1], prop_item[2])
    # overwrite the empty/default prop values of new nodes with the prop values retrieved from db
    for n, n_db_prop_sublist in zip(db_graph.getNodes(), db_graph_node_list):
        for empty_prop, n_db_prop in zip(db_graph_props, n_db_prop_sublist):
            empty_prop[n] = n_db_prop
    # prepare edge connectivity data to add edges:
    for i, sublist in enumerate(db_graph_edgeconn_list):
        # convert node id strings to [source id, target id] lists from extracted edge connectivity list:
        # [[edge_id, "source node id string", "target node id string"], ..] -> [[edge 1, [int(src id), int(tgt id)]], ..]
        db_graph_edgeconn_list[i] = [sublist[0], [int(sublist[1]), int(sublist[2])]]
    # create dict with keys = node ids and values = node objects
    node_id_to_obj_dict = dict()
    for n_id, n in enumerate(list(db_graph.getNodes())):
        node_id_to_obj_dict[n_id] = n
    # go through edge connectivity list and replace node ids in tuples with node objects
    for i, sublist in enumerate(db_graph_edgeconn_list):
        for item in sublist:
            if isinstance(item, list):
                # get node objects that correspond to the node id
                node_objs = [node_id_to_obj_dict[e] for e in item]
                # convert source, target lists to tuples for addEdges function
                node_objs_tup = tuple(node_objs)
                db_graph_edgeconn_list[i] = node_objs_tup
                # delete first item (edge id) - which is an int
                # if isinstance(item, int):
                #    sublist.remove(item)
    # add edges to db_graph - function tlp.addEdges takes scr and tgt node tuples as parameter
    db_graph.addEdges(db_graph_edgeconn_list)
    # TODO: also add props to edges (also save them in SQL db)
    # convert data in db_graph_edge_list back to tulip types
    for row_sublist in db_graph_edge_list:
        # remove ids (first entry) from sublists of db_graph_edge_list to have equal amount of props again
        row_sublist.pop(0)
        # get a list of column headers (which will have same length and order as row list!)
        # if I want an entry of row_sublist: row_sublist[col_header_list.index("viewColor")] = viewColor data
        # go through each sublist of db retreived list and change entries of row list back to tulip specific types if
        #  necessary
        # since first column gets removed, the index needs to be always shifted back by 1
        for i, prop_item in enumerate(row_sublist):
            if prop_item == row_sublist[edge_col_headers.index("viewColor")-1] or \
                            prop_item == row_sublist[edge_col_headers.index("viewBorderColor")-1] or \
                            prop_item == row_sublist[edge_col_headers.index("viewLabelColor")-1] or \
                            prop_item == row_sublist[edge_col_headers.index("viewLabelBorderColor")-1]:
                # it needs to be converted to tulip.Color
                row_sublist[i] = tlp.Color(prop_item[0], prop_item[1], prop_item[2], prop_item[3])
            elif prop_item == row_sublist[edge_col_headers.index("viewSize")-1] or \
                            prop_item == row_sublist[edge_col_headers.index("viewSrcAnchorSize")-1] or \
                            prop_item == row_sublist[edge_col_headers.index("viewTgtAnchorSize")-1]:
                # it needs to be converted to tulip.Size
                row_sublist[i] = tlp.Size(prop_item[0], prop_item[1], prop_item[2])
            #elif prop_item == row_sublist[edge_col_headers.index("viewLayout")-1] and \
            #     len(row_sublist[edge_col_headers.index("viewLayout") - 1]) > 0:
                # it needs to be converted to tulip.Layout (but only if the coords list is not empty)
            #    print(row_sublist[edge_col_headers.index("viewLayout")-1])
            #    row_sublist[i] = tlp.Vec3f(prop_item[0], prop_item[1], prop_item[2])
    # overwrite the empty/default prop values of new edge with the prop values retrieved from db
    for e, e_db_prop_sublist in zip(db_graph.getEdges(), db_graph_edge_list):
        for empty_prop, e_db_prop in zip(db_graph_props, e_db_prop_sublist):
            if empty_prop.getName() != "viewLayout":
                empty_prop[e] = e_db_prop
    # save reconstructed db graph
    recon_db_graph_path = os.path.join(output_path, 'reconstructed_graph_from_db.tlp')
    tlp.saveGraph(db_graph, recon_db_graph_path)

    return db_graph


#if __name__ == '__main__':
# execute this if in a console (don't execute if the script is imported as a module)
