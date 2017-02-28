"""
CBMPy: CBConfig module
======================
PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net)
Copyright (C) 2009-2016 Brett G. Olivier, VU University Amsterdam, Amsterdam, The Netherlands

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

Author: Brett G. Olivier
Contact email: bgoli@users.sourceforge.net
Last edit: $Author: bgoli $ ($Id: CBConfig.py 561 2017-01-24 13:36:09Z bgoli $)

"""
## gets rid of "invalid variable name" info
# pylint: disable=C0103
## gets rid of "line to long" info
# pylint: disable=C0301
## use with caution: gets rid of module xxx has no member errors (run once enabled)
# pylint: disable=E1101


# preparing for Python 3 port
#from __future__ import division, print_function
#from __future__ import absolute_import
#from __future__ import unicode_literals

# release

__CBCONFIG__ = {'VERSION_MAJOR' : 0,
                'VERSION_MINOR' : 1,
                'VERSION_MICRO' : 0,
                'VERSION_STATUS' : 'alpha',
                'VERSION' : '',
                'DEBUG' : False
               }

# default options dictionary controlling CBMViz (group based layout with flux coloring)
cfg_default = {"AUTO_ALIAS": False,
               # apply aliasing to every node of the graph with a certain connectivity
               # and option to alias biomass (result in main graph)
               # if auto alias is active the automatically generated layout will still be
               # overwritten if any other layouting functions are active as well
               "AUTO_ALIAS_OPTIONS": {"spec_split": 6, "reac_split": 4},  # threshold values for species and reac
                                                                          # splitting

               "PATHWAY_LISTS": True,    # layout annotated groups of sbml file as seperate pathways
               "SMALL_THRSH": 30,        # the threshold value of connections a species must have to be a small sps
               "LAYOUT_PATHWAYS": True,  # needs PATHWAY_LISTS true, so that necessary node lists are generated
               "LABEL_TYPE": "id",       # can either be "id" or "name"
               "LABEL_FONT_SIZE": 18,     # default value: 18
               "ADD_COFACTORS": False,
               #"ADD_COFACTORS": ["Respiratory chain",    # add nodes for small sps/cofactors after they have been
               #                  "Citric Acid Cycle"],   # deleted in path_dict

               "TRANSPORT_RADIAL": False,                   # apply Tree Radial algorithm to transport reactions in
                                                            # output_whole
               "METABOLITE_NODE_SIZE": 2.5, # for big nodes: set to 10
               "REACTION_NODE_SIZE": 1, # for big nodes: set to 4
               "COFACTOR_NODE_SIZE": 1.5, # for big nodes: set to 7

               "SELECTED_PATHWAYS": ["objective reaction",  # pathways to be saved in output_selected
                                     "Transport"],
               "GRID_OUTPUT": ["Glycolysis",                # align pathways on a grid, control position
                               "TCA cycle",
                               "Respiratory chain",
                               "Aminosugar Metabolism"],

               "EXPORT_TLP": True,   # save generated graphs as .tlp files
               "EXPORT_GML": True,   # export generated graphs as .gml files
               "EXPORT_SVG": True,   # export generated graphs as .svg files
               "EXPORT_JSON": True,  # export generated graphs as .json files
               "EXPORT_SBML": True,  # export generated graphs as .sbml files (in layout plugin)
               "EXPORT_SBGN": True,  # export generated graphs as .sbgn files

               "COLOR_FBA": True,              # color the reaction nodes and in and outgoing edges
               "COLOR_FBA_COLORBLIND": False,  # based on flux values obtained from FBA
               "COLOR_PATHWAYS": False,        # pathways specific coloring
               "COLOR_GAPS": False             # highlight gaps (in red), reactions that exist in both (all) models =
                                               # black
               }

# automatic layouting options dictionary (if no groups are defined)
cfg_auto = {"AUTO_ALIAS": False,
            # apply aliasing to every node of the graph with a certain connectivity
            # and option to alias biomass (result in main graph)
            # if auto alias is active the automatically generated layout will still be
            # overwritten if any other layouting functions are active as well
            "AUTO_ALIAS_OPTIONS": {"spec_split": 6, "reac_split": 4},  # threshold values for species and reac
                                                                       # splitting

            "PATHWAY_LISTS": True,    # layout annotated groups of sbml file as seperate pathways
            "SMALL_THRSH": 30,        # the threshold value of connections a species must have to be a small sps
            "LAYOUT_PATHWAYS": True,  # needs PATHWAY_LISTS true, so that necessary node lists are generated
            "LABEL_TYPE": "id",       # can either be "id" or "name"
            "LABEL_FONT_SIZE": 80,    # default value: 80
            "ADD_COFACTORS": False,
            #"ADD_COFACTORS": ["Respiratory chain",    # add nodes for small sps/cofactors after they have been
            #                  "Citric Acid Cycle"],   # deleted in path_dict

            "TRANSPORT_RADIAL": True,                    # apply Tree Radial algorithm to transport reactions in
                                                         # output_whole

            "METABOLITE_NODE_SIZE": 2.5, # for big nodes: set to 10
            "REACTION_NODE_SIZE": 1, # for big nodes: set to 4
            "COFACTOR_NODE_SIZE": 1.5, # for big nodes: set to 7

            "SELECTED_PATHWAYS": ["objective reaction",  # pathways to be saved in output_selected
                                  "Transport"],
            "GRID_OUTPUT": ["Glycolysis",                # align pathways on a grid, control position
                            "TCA cycle",
                            "Respiratory chain",
                            "Aminosugar Metabolism"],

            "EXPORT_TLP": True,   # save generated graphs as .tlp files
            "EXPORT_GML": True,   # export generated graphs as .gml files
            "EXPORT_SVG": True,   # export generated graphs as .svg files
            "EXPORT_JSON": True,  # export generated graphs as .json files
            "EXPORT_SBML": True,  # export generated graphs as .sbml files (in layout plugin)
            "EXPORT_SBGN": True,  # export generated graphs as .sbgn files

            "COLOR_FBA": True,              # color the reaction nodes and in and outgoing edges
            "COLOR_FBA_COLORBLIND": False,  # based on flux values obtained from FBA
            "COLOR_PATHWAYS": False,        # pathways specific coloring
            "COLOR_GAPS": False             # highlight gaps (in red), reactions that exist in both (all) models =
                                           # black
            }

def current_version():
    """
    Return the current CBMPy version as a string

    """
    return '{}.{}.{}'.format(__CBCONFIG__['VERSION_MAJOR'], __CBCONFIG__['VERSION_MINOR'], __CBCONFIG__['VERSION_MICRO'])

def current_version_tuple():
    """
    Return the current CBMPy version as a tuple (x, y, z)

    """
    return (__CBCONFIG__['VERSION_MAJOR'], __CBCONFIG__['VERSION_MINOR'], __CBCONFIG__['VERSION_MICRO'])

__CBCONFIG__['VERSION'] = current_version()


