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

#
# preparing for Python 3 port
#from __future__ import division, print_function
#from __future__ import absolute_import
#from __future__ import unicode_literals

__SILENT_START__ = False

from MetaViz import *
import CBConfig
__version__ = CBConfig.current_version()


if not __SILENT_START__:
    print("#################################################################################################")
    print("# CBMViz - ({}) - Visualization Pipeline for Metabolic Networks                              #".format(__version__))
    print("# (c) Nicolas Huber, Brett G. Olivier & Bas Teusink, Vrije Universiteit Amsterdam, 2016-2017    #")
    print("#                                                                                               #")
    print("# CBMViz is developed as part of the BeBasic MetaToolKit Project                                #")
    print("# Generate a structured network map using only a SBML file (Supports SBML L3V1, FBC and Groups) #")
    print("# Graph object and can be exported as SVG, JSON, GML, new SBML with layout plugin as well       #")
    print("# as SBGN. Graph properties are also stored in a SQLite database.                               #")
    print("#                                                                                               #")
    print("# Necessary  dependencies for this pipeline (on Wind64 Anaconda)                                #")
    print("# ==============================================================                                #")
    print("# Install packages: pip install geojson tarjan libsbgnpy tulip-python tulipgui-python           #")
    print("# Install CBMPy using: conda install -c bgoli -c sbmlteam cbmpy                                 #")
    print("#                                                                                               #")
    print("# CBMViz contains code from Mimoza version X.X.X (http://url) paper doi or pubmed id            #")
    print("#                                                                                               #")
    print("# CBMViz is distributed under the GNU GPL v 3.0 licence, see                                    #")
    print("# LICENCE_GPLv3.txt (supplied with this release) for details                                    #")
    print("#################################################################################################")

