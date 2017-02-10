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
    print("########################################################################################################################")
    print("#  CBMViz - ({}) - Visualization Pipeline for Metabolic Networks                                                  #".format(__version__))
    print("#  Nicolas Huber, Dr. Brett Olivier, Prof. Bas Teusink - VU Amsterdam, 2016 - 2017                                     #")
    print("#  Generate a structured network map using only a SBML file (idealy with defined groups) as input. Resulting graph is  #")
    print("#  initially created as Tulip graph object and can be exported as SVG, JSON, GML, new SBML with layout plugin as well  #")
    print("#  as SBGN. Graph properties are also stored in a SQLite database.                                                     #")
    print("#  Necessary  dependencies for this pipeline (on Wind64 Anaconda):                                                     #")
    print("#  pip install cbmpy geojson tarjan libsbgnpy jinja2 tulip-python tulipgui-python                                             #")
    print("#  CBMViz uses parts of the Mimoza software version X.X.X (http://url) paper doi or pubmed id                          #")
    print("########################################################################################################################")

#
