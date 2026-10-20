"""
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

"""

from setuptools import setup

mydata_files = []

STATUS = 'alpha'

mypackages= ['cbmviz', 'cbmviz.mod_sbml', 'cbmviz.mod_sbml.annotation','cbmviz.mod_sbml.annotation.chebi',
             'cbmviz.mod_sbml.annotation.gene_ontology','cbmviz.mod_sbml.annotation.kegg','cbmviz.mod_sbml.annotation.pts',
             'cbmviz.mod_sbml.data','cbmviz.mod_sbml.onto','cbmviz.mod_sbml.sbml','cbmviz.mod_sbml.sbml.compartment',
             'cbmviz.mod_sbml.serialization','cbmviz.mod_sbml.utils']

setup(
 package_dir = {'cbmviz' : 'src'},
 packages = mypackages,
 data_files = mydata_files,
 name = "cbmviz",
 description = "CBMViz: Visualising Genome Scale Metabolic Pathways",
 #version = "0.1.0.{}".format(STATUS),
 version = "0.1.1",
 maintainer = 'Brett G. Olivier',
 author = 'Nicolas Huber',
 maintainer_email = 'bgoli@users.sourceforge.net',
 url = "https://github.com/SystemsBioinformatics",
 #download_url = "http://sourceforge.net/projects/cbmpy/files/release/0.7/",
 license = "GNU General Public License (GPL)",
 keywords = "computational systems biology, visualisation, genome scale models, sbml, constraint-based modelling, fbc, linear programming, groups" ,
 python_requires = '>=3.10',
 install_requires = [
     'cbmpy',
     'geojson',
     'jinja2',
     # libsbgnpy 0.4 and newer have a different API
     'libsbgnpy>=0.3.3,<0.4',
     'numpy',
     'openpyxl',
     'pandas',
     'python-libsbml',
     'swiglpk',
     'sympy',
     'tarjan',
 ],
 # the Tulip 6.1 Python bindings (tulip, tulipgui) are not on PyPI and have to be built from source,
 # see the README in the repository root
 package_data = {'cbmviz.mod_sbml': ['data/*']},
 platforms = ["Windows", "Linux"],
 classifiers = [
 'Development Status :: 4 - Beta',
 'Environment :: Console',
 'Intended Audience :: End Users/Desktop',
 'Intended Audience :: Science/Research',
 'License :: OSI Approved :: GNU General Public License (GPL)',
 'Natural Language :: English',
 'Operating System :: OS Independent',
 'Programming Language :: Python :: 3',
 'Pr