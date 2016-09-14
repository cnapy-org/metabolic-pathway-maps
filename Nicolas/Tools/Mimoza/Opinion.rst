.. include:: <isonum.txt>

###########
Mimoza Test
###########

* lots of dependencies needed to run mimoza locally (from ReadMe):

::

	Dependencies
	==================

	- *Mimoza* uses `libSBML <http://sbml.org/Software/libSBML>`_ library for python with 	 the groups and layout extensions.
	Install it:
	* sudo pip install python-libsbml-experimental

	- *Mimoza* uses *Model Generalization* library for python to produce generalized views of the model.
	To install it, please download and extract the source files from `gforge.inria.fr <http://metamogen.gforge.inria.fr>`_.
	From the directory where you have extracted it, execute:
	* sudo python setup.py
	
	- *Mimoza* uses *Tulip 4.4* library for python to layout metabolic networks.
	To install it, please download and extract the source files from `sourceforge.net <https://sourceforge.net/projects/auber/files/tulip/tulip-4.4.0/>`_.
	From the directory where you have extracted it, execute:
	* mkdir build
	* cd build
	* cmake ..
	* make
	* sudo make install

	- *Mimoza* uses `SymPy <http://www.sympy.org/>`_, `Python bindings for geojson <https://pypi.python.org/pypi/geojson/>`_
	and `Jinja2 <http://jinja.pocoo.org/>`_:
	* sudo pip install sympy 
	* sudo pip install geojson 
	* sudo pip install Jinja2 

	If you want to have export of your maps in `SBGN PD<http://www.sbgn.org/>`_, install `libSBGN bindings for Python <https://github.com/matthiaskoenig/libsbgn-python>`_:
	* sudo pip install libsbgnpy

*******
Errors:
*******

::
  
	C:\Users\Nicolas\Anaconda2\lib\distutils\dist.py:267: UserWarning: Unknown distribution option: 'platform'
	warnings.warn(msg)
	usage: setup.py [global_opts] cmd1 [cmd1_opts] [cmd2 [cmd2_opts] ...]
	or: setup.py --help [cmd1 cmd2 ...]
	or: setup.py --help-commands
	or: setup.py cmd --help

	error: no commands supplied

Update 12/09/2016:
==================

* previous error was because a command was missing |rarr| enter "python setup.py build" in console:

::

	File "C:\Users\Nicolas\Anaconda2\lib\distutils\util.py", line 126, in convert_path
		raise ValueError, "path '%s' cannot end with '/'" % pathname
	ValueError: path 'sbml_generalization/data/' cannot end with '/'

|rarr| mod_gen does not run on windows systems! need virtual machine with linux (e.g. ubuntu) in virtual box to test

* error after compiling tulip and installing all other described dependencies:

::

	File "/home/user/Mimoza/sbml_vis/converter/sbml2tlp.py", line 2, in <module>
		from tulip import tlp
	ImportError: cannot import name tlp

|rarr| answer from author:

* "the import error 'from tulip import tlp' seems to be due the missing package 'tulip-python ' that you can install (at least on linux and mac) via 'sudo pip install tulip-python'"

Update 13/09/2016:
==================

* tulip built and compiled again after "sudo apt-get install graphviz" for doxygen_dot:executable:
	|rarr| same ImportError (tlp from tulip)
	|rarr| sudo apt-get install tulip-python

Update 13/09/2016:
==================

* new error on ubuntu in virtual machine, after fixing tlp import error:

::

	File "/home/nicolas/Mimoza/sbml_vis/converter/sbml2tlp.py", line 5, in <module>
		from sbml_generalization.sbml.compartment.compartment_positioner import get_comp2go, comp2level
	ImportError: No module named compartment.compartment_positioner

|rarr| installed mod_sbml (manually, pip install did not work)

* additionally required libraries:

	* tarjan (pip install tarjan)
	* libsbgnpy (pip install libsbgnpy)
	* lxml (pip install lxml)

|rarr| after allocating 1gb more ram (up to 2gb) to the virtual machine the call "python main.py --model "model_path" --verbose" got no longer killed

* new error:

::

	File "/usr/local/lib/python2.7/dist-packages/SBMLModelGeneralization-1.0-py2.7.egg/sbml_generalization/annotation/annotate_with_chebi.py", line 117, in get_species_to_chebi
		t_ids = chebi.get_ids_by_name(name)
	AttributeError: Ontology instance has no attribute 'get_ids_by_name'

Update 14/09/2016:
==================

* AttributeError from 13/09/2016:

::

	File "/usr/local/lib/python2.7/dist-packages/SBMLModelGeneralization-1.0-py2.7.egg/sbml_generalization/annotation/annotate_with_chebi.py", line 117, in get_species_to_chebi
		t_ids = chebi.get_ids_by_name(name)
	AttributeError: Ontology instance has no attribute 'get_ids_by_name'

|rarr| suggestion by author: new version of model_gen from github
|rarr| success!


*****************
How do they draw?
*****************

* python based software; they use the Tulip python API (Tulip is originally written in C++) 

*******
Opinion
*******

* since I was not able to install the Model Generalization library, I was not able to test Mimoza yet but idea of library sounds promising - look at the code!
 