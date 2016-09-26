.. include:: <isonum.txt>

######################################
Test MIMOZA with different model files
######################################

************************
Ecoli_iAF1260.glc.l3.xml
************************

	|rarr| partial success, id: aaccc933ca825cee2a9e80e08d47c3f2

* after saving the model the process got killed (by kernel? - model too big?):

::

	.
	.
	2016-09-15 13:25:26: tlp nodes -> geojson features
	2016-09-15 13:28:05: serializing layout
	2016-09-15 13:28:09: saving to /home/nicolas/models_to_test/aaccc933ca825cee2a9e80e08d47c3f2/CBMPY_Model_with_groups.xml
	Getötet

* last lines should be:

::
	
	.	
	.
	2016-09-15 13:33:37: create html
	2016-09-15 13:33:37: Rendering the model...
	2016-09-15 13:33:37: Rendering the mini model...	

******************
Hpylori_iIT341.xml
******************

	|rarr| success, id: 4da52117a280175f83efe5be7e5b38c0

*********************
Hpylori_iIT341.v1.xml
*********************

	|rarr| success, id: b7f382545906454ed0033da2e2565380
	
************************	
Hpylori_iIT341.cobra.xml
************************

	|rarr| success, id: b86d3d89c093d2c02221bba44abc0f75

**************
iAF1260.v1.xml
**************

* error:

::

	2016-09-15 13:53:42: Generalizing the model...
	2016-09-15 13:54:35: mapping species to ChEBI
	Traceback (most recent call last):
  File "main.py", line 18, in <module>
    process_sbml(params.model, params.verbose, ub_ch_ids=None, web_page_prefix=None, generalize=True, log_file=None)
  File "/home/nicolas/mimoza-master/sbml_vis/mimoza_pipeline.py", line 96, in process_sbml
    generalize_model(groups_sbml, gen_sbml, sbml, chebi, ub_chebi_ids=ub_ch_ids)
  File "build/bdist.linux-x86_64/egg/sbml_generalization/generalization/sbml_generalizer.py", line 52, in generalize_model
  File "/usr/local/lib/python2.7/dist-packages/mod_sbml-0.1.1-py2.7.egg/mod_sbml/onto/__init__.py", line 121, in filter_ontology
    terms_to_keep |= onto.get_sub_tree(term, relationships=relationships)
  File "/usr/local/lib/python2.7/dist-packages/mod_sbml-0.1.1-py2.7.egg/mod_sbml/onto/obo_ontology.py", line 323, in get_sub_tree
    | self.get_equivalents(t, None, 0, relationships) | {t}
  File "/usr/local/lib/python2.7/dist-packages/mod_sbml-0.1.1-py2.7.egg/mod_sbml/onto/obo_ontology.py", line 338, in get_generalized_descendants
    set())
  File "/usr/local/lib/python2.7/dist-packages/mod_sbml-0.1.1-py2.7.egg/mod_sbml/onto/obo_ontology.py", line 337, in <genexpr>
    set()) for it in terms),
	AttributeError: 'NoneType' object has no attribute 'get_id'
	
***********
iAF1260.xml
***********

* error:

::

  File "/usr/local/lib/python2.7/dist-packages/mod_sbml-0.1.1-py2.7.egg/mod_sbml/onto/obo_ontology.py", line 337, in <genexpr>
    set()) for it in terms),
	AttributeError: 'NoneType' object has no attribute 'get_id'

***********************
Ecoli_iJR904.glc.l3.xml
***********************

	|rarr| success, id: e8c68d01c5a70ab13aca816be10e974f

*****************************
Lactococcus_lactis_MG1363.xml
*****************************

	|rarr| success, id: 3900dbddaadbae9718903d5f85c949d9

*****************
LBulg.lcts.l3.xml
*****************

	|rarr| success, id: e52444079b5e775e82a94cb69af03647

******************
Saureus_iSB619.xml
******************

	|rarr| success, id: c2a961fa1a229a5d01378e30c3130c67

**********************
Scerevisiae_iND750.xml
**********************

	|rarr| success, id: 8c1296044fcc628a533453130719c700

