# coding: utf-8
import cbmpy
cmod = cbmpy.readSBML3FBC('iJO1366.xml')
cmod.getReaction('R_THRt2rpp')
cmod.getReaction('R_THRt2rpp').getStoichiometry()
cmod.getReaction('R_THRt2pp').getStoichiometry()
cbmpy.doFBA(cmod)
cmod.getReaction('R_THRt2pp').getValue()
cmod.getReaction('R_URIt2pp_copy1').getValue()
cmod.getReaction('R_URIt2pp_copy2').getValue()
cmod.getReaction('R_URIt2pp_copy2').getStoichiometry()
cmod.getReaction('R_URIt2pp_copy1').getStoichiometry()
get_ipython().magic(u'pinfo cmod.setAllFluxBounds')
cmod.getReaction('R_APH80').getStoichiometry()
cmod.getReaction('R_APH180').getStoichiometry()
cmod.getSpecies('M_ocdcap_c')
cmod.getSpecies('M_ocdcap_c').isReagentOf()
cmod.getSpecies('M_ocdca_c').isReagentOf()
cmod.getReaction('R_AGPAT180').getStoichiometry()
cmod.getSpecies('M_ocdcaACP_c').isReagentOf()
cmod.changeAllFluxBoundsWithValue(-1000, cbmpy.NINF)
cmod.changeAllFluxBoundsWithValue(1000, cbmpy.INF)
cmod.changeAllFluxBoundsWithValue(-999999, cbmpy.NINF)
cmod.changeAllFluxBoundsWithValue(999999, cbmpy.INF)
a = [-1,2,0,1000,1333,-1112]
a
max(abs(a))
import numpy
n = numpy.array(a)
n
man(n)
max(n)
max(abs(n))
a.append(-444000)
max(abs(n))
n = numpy.array(a)
max(abs(n))
b = max(abs(n))
a += [b,-b]
a
n
mina = min(n)
maxa = max(n)
mina
maxa
if abs(min(a)) > abs(max(a)):
    a.append(abs(min(a)))
else:
    a.append(-max(a
    ))
    
a
a = [-1,2,0,1000,1333,-1112]
if abs(min(a)) > abs(max(a)):
    a.append(abs(min(a)))
else:
    a.append(-max(a
    ))
    
a
