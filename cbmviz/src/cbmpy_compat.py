"""
Workarounds that let cbmpy (0.8.8.1) run with current versions of its dependencies. Each one is only applied if the
problem is detected, so this module does nothing once cbmpy or its dependencies are fixed. cbmpy looks up these
names at call time, so patching the modules after cbmpy has been imported is enough.

1. NumPy 2 removed numpy.Inf, numpy.NINF, numpy.NaN and numpy.alltrue, which cbmpy still uses (e.g. when writing
   SBML or setting flux bounds). They are restored as aliases of numpy.inf, -numpy.inf, numpy.nan and numpy.all.

2. swiglpk.as_intArray / swiglpk.as_doubleArray are broken in swiglpk 5.0.12 and 5.0.13 on Windows with Python 3.13
   (the only versions with wheels for it): the arrays they return do not contain the list values at the 1-based
   positions GLPK reads, so the GLPK solver used by cbmpy aborts the whole process
   ("glp_set_mat_row: ... column index out of range"). They are replaced by pure Python versions that fill the
   arrays element by element.
"""

import numpy

try:
    import swiglpk
except ImportError:
    swiglpk = None

for _name, _value in (('Inf', numpy.inf), ('NINF', -numpy.inf), ('NaN', numpy.nan), ('alltrue', numpy.all)):
    try:
        getattr(numpy, _name)
    except AttributeError:
        setattr(numpy, _name, _value)


def as_intArray(values):
    """Return a GLPK int array with values[i] at position i + 1 (GLPK arrays are 1-based)."""
    arr = swiglpk.intArray(len(values) + 1)
    for i, value in enumerate(values):
        arr[i + 1] = int(value)
    return arr


def as_doubleArray(values):
    """Return a GLPK double array with values[i] at position i + 1 (GLPK arrays are 1-based)."""
    arr = swiglpk.doubleArray(len(values) + 1)
    for i, value in enumerate(values):
        arr[i + 1] = float(value)
    return arr


def _swiglpk_arrays_ok():
    int_arr = swiglpk.as_intArray([7, 8])
    double_arr = swiglpk.as_doubleArray([0.5, 1.5])
    return int_arr[1] == 7 and int_arr[2] == 8 and double_arr[1] == 0.5 and double_arr[2] == 1.5


if swiglpk is not None and not _swiglpk_arrays_ok():
    swiglpk.as_intArray = as_intArray
    swiglpk.as_doubleArray = as_doubleArray
