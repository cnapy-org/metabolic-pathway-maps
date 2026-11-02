"""escher_handlayout - toolkit for hand-specified (coordinate-by-coordinate) layouts of Escher maps.

The toolkit only does *mechanical* work: it places markers and cofactor dots relative to
primary-metabolite coordinates written by hand in a spec file, computes Bezier control points,
places labels, packs detached pieces, runs QA checks and writes valid Escher JSON.
It contains no automatic graph-layout algorithm on purpose."""
__version__ = "1.0.0"
from .engine import Layout
