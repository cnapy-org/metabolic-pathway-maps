# metabolic-pathway-maps
Visualising genome scale models and analysis

Contact Brett for more information

## Installation

CBMViz requires Python 3.10 or newer and the Python bindings of
[Tulip](https://github.com/Tulip-Dev/tulip) 6.1 (`tulip` and `tulipgui`).
Install the package and its other dependencies from the `cbmviz` directory:

```bash
python -m pip install ./cbmviz
```

The Tulip 6.1 bindings are not available on PyPI and have to be built from
source, with the same Python version that runs CBMViz. After building and
installing Tulip, make its Python package importable, for example with a `.pth`
file in the `site-packages` folder of your (virtual) environment. On Windows
with a Tulip build made with MSYS2 (UCRT64), the file also has to put the
folders with the Tulip and MSYS2 DLLs on `PATH`. Tulip registers every `PATH`
folder as a DLL directory and fails on folders that can not be read, so these are
left out:

```
<tulip install dir>\lib\tulip\python
import os; os.environ['PATH'] = ';'.join(p for p in [r'<tulip install dir>\bin', r'<msys64 dir>\ucrt64\bin'] + os.environ['PATH'].split(';') if p and os.access(p, os.R_OK))
```

The bundled example scripts in `cbmviz-test/` require their corresponding
model files.

(C)Brett G. Olivier, Vrije Universteit Amsterdam, Amsterdam, 2016. All rights reserved.
