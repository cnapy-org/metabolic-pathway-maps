# metabolic-pathway-maps
Visualising genome scale models and analysis

Contact Brett for more information

## Python 3.10 installation

CBMViz requires Python 3.10 or newer.  Install the core dependencies and the
package from the `cbmviz` directory:

```bash
python3.10 -m pip install .
```

The map-rendering pipeline additionally needs the Tulip and SBGN bindings,
which are distributed separately.  Install them with:

```bash
python3.10 -m pip install '.[visualization]'
```

The bundled example scripts in `cbmviz-test/` require their corresponding
model files and optional visualization dependencies.

(C)Brett G. Olivier, Vrije Universteit Amsterdam, Amsterdam, 2016. All rights reserved.
