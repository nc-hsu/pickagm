# Installation

## Requirements

- **Python 3.12** or later.
- **OpenQuake engine 3.26** or later (`openquake.engine` on PyPI), which
  provides `openquake.hazardlib`.
- NumPy, pandas, SciPy and h5py. pip installs these automatically.

## 1. Install OpenQuake

Install OpenQuake **before** pickagm.

OpenQuake depends on GDAL, a library for geographic data. On PyPI, GDAL is
only available as source code, so pip would try to compile it, which
usually fails. OpenQuake therefore publishes requirement files that
download ready-built packages ("wheels") from its own server, with versions
that are known to work together. Pick the file for your operating system;
the examples below are for OpenQuake 3.26 and Python 3.12.

Linux:

```console
pip install -r https://raw.githubusercontent.com/gem/oq-engine/engine-3.26/requirements-py312-linux64.txt
pip install openquake.engine
```

Windows:

```console
pip install -r https://raw.githubusercontent.com/gem/oq-engine/engine-3.26/requirements-py312-win64.txt
pip install openquake.engine
```

macOS (Apple silicon):

```console
pip install -r https://raw.githubusercontent.com/gem/oq-engine/engine-3.26/requirements-py312-macos_arm64.txt
pip install openquake.engine
```

For other set-ups, see the
[OpenQuake installation guide](https://github.com/gem/oq-engine#installation).

:::{note}
OpenQuake's requirement files pin exact versions of many packages
(for example NumPy, pandas and SciPy). A fresh virtual environment avoids
clashes with packages you already have:

```console
python -m venv .venv
```
:::

## 2. Install pickagm

```console
pip install pickagm
```

To reproduce results exactly, pin the version:

```console
pip install pickagm==0.1.0
```

The download is about 24 MB, because the package includes the ESM flatfile
(see {doc}`data`).

## Check the installation

This loads a bundled correlation model and the first five records of the
bundled flatfile. If both work, the data files were installed correctly.

```python
import numpy as np
from pickagm import corrmodels, eqdbases

model = corrmodels.CORR_MODELS["clemett_asc"](incl_SA_Ts=np.array([0.1, 1.0]))
print(model.rho["total"])

print(eqdbases.load_esm_flatfile(5).shape)   # (5, 329)
```

## Optional extras

- `pip install pickagm[test]` installs pytest, to run the smoke test in the
  repository's `tests/` folder.
- `pip install pickagm[docs]` installs Sphinx and the other tools used to
  build this documentation.
