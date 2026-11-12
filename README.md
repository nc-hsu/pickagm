# pickagm

> **Warning: pickagm is under heavy development and may change at any time.**
> Functions, arguments and results can change between releases without
> notice. If you need exactly the same results later, pin the version,
> for example `pickagm==0.1.0`.

pickagm is a Python library for selecting earthquake ground-motion records
for nonlinear structural analysis. It selects records that are
*hazard-consistent*: the selected set matches the ground shaking that the
seismic hazard at a site predicts, not only a single target spectrum.

It does this with the **generalised conditional intensity measure (GCIM)**
approach. Given the hazard at one intensity level (for example a spectral
acceleration at a chosen period) and the earthquake scenarios that cause
it, pickagm:

1. builds target distributions of other intensity measures (spectral
   accelerations at other periods, significant duration, Arias intensity,
   and so on), using ground-motion models from
   [OpenQuake](https://github.com/gem/oq-engine) and empirical correlations
   between intensity measures;
2. draws random target realisations from those distributions;
3. scales the records in a database to the conditioning intensity level and
   picks, for each realisation, the record that fits it best;
4. checks how well the selected set matches the targets.

The [ESM flatfile](https://doi.org/10.13127/esm/flatfile.1.0) (intensity
measures and metadata of European and Middle Eastern strong-motion records)
is included with the package as the record database.

## Documentation

<https://pickagm.readthedocs.io>

## Installation

pickagm needs **Python 3.12 or later** and **OpenQuake engine 3.26 or
later**.

Install OpenQuake first. OpenQuake depends on GDAL, which PyPI only offers
as source code, so a plain `pip install openquake.engine` usually fails.
OpenQuake publishes pinned requirement files that download ready-built
wheels instead. For example, for OpenQuake 3.26 and Python 3.12 on Linux:

```console
pip install -r https://raw.githubusercontent.com/gem/oq-engine/engine-3.26/requirements-py312-linux64.txt
pip install openquake.engine
```

On Windows use `requirements-py312-win64.txt` instead, and on macOS
`requirements-py312-macos_arm64.txt` (Apple silicon). See the
[OpenQuake installation guide](https://github.com/gem/oq-engine#installation)
for other options.

Then install pickagm:

```console
pip install pickagm
```

The package is about 24 MB, because it includes the ESM flatfile.

## A first look

Load one of the bundled correlation models and the first rows of the
bundled ESM flatfile:

```python
import numpy as np
from pickagm import corrmodels, eqdbases

# Correlation model for active shallow crust, interpolated to three periods
# (seconds). The periods must be a numpy array.
model = corrmodels.CORR_MODELS["clemett_asc"](
    incl_nonSA_imts=["PGA"],
    incl_SA_Ts=np.array([0.1, 0.35, 1.0]),
)
print(model.rho["total"])   # total-residual correlation matrix (pandas DataFrame)

# The first 5 records of the ESM flatfile.
flatfile = eqdbases.load_esm_flatfile(5)
print(flatfile[["event_id", "station_code", "Mw", "JB_dist"]])
```

A full worked example of a record selection is in the
[documentation](https://pickagm.readthedocs.io/en/latest/example.html).

## How to cite

If you use pickagm, please cite the software and the version you used. The
repository's citation metadata is in
[`CITATION.cff`](https://github.com/nc-hsu/pickagm/blob/main/CITATION.cff);
GitHub shows it under "Cite this repository".

If your work uses the bundled data, please also credit the data sources
listed below.

## Licence

Copyright (C) Nicholas Clemett.

pickagm is licensed under **GPL-3.0-or-later**: you can redistribute it
and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version. It is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
[`LICENSE`](https://github.com/nc-hsu/pickagm/blob/main/LICENSE) file for the
full text.

## Data credits

Full details are in
[`DATA_LICENCES.md`](https://github.com/nc-hsu/pickagm/blob/main/DATA_LICENCES.md).

- **ESM flatfile** (`data/ESM_flatfile_SA.csv`, installed): Engineering
  Strong Motion database flatfile, published by INGV, licensed under
  [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Converted from
  semicolon-separated to comma-separated values; the values are unchanged.
  Lanzano et al. (2018), <https://doi.org/10.13127/esm/flatfile.1.0>, and
  Lanzano et al. (2019), *Bulletin of Earthquake Engineering* 17, 561–582,
  <https://doi.org/10.1007/s10518-018-0480-z>.
- **Correlation coefficients** (`im_corr_data/clemett_corr_coeffs.hdf5`,
  installed): the correlation model of Clemett and Gündel (paper in
  preparation), distributed under GPL-3.0-or-later.
- **Baker and Bradley (2017) correlations** (`im_corr_data/BB17.csv`): kept in
  the source repository only, **not installed**. From
  [bakerjw/NGAW2_correlations](https://github.com/bakerjw/NGAW2_correlations),
  MIT licence, Copyright (c) 2016 Jack Baker.
