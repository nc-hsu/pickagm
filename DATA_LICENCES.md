# Data credits and licences

The pickagm source code is licensed under GPL-3.0-or-later (see `LICENSE`).
The repository also holds three data files. This file says where each one
comes from, under which licence it is shared, and whether it is installed
with the package.

| File (under `src/pickagm/`) | Installed with the package? | Licence |
|---|---|---|
| `data/ESM_flatfile_SA.csv` | Yes | CC BY 4.0 |
| `im_corr_data/clemett_corr_coeffs.hdf5` | Yes | GPL-3.0-or-later |
| `im_corr_data/BB17.csv` | No (repository only) | MIT |

---

## 1. ESM flatfile: `data/ESM_flatfile_SA.csv`

This is the Engineering Strong Motion (ESM) flatfile published by the
Istituto Nazionale di Geofisica e Vulcanologia (INGV) in 2018, converted
from semicolon-separated to comma-separated values. The values themselves
are unchanged.

- **Licence:** Creative Commons Attribution 4.0 International (CC BY 4.0),
  <https://creativecommons.org/licenses/by/4.0/>.
- **Dataset:** Lanzano, G., Luzi, L., Russo, E., Felicetta, C., D'Amico,
  M. C., Sgobba, S., and Pacor, F. (2018). *Engineering Strong Motion
  Database (ESM) flatfile* [Data set]. Istituto Nazionale di Geofisica e
  Vulcanologia (INGV). <https://doi.org/10.13127/esm/flatfile.1.0>
- **Paper describing it:** Lanzano, G., Sgobba, S., Luzi, L., Puglia, R.,
  Pacor, F., Felicetta, C., D'Amico, M., Cotton, F., and Bindi, D. (2019).
  The pan-European Engineering Strong Motion (ESM) flatfile: compilation
  criteria and data statistics. *Bulletin of Earthquake Engineering*,
  17(2), 561–582. <https://doi.org/10.1007/s10518-018-0480-z>

If you use results that depend on this file, please cite the dataset and
the paper above.

## 2. Correlation coefficients: `im_corr_data/clemett_corr_coeffs.hdf5`

Empirical correlation coefficients between ground-motion intensity
measures for four European tectonic settings (active shallow crust,
subduction interface, subduction slab and Vrancea), used by the
`ClemettCorrelationModel*` classes in `pickagm.corrmodels`. They come from
the correlation model of Clemett and Gündel (paper in preparation).

This file holds the original coefficients of that model. A later release
of pickagm will replace them.

- **Licence:** distributed under the package licence, GPL-3.0-or-later.

## 3. Baker and Bradley (2017) correlations: `im_corr_data/BB17.csv`

Correlation coefficients from Baker and Bradley (2017), as published in the
GitHub repository <https://github.com/bakerjw/NGAW2_correlations>.

This file is kept in the source repository only. **It is not installed with
the package**, because no pickagm module reads it.

- **Paper:** Baker, J. W., and Bradley, B. A. (2017). Intensity measure
  correlations observed in the NGA-West2 database, and dependence of
  correlations on rupture and site parameters. *Earthquake Spectra*, 33(1),
  145–156. <https://doi.org/10.1193/060716EQS095M>
- **Licence:** MIT. The copyright and permission notice of the source
  repository, copied verbatim:

```text
The MIT License (MIT)

Copyright (c) 2016 Jack Baker

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
