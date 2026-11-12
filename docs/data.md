# Bundled data

The package installs two data files. A third file is kept in the source
repository only. The full credits and licence texts are in
[`DATA_LICENCES.md`](https://github.com/nc-hsu/pickagm/blob/main/DATA_LICENCES.md).

| File (inside the `pickagm` package) | Installed? | Licence | Read by |
|---|---|---|---|
| `data/ESM_flatfile_SA.csv` | yes | CC BY 4.0 | `eqdbases.load_esm_flatfile` |
| `im_corr_data/clemett_corr_coeffs.hdf5` | yes | GPL-3.0-or-later | the `ClemettCorrelationModel*` classes |
| `im_corr_data/BB17.csv` | no | MIT | nothing (repository only) |

## ESM flatfile

The Engineering Strong Motion (ESM) flatfile is a table of strong-motion
recordings from Europe and the Middle East, published by the Istituto
Nazionale di Geofisica e Vulcanologia (INGV) in 2018. Each row is one
recording. The columns give the earthquake (for example `event_id`, `Mw`,
`ev_depth_km`), the station (`station_code`, `vs30_m_sec`, …), distances
(`epi_dist`, `JB_dist`, `rup_dist`, …), processing information, and
intensity measures for each component (`U_…`, `V_…`, `W_…`) and for the
rotated horizontal measures (`rotD50_…`, `rotD100_…`, `rotD00_…`).

The bundled copy has 23,014 records and 329 columns. It is the original
flatfile converted from semicolon-separated to comma-separated values; the
values are unchanged. pickagm reads it with:

```python
from pickagm import eqdbases

flatfile = eqdbases.load_esm_flatfile(None)   # all rows; give a number to read fewer
```

The IM columns start at `U_pga`. How their names and units are converted
for selection is described in {doc}`conventions`.

Licence: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). If your
results depend on it, please cite:

- Lanzano, G., Luzi, L., Russo, E., Felicetta, C., D'Amico, M. C.,
  Sgobba, S., and Pacor, F. (2018). *Engineering Strong Motion Database
  (ESM) flatfile* [Data set]. Istituto Nazionale di Geofisica e
  Vulcanologia (INGV). <https://doi.org/10.13127/esm/flatfile.1.0>
- Lanzano, G., Sgobba, S., Luzi, L., Puglia, R., Pacor, F., Felicetta, C.,
  D'Amico, M., Cotton, F., and Bindi, D. (2019). The pan-European
  Engineering Strong Motion (ESM) flatfile: compilation criteria and data
  statistics. *Bulletin of Earthquake Engineering*, 17(2), 561–582.
  <https://doi.org/10.1007/s10518-018-0480-z>

## Correlation coefficients

`clemett_corr_coeffs.hdf5` holds the coefficients of the correlation models
described in {doc}`correlation-models`, from Clemett and Gündel (paper in
preparation). Version 0.1.0 contains the original coefficients of that
model; a later release will replace them. The file is distributed under
the package licence, GPL-3.0-or-later.

## Baker and Bradley (2017) correlations (not installed)

`BB17.csv` holds the correlation coefficients of Baker and Bradley (2017),
from the repository
[bakerjw/NGAW2_correlations](https://github.com/bakerjw/NGAW2_correlations)
(MIT licence, Copyright (c) 2016 Jack Baker). It is kept in the pickagm
source repository, but no module reads it, so it is not installed with the
package.

- Baker, J. W., and Bradley, B. A. (2017). Intensity measure correlations
  observed in the NGA-West2 database, and dependence of correlations on
  rupture and site parameters. *Earthquake Spectra*, 33(1), 145–156.
  <https://doi.org/10.1193/060716EQS095M>
