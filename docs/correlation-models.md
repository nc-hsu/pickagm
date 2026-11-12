# Correlation models

The target distributions depend on how strongly the IMs are correlated with
each other (see {doc}`workflow`). pickagm ships one family of empirical
correlation models, for four European tectonic settings, in
`pickagm.corrmodels`.

## The registry: `CORR_MODELS`

`CORR_MODELS` is a dictionary from a short key to a model class:

| Key | Class | Setting | `str(model)` |
|---|---|---|---|
| `"clemett_asc"` | `ClemettCorrelationModelAsc` | active shallow crust | `Clemett26ASC` |
| `"clemett_sinter"` | `ClemettCorrelationModelSInter` | subduction interface | `Clemett26SInter` |
| `"clemett_sslab"` | `ClemettCorrelationModelSSlab` | subduction slab (in-slab) | `Clemett26SSlab` |
| `"clemett_vrancea"` | `ClemettCorrelationModelVrancea` | Vrancea deep seismicity | `Clemett26Vrancea` |

The models are from Clemett and Gündel (paper in preparation).

## Building a model

Creating a model reads its coefficients from the bundled HDF5 file and
interpolates them to the IMs you ask for:

```python
import numpy as np
from pickagm.corrmodels import CORR_MODELS

model = CORR_MODELS["clemett_asc"](
    incl_nonSA_imts=["PGA", "RSD595"],      # optional
    incl_SA_Ts=np.array([0.1, 0.35, 1.0]),  # SA periods in seconds
)
rho = model.rho["total"]   # pandas DataFrame, labelled by IM name
```

- `incl_SA_Ts`: the spectral periods you need. **Always give it**; a model
  created without it fails. Periods between the stored ones are
  interpolated linearly. Periods above 8 s raise a `ValueError`. Periods
  below the shortest stored period (0.025 s) are interpolated between PGA
  and SA(0.025); PGA is then included in the matrix.
- `incl_nonSA_imts`: a list of non-SA IM names to include, from `"IA"`,
  `"RSD595"`, `"AvgSA([0, 3])"`, `"AvgSA([0, 6])"` and `"PGA"`. A name that
  is not in the model raises a `ValueError`. If you leave it out, the
  matrix contains only SA (plus PGA in the case above).

`model.rho` is a dictionary with one matrix per type of residual. The keys
follow the notation of Al Atik et al. (2010), as noted in the code:

| Key | Residual |
|---|---|
| `"total"` | total residual |
| `"dB_e"` | between-event residual |
| `"dS2S_s"` | between-site (site-to-site) residual |
| `"dWS_es"` | within-event, within-site residual |

Every matrix has ones on the diagonal. Its rows and columns hold the
requested non-SA IMs first, in the order they are stored in the file
(`IA`, `RSD595`, `AvgSA([0, 3])`, `AvgSA([0, 6])`, `PGA`), then the SA
periods in the order you gave them. SA labels use the period rounded to
three decimals, e.g. `"SA(0.35)"`.

To use a model in a selection, put one matrix per tectonic region type in
a correlation map (see {doc}`conventions`).

## What the HDF5 file holds

The coefficients are in `pickagm/im_corr_data/clemett_corr_coeffs.hdf5`.
The file has one group per setting (`asc`, `sinter`, `sslab`,
`vrancea`). Each group has four datasets, one per residual type (`total`,
`dB_e`, `dS2S_s`, `dWS_es`). Each dataset is a 38 × 38 matrix of
correlation coefficients, with an attribute `imts` (a JSON list) that names
the IM of each row and column, as `[name, period]` pairs:

- 5 non-SA IMs: `IA`, `RSD595`, `AvgSA([0, 3])`, `AvgSA([0, 6])`, `PGA`;
- 33 spectral periods from 0.025 s to 8 s.

`corrmodels.load_corr_coeff_data_from_hdf5(path, group)` reads one group
into a dictionary of NumPy arrays plus the IM list.

:::{note}
The file in version 0.1.0 holds the original coefficients of the model. A
later release will replace them. If you need the 0.1.0 values, pin
`pickagm==0.1.0`.
:::

## Conditional correlations

`conditional_correlation_matrix(corr_df, conditioning_imt)` takes a matrix
and the name of the conditioning IM and returns the correlations of the
other IMs given that IM (formula in {doc}`workflow`). `gcim_simulation`
uses it internally.
