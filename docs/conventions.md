# Conventions

This page collects the naming rules and data shapes that the functions
expect. Most problems when setting up a selection come from a name that
does not match, so it is worth reading once.

## Intensity measure names

pickagm names intensity measures (IMs) exactly as OpenQuake prints its IM
types (`openquake.hazardlib.imt`). The name is the string you get with
`str(imt)` or `imt.string`:

| OpenQuake object | Name used in pickagm |
|---|---|
| `SA(0.1)` | `"SA(0.1)"` |
| `SA(1.0)` | `"SA(1.0)"` |
| `PGA()` | `"PGA"` |
| `PGV()` | `"PGV"` |
| `RSD595()` | `"RSD595"` |
| `CAV()` | `"CAV"` |
| `IA()` | `"IA"` |
| `AvgSA([0, 3])` | `"AvgSA([0, 3])"` |

These names label the rows and columns of correlation matrices, the
columns of the record database, the index of the realisations from
`gcim_simulation`, and the keys of the target CDFs.

Functions that take IMs as arguments (for example `selection_imts` and
`conditioning_imt`) take OpenQuake IM **objects**, not strings.

The *type* of an IM, without the period, is `imt.name`: `"SA"` for every
spectral acceleration, `"AvgSA"` for every AvgSA. This is the key used in
the GMM map (below).

## Names and units of the ESM flatfile

`eqdbases.filter_imt_data_esm` renames the ESM spectral acceleration
columns to OpenQuake names and converts their units:

| ESM column (per component) | Becomes | Units |
|---|---|---|
| `T0_010`, `T1_000`, … (SA at 0.01 s, 1.0 s, …) | `SA(0.01)`, `SA(1.0)`, … | cm/s² divided by 981, giving g |
| `CAV` | `CAV` | unchanged |

The module also defines the names and units meant for the other ESM IMs:

- `eqdbases.esm_nonSA_im_translations`: `pga` → `PGA`, `pgv` → `PGV`,
  `T90` → `RSD595`, `CAV` → `CAV`, `ia` → `IA`, `pgd` → `PGD`;
- `eqdbases.esm_unit_conversions`: `PGA`, `SA` and `AvgSA` divided by 981
  (cm/s² to g), `IA` divided by 100 (cm/s to m/s).

:::{warning}
In version 0.1.0 these non-SA renames are **not applied** by
`filter_imt_data_esm`: the ESM columns keep their own names (`pga`, `pgv`,
`T90`, `ia`, …), so a request for `PGA()`, `PGV()`, `RSD595()` or `IA()`
finds nothing and the IM is left out without a warning. Only SA and `CAV`
(whose ESM name is already `CAV`) come through. If you need other IMs,
rename and convert those columns yourself before building the record
database. See {doc}`limitations`.
:::

If you ask for a spectral acceleration at a period that the flatfile does
not have, it is interpolated linearly in log(period)–log(SA) between the
neighbouring periods (`dfops.interpolate_for_sa`). The flatfile has periods
from 0.01 s to 10 s; periods outside that range raise a `ValueError`.

## Record components: `U`, `V`, `W`

The ESM flatfile stores each recording's three components under the
prefixes `U`, `V` and `W`. In the bundled flatfile, `W` is always the
vertical channel (channel code ending in `Z`), and `U` and `V` are the two
horizontal channels (mostly east and north). `filter_imt_data_esm` keeps
all three; which ones you use for selection is up to you. The `rotD50`,
`rotD100` and `rotD00` columns of the flatfile are not used by
`filter_imt_data_esm`.

The scaling and scoring functions treat every row of the record database
as one single-component record.

## Tectonic region types

The tectonic region type (TRT) is a plain string, such as
`"Active Shallow Crust"`. The same strings must be used in the
disaggregation table (`TRT` column), in the `trt` field of the contexts,
as keys of the GMM map and the correlation map, and (if you use it) in the
`("metadata", "trt")` column of the record database. pickagm does not
check them against OpenQuake's list.

## `GmmMap` and `CorrMap`

These two type aliases are defined in `pickagm.typing`:

```python
GmmMap = dict[TRT, dict[IMTName, GMPE]]
CorrMap = dict[TRT, pd.DataFrame]
```

**`GmmMap`** says which ground-motion model to use for each tectonic region
type and each IM type. The outer key is the TRT; the inner key is
`imt.name` (`"SA"`, `"PGA"`, `"RSD595"`, `"AvgSA"`, …); the value is an
OpenQuake GMPE object (or `avgSA.indirect_AvgSA_GMPE`). It must contain the
conditioning IM type as well as every selection IM type:

```python
gmm_map = {
    "Active Shallow Crust": {
        "SA": BooreEtAl1997ArbitraryHorizontalUnspecified(),
    },
}
```

**`CorrMap`** gives one correlation matrix per TRT, as a pandas DataFrame
whose index and columns are IM names. It must contain the conditioning IM
and every selection IM. The matrices from the bundled models are in the
`rho` attribute of a correlation model:

```python
model = corrmodels.CORR_MODELS["clemett_asc"](incl_SA_Ts=periods)
corr_map = {"Active Shallow Crust": model.rho["total"]}
```

## Scaling

Scaling a record by an amplitude factor $\alpha$ changes its IMs by
different amounts. pickagm describes this with an exponent $\gamma$ per IM
type (`distributions.SCALING_EXPONENTS`): the IM is multiplied by
$\alpha^{\gamma}$.

| IM type | $\gamma$ | Why |
|---|---|---|
| `SA`, `PGA`, `PGV`, `CAV`, `AvgSA` | 1 | proportional to amplitude |
| `IA` | 2 | proportional to amplitude squared |
| `RSD595` | 0 | a duration: does not change with amplitude |

The factor that brings the conditioning IM $c$ of a record to the hazard
level $x_c$ is

$$
\alpha = \left(\frac{x_c}{IM_c}\right)^{1/\gamma_c}
$$

and every other IM $j$ of the record is then multiplied by
$\alpha^{\gamma_j}$ (`selection.scale_records_in_db`). Because
$\gamma = 0$ for RSD595, a duration cannot be the conditioning IM. An IM
whose name does not start with one of the names in the table raises
`ValueError("Not all IMTs were scaled")`.
