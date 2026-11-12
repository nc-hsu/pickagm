# Known limitations

pickagm 0.1.0 is an early version. These are the limits we know about,
found by reading the code and by running it. They are listed so that you
are not surprised; they may change in later releases.

## Not implemented

These raise `NotImplementedError`:

- **Weighting of GMMs.** `avgSA.indirect_AvgSA_GMPE(..., gmm_weights=...)`
  with any weights. If you pass several GMMs, all of them are evaluated in
  turn into the same output arrays, so only the last one's values are kept;
  pass a single GMM.
- **Weighting of periods in AvgSA.** `avgSA.compute_ln_AvgSA` and
  `avgSA.compute_std_dev_ln_AvgSA` with `weights`. AvgSA always uses equal
  weights.

## Things that fail without a clear message

- **Correlation models need `incl_SA_Ts`.** Creating a
  `ClemettCorrelationModel*` without SA periods fails with a `TypeError` or
  `UnboundLocalError`.
- **`EmpiricalCorrelationModel.get_correlation`** fails with an
  `AttributeError`. Read values from the matrices in `model.rho` instead,
  e.g. `model.rho["total"].loc["SA(0.1)", "SA(1.0)"]`.
- **`filter_imt_data_esm` needs IM columns only, and no missing values.**
  Pass the IM part of the flatfile (from the `U_pga` column on, e.g. with
  `extract_im_data`), not the whole flatfile, and drop records with
  missing values in the components you use. Otherwise it fails with a
  `TypeError` or with `ValueError: Array must not contain infs or nans`.
- **A duration cannot be the conditioning IM.** RSD595 has a scaling
  exponent of 0 (see {doc}`conventions`), so scaling to it divides by
  zero.
- **`avgSA.indirect_AvgSA_GMPE` needs `avg_periods`.** Without them, the
  AvgSA IM itself is passed on to the SA model.

## Things that happen silently

- **IMs missing from the database are skipped.** `filter_imt_data_esm` only
  returns the IMs it can find in the flatfile (or, for SA, interpolate).
  An IM that the flatfile does not have, such as `RSD575`, is left out
  without a warning, as its own docstring says. Check the columns of the
  result before you go on.
- **Only SA and CAV come through `filter_imt_data_esm`.** The renames from
  the ESM names `pga`, `pgv`, `T90` and `ia` to `PGA`, `PGV`, `RSD595` and
  `IA` are defined in `eqdbases` but are not applied in 0.1.0: the code
  passes tuple keys to pandas' `DataFrame.rename`, which pandas ignores for
  multi-level columns. So `PGA()`, `PGV()`, `RSD595()` and `IA()` are
  skipped like a missing IM, and their unit conversions never happen. To
  use them, rename and convert those columns yourself (see
  {doc}`conventions` for the intended names and units).
- **A bound of 0 counts as "no bound".** In `preselection` and the
  `filter_on_*` functions, a lower or upper bound of `0` is treated like
  `None`.
- **Missing metadata removes the record.** A record whose magnitude,
  distance or $V_{S30}$ is missing (`NaN`) fails the bound check on that
  value, so preselection drops it whenever that bound is used.
- **Empty picks.** If no record is allowed for a realisation,
  `pick_records` returns an all-`NaN` row for it rather than raising an
  error.

## Scope

- Only a full disaggregation by tectonic region type, magnitude, distance
  and epsilon is supported by `gcim_distributions` and `gcim_simulation`.
- The disaggregation distance is used as the Joyner–Boore distance `rjb`.
- The bundled database is the 2018 ESM flatfile. Other databases can be
  used, but you have to bring them into the layout described in
  {doc}`workflow`.
- There are no docstrings yet for most functions, so the
  {doc}`api/index` shows signatures only.
