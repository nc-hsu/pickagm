# How a selection works

This page walks through a record selection in the order the steps happen,
and names the pickagm function for each step. The {doc}`example` puts all
of them together in one script.

Overview:

1. **Hazard input**: a disaggregation table and one hazard level.
2. **Rupture contexts** (`pickagm.contexts`): turn the table into the inputs
   that ground-motion models need.
3. **Targets** (`pickagm.distributions`): build the target distributions
   and draw random realisations from them.
4. **Record database** (`pickagm.eqdbases`): load the records and keep the
   intensity measures you need.
5. **Preselection and scaling** (`pickagm.selection`): remove unsuitable
   records and scale the rest to the hazard level.
6. **Scoring and picking** (`pickagm.selection`): compare every record with
   every realisation and pick the best match for each.
7. **Ensemble checks** (`pickagm.distributions`): test how well the picked
   set matches the targets.

## 1. Hazard input

pickagm starts from the output of a probabilistic seismic hazard analysis
(for example from OpenQuake):

- the **conditioning intensity measure** (IM) and its **level**: for
  example SA(0.35 s) = 0.3 g. This is called `conditioning_imt` and `imtl`
  (or `conditioning_iml`) in the code;
- a **disaggregation table** (a pandas DataFrame, `disagg_dst`) that says
  which earthquakes cause that level of shaking. It needs the columns:

  | Column | Meaning |
  |---|---|
  | `TRT` | tectonic region type, e.g. `"Active Shallow Crust"` |
  | `Mag` | magnitude |
  | `Dist` | distance (used as the Joyner–Boore distance `rjb`) |
  | `Eps` | epsilon bin |
  | `P(m\|X=x)` | probability of the scenario given that the level *occurs* |
  | `P(m\|X>x)` | probability of the scenario given that the level is *exceeded* |

  `P(m|X=x)` is always needed, because rows where it is zero are dropped.
  `P(m|X>x)` is only read when you choose exceedance weights (see the
  `occurence` argument below).

The `Eps` column is grouped away: pickagm sums the rows that share the same
`TRT`, `Mag` and `Dist`, and works out epsilon itself from the
ground-motion model and the hazard level.

## 2. Rupture contexts: `pickagm.contexts`

A ground-motion model (GMM) needs the magnitude, distance, site conditions
and so on of each scenario. In OpenQuake these are passed as a NumPy record
array called a *context*.

- `get_rup_ctxs(disagg_type, rup_samples)` turns each row of the
  disaggregation table into a dictionary. It renames the columns:
  `Mag` → `mag`, `Dist` → `rjb`, `Lon` → `lon`, `Lat` → `lat`,
  `TRT` → `trt`, `Eps` → `eps`.
- `get_site_rup_ctxs(disagg_type, disagg_dist, ctx_builder)` passes each
  dictionary to your **context builder** and joins the results into one
  record array.

The context builder is an object you write. It must have a method
`ctx(rup_params)` that takes one of those dictionaries and returns a NumPy
record array with every parameter your GMMs need (site parameters such as
`vs30` included) **and** a `trt` field. Return one row per scenario. The
expected shape of the object is described by the `SiteRupContextBuilder`
protocol in `pickagm.contexts`.

## 3. Targets: `pickagm.distributions`

### Conditional distributions

For each scenario $k$ and each selection IM $j$, the GMM gives the mean
$\mu_{j,k}$ and standard deviation $\sigma_{j,k}$ of $\ln IM_j$. First,
pickagm works out how far the hazard level is above the GMM median for the
conditioning IM $c$, in standard deviations:

$$
\varepsilon_k = \frac{\ln x_c - \mu_{c,k}}{\sigma_{c,k}}
$$

Knowing that the conditioning IM is high (or low) tells us something about
the other IMs, through their correlation $\rho_{j,c}$ with it. The
conditional mean and standard deviation of $\ln IM_j$ are:

$$
\mu_{j,k \mid c} = \mu_{j,k} + \rho_{j,c}\,\sigma_{j,k}\,\varepsilon_k,
\qquad
\sigma_{j,k \mid c} = \sigma_{j,k}\sqrt{1 - \rho_{j,c}^2}
$$

This is done by `conditional_im_distributions`. The correlations come from
a correlation matrix per tectonic region type (see {doc}`conventions` and
{doc}`correlation-models`).

### Target distributions: `gcim_distributions`

Each scenario gives a lognormal distribution for each IM. The target is the
mix of these distributions, weighted by the scenario probabilities $w_k$
(`P(m|X=x)` if `occurence=True`, `P(m|X>x)` if `occurence=False`):

$$
F_j(x) = \sum_k w_k \, \Phi\!\left(\frac{\ln x - \mu_{j,k \mid c}}{\sigma_{j,k \mid c}}\right)
$$

where $\Phi$ is the standard normal cumulative distribution function.

```python
df_stats, pdfs, cdfs = gcim_distributions(
    disagg_dst, imtl, gmm_map, selection_imts, conditioning_imt,
    corr_model_map, ctx_builder, occurence=True)
```

- `df_stats`: one row per selection IM, with the mean, standard deviation
  (`sigma`) and percentiles (by default `p16`, `p50`, `p84`) of $\ln IM$.
- `pdfs`, `cdfs`: dictionaries keyed by IM name. Each value is a two-column
  array: IM values (in linear units) and the probability density or
  cumulative probability.

### Random realisations: `gcim_simulation`

A selection needs one target per record. `gcim_simulation` draws
`n_samples` scenarios from the disaggregation table (with replacement,
weighted by their probability), then draws one random set of IM values for
each scenario. The values for the different IMs are drawn together, so
that they keep the right correlation with each other given the
conditioning IM:

$$
\rho_{ij \mid c} = \frac{\rho_{ij} - \rho_{ic}\,\rho_{jc}}
{\sqrt{\left(1 - \rho_{ic}^2\right)\left(1 - \rho_{jc}^2\right)}}
$$

(`corrmodels.conditional_correlation_matrix`). Each realisation is
$\exp(\mu_{\mid c} + \sigma_{\mid c}\, L u)$, where $u$ is a vector of
standard normal random numbers and $L$ the Cholesky factor of the
conditional correlation matrix (`generate_simulation`).

```python
sims, mu_cond, sig_cond, ctxs, eps = gcim_simulation(
    n_samples, disagg_dst, imtl, gmm_map, selection_imts, conditioning_imt,
    corr_model_map, ctx_builder, occurence=True,
    rup_rng=np.random.default_rng(1), sim_rng=np.random.default_rng(2))
```

- `sims`: DataFrame, one row per selection IM (labelled by IM name) and one
  column per realisation (`0 … n_samples-1`), in linear units.
- `mu_cond`, `sig_cond`: arrays (IMs × realisations) of the conditional
  mean and standard deviation of $\ln IM$.
- `ctxs`: the record array of the sampled scenarios (`ctxs.trt` gives the
  tectonic region type of each realisation).

Two separate random number generators are used: `rup_rng` for sampling the
scenarios and `sim_rng` for the IM values.

## 4. Record database: `pickagm.eqdbases`

- `load_esm_flatfile(n_rows)` reads the bundled ESM flatfile (the first
  `n_rows` rows; `None` reads all of them).
- `extract_im_data(df, idx)` keeps the columns from position `idx` onwards
  (the IM columns); `extract_metadata(df, idx)` keeps the columns before it.
- `filter_imt_data_esm(im_df, imts)` turns the ESM IM columns into the
  form pickagm uses: it renames them to OpenQuake IM names, keeps only the
  IMs in `imts`, interpolates spectral accelerations at periods the
  flatfile does not have, and converts units. The result has two column
  levels: the component (`U`, `V`, `W`) and the IM name. In 0.1.0 this
  works for SA and CAV only; see {doc}`conventions` and
  {doc}`limitations`.

pickagm does not build the full record database for you. For the
selection functions, the database is a DataFrame with **one row per record
component** and two column levels:

- `"metadata"`: at least `mag`, `rjb`, `vs30` (for preselection),
  `event_id`, `station_code`, `location_code` (to avoid picking two
  components of the same recording), and `trt` (if you match tectonic
  region types). `max_usable_T` is needed only if you filter on it.
- `"ims"`: one column per IM, named as in OpenQuake (`"SA(0.1)"`, …),
  including the conditioning IM.

The {doc}`example` shows one way to build it from the ESM flatfile.

## 5. Preselection and scaling: `pickagm.selection`

- `preselection(df, magnitude_bounds, distance_bounds, vs30_bounds, usable_T)`
  keeps records whose magnitude, distance and $V_{S30}$ lie within the
  bounds, and (if `usable_T` is given) whose `max_usable_T` is at least
  `usable_T`. Each bound is a `(lower, upper)` pair; use `None` for no
  limit on one side, or pass `None` instead of the pair to skip that
  filter. The same filters exist one at a time: `filter_on_m`,
  `filter_on_d`, `filter_on_vs30`, `filter_by_usable_T`.
- `scale_records_in_db(df, conditioning_imt, conditioning_iml)` scales every
  record so that its conditioning IM equals the hazard level. It adds the
  scale factor as `("metadata", "alpha")` and the scaled IMs as a new
  column group `"ims_scaled"`. How each IM scales is described in
  {doc}`conventions`.
- `filter_on_sf(df, sf_bounds)` then removes records whose scale factor is
  outside `sf_bounds`, e.g. `(0.25, 4.0)`.

## 6. Scoring and picking: `pickagm.selection`

Every record is compared with every realisation. The score is a weighted
sum of squared differences in log space, each divided by the conditional
standard deviation, so that IMs with a wide target spread count less:

$$
\text{SSE}_{r,n} = \sum_j w_j
\left(\frac{\ln s_{j,n} - \ln a_{j,r}}{\sigma_{j,n \mid c}}\right)^2
$$

where $a_{j,r}$ is the scaled IM $j$ of record $r$ and $s_{j,n}$ the value
in realisation $n$. A lower score is a better match.

- `set_weights(weights, n)` returns the IM weights $w_j$, normalised to sum
  to one. `None` gives equal weights.
- `calculate_sse_scores(df, simulations, sig_cond, weights, conditioning_imt)`
  returns a DataFrame of scores (records × realisations). The conditioning
  IM is left out, since after scaling it is the same for every record.
- `pick_records(df, scores, max_count, sim_trts, ok_trt_matches)` goes
  through the realisations in order and, for each, picks the record with
  the lowest score that is still allowed:
  - the same recording (same `event_id`, `station_code` and
    `location_code`), in any component, cannot be picked twice;
  - at most `max_count` records may come from one event;
  - if `sim_trts` (the tectonic region type of each realisation, e.g.
    `list(ctxs.trt)`) and `ok_trt_matches` (a dictionary from a region type
    to the list of record region types allowed for it) are both given,
    only records with an allowed `("metadata", "trt")` are considered.

  If no record is allowed for a realisation, that row of the result is
  empty (all `NaN`).

## 7. Ensemble checks: `pickagm.distributions`

- `ensemble_ks_test(selected["ims_scaled"], cdfs, p_val=0.05)` runs a
  one-sample Kolmogorov–Smirnov (KS) test for each IM, comparing the
  scaled IMs of the picked records with the target CDF from
  `gcim_distributions`. It returns whether all IMs pass, the test results
  per IM, and the list of IMs that fail (p-value below `p_val`). IMs
  without a target CDF, such as the conditioning IM, are skipped.
- `ensemble_r_score(ks_results, weights)` sums the weighted squared KS
  statistics into one number for the whole set (lower is better); the
  code cites Bradley (2012) for this score.

Other helpers for looking at a picked set: `ensemble_ecdfs` (empirical CDFs
for plotting), `ensemble_quantiles`, and `ensemble_ks_bounds` (the KS
acceptance band around each target CDF).

## AvgSA as an intensity measure: `pickagm.avgSA`

Average spectral acceleration (AvgSA) is the geometric mean of SA over a
set of periods $T_1, \dots, T_n$. `indirect_AvgSA_GMPE` wraps one or more
OpenQuake GMMs for SA and behaves like a GMM for AvgSA, so it can be put in
the GMM map under the name `"AvgSA"`. It computes

$$
\ln \text{AvgSA} = \frac{1}{n}\sum_{i=1}^{n} \ln SA(T_i),
\qquad
\sigma_{\ln \text{AvgSA}} = \frac{1}{n}\sqrt{\sum_{i=1}^{n}\sum_{j=1}^{n}
\rho_{ij}\,\sigma_i\,\sigma_j}
$$

following Kohrangi et al. (2018), as noted in the code. You give the
periods (`avg_periods`) and either the total correlation matrix
(`rho_total`) or the between- and within-event matrices (`rho_deltaB`,
`rho_deltaW`) for those periods. A period of 0 is treated as PGA.
