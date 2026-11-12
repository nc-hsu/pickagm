# Worked example

This example selects 20 records for one earthquake scenario. It is small on
purpose, so that each step is easy to follow; a real selection would use a
full disaggregation from a hazard analysis.

- Ground-motion model: Boore et al. (1997), as implemented in OpenQuake
  (`BooreEtAl1997ArbitraryHorizontalUnspecified`).
- Conditioning IM: SA(0.35 s), at a hazard level of 0.3 g.
- Selection IMs: SA at 20 periods from 0.1 s to 2.0 s.
- Correlation model: `clemett_asc` (active shallow crust).
- Record database: the bundled ESM flatfile.

The code below is shown, not run, when this documentation is built. It was
run from start to end with pickagm 0.1.0 and OpenQuake 3.26.2 on Python
3.12; the printed results shown are from that run. Random numbers, and
therefore the picked records, can differ with other package versions.
Paste the blocks in order into one script or notebook.

## 1. Imports and scenario

The disaggregation table has a single scenario: magnitude 5.38 at 14.9 km,
with all of the probability.

```python
import numpy as np
import pandas as pd
from openquake.hazardlib.imt import SA
from openquake.hazardlib.gsim.boore_1997 import (
    BooreEtAl1997ArbitraryHorizontalUnspecified,
)

from pickagm import corrmodels, distributions, eqdbases, selection

TRT = "Active Shallow Crust"
conditioning_imt = SA(0.35)
periods = np.round(np.linspace(0.1, 2.0, 20), 2)
selection_imts = [SA(t) for t in periods]
imtl = 0.3   # hazard level of SA(0.35 s), in g

disagg = pd.DataFrame({
    "TRT": [TRT], "Mag": [5.38], "Dist": [14.9], "Eps": [1.13],
    "P(m|X=x)": [1.0], "P(m|X>x)": [1.0],
})
```

## 2. Context builder, GMM map and correlation map

The context builder adds the site's $V_{S30}$ to each scenario and returns
the record array that the Boore et al. (1997) model needs (`mag`, `rjb`,
`vs30`), plus the `trt` field that pickagm uses.

```python
class SimpleContextBuilder:
    def __init__(self, n_ctxs, site_params):
        self.n_ctxs = n_ctxs
        self.site_params = site_params

    def ctx(self, rup_params):
        ctx = np.recarray(
            self.n_ctxs,
            dtype=[("mag", "f8"), ("rjb", "f8"), ("vs30", "f8"),
                   ("eps", "f8"), ("trt", "U32")],
        )
        ctx.mag = rup_params["mag"]
        ctx.rjb = rup_params["rjb"]
        ctx.eps = rup_params["eps"]
        ctx.trt = rup_params["trt"]
        ctx.vs30 = self.site_params["vs30"]
        return ctx


ctx_builder = SimpleContextBuilder(1, {"vs30": 800.0})

gmm = BooreEtAl1997ArbitraryHorizontalUnspecified()
gmm_map = {TRT: {"SA": gmm}}

# The correlation matrix must contain the conditioning period as well.
all_periods = np.sort(np.append(periods, conditioning_imt.period))
corr_model = corrmodels.CORR_MODELS["clemett_asc"](incl_SA_Ts=all_periods)
corr_map = {TRT: corr_model.rho["total"]}
```

## 3. Target distributions and realisations

```python
df_stats, pdfs, cdfs = distributions.gcim_distributions(
    disagg, imtl, gmm_map, selection_imts, conditioning_imt, corr_map,
    ctx_builder, occurence=True)
print(df_stats.head())

n_records = 20
sims, mu_cond, sig_cond, ctxs, eps = distributions.gcim_simulation(
    n_records, disagg, imtl, gmm_map, selection_imts, conditioning_imt,
    corr_map, ctx_builder, occurence=True,
    rup_rng=np.random.default_rng(1), sim_rng=np.random.default_rng(2))
print(sims.shape, mu_cond.shape, sig_cond.shape)
```

Output (the statistics are of $\ln IM$, with IM in g):

```text
             mean     sigma       p16       p50       p84
SA(0.1) -1.394871  0.283798 -1.677097 -1.394871 -1.112645
SA(0.2) -1.007867  0.213382 -1.220066 -1.007867 -0.795668
SA(0.3) -1.105022  0.111710 -1.216112 -1.105022 -0.993931
SA(0.4) -1.354808  0.106406 -1.460625 -1.354808 -1.248991
SA(0.5) -1.652349  0.191506 -1.842794 -1.652349 -1.461904
(20, 20) (20, 20) (20, 20)
```

The spread (`sigma`) is smallest next to the conditioning period of 0.35 s,
because those periods are most strongly correlated with it. The epsilon of
the scenario, worked out from the GMM and the hazard level, is about 1.74
(it is in `eps`).

## 4. Record database

`filter_imt_data_esm` needs only the IM columns, with no missing values
(see {doc}`limitations`). The selection functions then need one row per
record component, with a `"metadata"` and an `"ims"` column group (see
{doc}`workflow`). Here the two horizontal components, `U` and `V`, are used
as separate records.

```python
flatfile = eqdbases.load_esm_flatfile(None)
idx_im_data = flatfile.columns.get_loc("U_pga")
im_df = eqdbases.extract_im_data(flatfile, idx_im_data)
im_df = im_df[[c for c in im_df.columns if c[:2] in ("U_", "V_", "W_")]]
im_df = im_df.dropna()
im_df = eqdbases.filter_imt_data_esm(im_df, selection_imts + [conditioning_imt])

meta = pd.DataFrame({
    "event_id": flatfile["event_id"],
    "station_code": flatfile["station_code"],
    "location_code": flatfile["location_code"],
    "mag": flatfile["Mw"],
    "rjb": flatfile["JB_dist"],
    "vs30": flatfile["vs30_m_sec"],
    "trt": TRT,
})
rows = []
for comp in ["U", "V"]:
    part = pd.concat({"metadata": meta.assign(component=comp),
                      "ims": im_df[comp]}, axis=1)
    part.index = [f"{i}_{comp}" for i in part.index]
    rows.append(part)
db = pd.concat(rows)
db = db[db["ims"].notna().all(axis=1) & (db["ims"] > 0).all(axis=1)]
print(db.shape)
```

```text
(46026, 29)
```

Every record is labelled with the same tectonic region type here, to keep
the example short. In a real selection you would assign each record its own
region type.

## 5. Preselection and scaling

```python
db = selection.preselection(db, magnitude_bounds=(4.5, 6.5),
                            distance_bounds=(None, 50.0),
                            vs30_bounds=(400.0, None), usable_T=None)
db = selection.scale_records_in_db(db, conditioning_imt, imtl)
db = selection.filter_on_sf(db, (0.25, 4.0))
print(db.shape)
```

```text
(135, 51)
```

135 record components are left. The database now also has the scale
factor `("metadata", "alpha")` and the scaled IMs under `"ims_scaled"`.

## 6. Scoring and picking

`sig_cond` is turned into a DataFrame labelled by IM name, so that its rows
are matched to the record columns by name.

```python
weights = selection.set_weights(None, len(selection_imts))
sig_cond_df = pd.DataFrame(sig_cond, index=[imt.string for imt in selection_imts],
                           columns=sims.columns)
scores = selection.calculate_sse_scores(db, sims, sig_cond_df, weights,
                                        conditioning_imt)
picked = selection.pick_records(db, scores, max_count=3,
                                sim_trts=list(ctxs.trt),
                                ok_trt_matches={TRT: [TRT]})
print(picked[[("metadata", "event_id"), ("metadata", "station_code"),
              ("metadata", "component"), ("metadata", "alpha")]].head())
```

```text
                      metadata
                      event_id station_code component     alpha
22045_V           TK-1999-0294         4101         V  1.374611
15608_U           IT-2009-0009          GSA         U  0.801855
22798_V           IR-2002-0035          DSJ         V  2.048045
13873_U           GR-1993-0027         AIGA         U  2.890963
11310_V  EMSC-20161030_0000029          AMT         V  0.280344
```

`max_count=3` allows at most three records from the same earthquake.

## 7. Ensemble checks

```python
passed, ks_results, failed = distributions.ensemble_ks_test(
    picked["ims_scaled"], cdfs, p_val=0.05)
r = distributions.ensemble_r_score(ks_results, weights)
print(passed, failed, r)
```

```text
False ['SA(0.3)', 'SA(0.4)'] 0.0434190065773271
```

With only 135 candidates, the picked set does not pass the KS test at
SA(0.3 s) and SA(0.4 s), the periods next to the conditioning period where
the target is narrowest. Wider preselection bounds give more candidates to
choose from; a real selection would typically try several sets and keep
the best.
