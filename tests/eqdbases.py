from openquake.hazardlib.imt import SA, RSD575, RSD595, PGA, PGV
import numpy as np

import pickagm.eqdbases as eqdb

independent_periods = np.round(np.linspace(0.1, 2.0, 20),2).tolist()
conditioning_period = 0.35
imts = [SA(t) for t in independent_periods] + [SA(conditioning_period)]
imts = sorted(imts, key=lambda x: x.period)


df = eqdb.load_esm_flatfile(20)
df = eqdb.filter_imt_data_esm(df, imts)


pass