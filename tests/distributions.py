import pickagm.distributions as imd
from openquake.hazardlib.imt import SA, RSD575, RSD595, PGA, PGV
from openquake.hazardlib.gsim.boore_1997 import (
    BooreEtAl1997ArbitraryHorizontalUnspecified
    )
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Selection Parameters
mag = 5.38           # mean magnitude contributing to hazard
rjb = 14.9          # mean distance contributing to hazard
eps = 1.13           # mean epsilon contributing to hazard
vs30 = 800.0

indpendent_periods = np.round(np.linspace(0.1, 2.0, 20),2).tolist()
        #   np.linspace(1.0, 2.0, 0).tolist() + \
        #   np.linspace(2, 5, 7).tolist() + \
        #   np.linspace(6, 10, 5).tolist()

conditioning_period = 0.35
conditioning_imt = SA(conditioning_period)

gmm_mapping = {"SA": BooreEtAl1997ArbitraryHorizontalUnspecified()}

corr_model = imd.CorrelationModel("BB17")
    # "BB17", exclude_SA_T_gt=SA(2.0), exclude_SA_T_lt=SA(0.1), 
    # exclude_imts=[RSD575(), RSD595(), PGA(), PGV()])

site_rup_ctx = np.recarray(1, dtype=[('mag', 'f4'), ('rjb', 'f4'), ('vs30', 'f4')])
site_rup_ctx.mag = np.array([mag], dtype='f4')
site_rup_ctx.rjb = np.array([rjb], dtype='f4')  
site_rup_ctx.vs30 = np.array([vs30], dtype='f4') 

# selection stuff
imts = [SA(t) for t in indpendent_periods] + [SA(conditioning_period)]
imts = sorted(imts, key=lambda x: x.period)
periods = np.array([imt.period for imt in imts])

# create the im distributions
im_cond_distribution, imt_distributions = imd.create_im_distributions(
    imts, corr_model, gmm_mapping, site_rup_ctx, conditioning_imt, eps
)

corr_mat, corr_vec = imd.correlation_matrix(
    imts, corr_model, conditioning_imt, incl_conditioning_imt=True)
Sigma_0 = imd.unconditional_covariance_matrix(imt_distributions, corr_mat)
Sigma_1 = imd.covariance_vector(
    imt_distributions, im_cond_distribution, corr_vec)
Sigma = imd.conditional_covariance_matrix(
    Sigma_0, Sigma_1, im_cond_distribution)
Mu = imd.conditional_mean_vector(imt_distributions)

cs = imd.get_cs(imt_distributions, im_cond_distribution)

realisations = imd.sample_im_distribution(Mu, Sigma, n=80)

# n = 2
# fig, ax = plt.subplots()
# for ii in range(r.shape[1]):
#     ax.loglog(periods, np.exp(r[:,ii]), color="0.8")
# ax.loglog(cs[:,0], np.exp(cs[:,1]), ls="-", color="b")
# ax.loglog(cs[:,0], np.exp(cs[:,1] + n * cs[:,2]), ls="-.", color="b")
# ax.loglog(cs[:,0], np.exp(cs[:,1] - n * cs[:,2]), ls="-.", color="b")
# plt.show()

ecdfs = imd.get_empirical_cdfs(realisations)
realisations["mu_samples"] = realisations.mean(axis=1) 
realisations["sigma_samples"] = realisations.std(axis=1) 

target_sigma = imd.conditional_sigma_vector(imt_distributions)

fig, ax = plt.subplots()
ax.loglog(periods, np.exp(Mu), ls="-", color="b")
ax.loglog(periods, np.exp(realisations["mu_samples"]), ls="-.", color="r")

fig, ax = plt.subplots()
ax.semilogx(periods, target_sigma, ls="-", color="b")
ax.semilogx(periods, realisations["sigma_samples"], ls="-.", color="r")

imtd = imt_distributions[5]

print(imtd)
cdf = imtd.cdf(3, cond=True)
ecdf = ecdfs.iloc[5].cdf
fig, ax = plt.subplots()
ax.plot(np.exp(cdf[:,0]), cdf[:,1], ls="-", color="b")
ax.plot(np.exp(ecdf.quantiles), ecdf.probabilities, ls="-", color="r")

plt.show()

pass