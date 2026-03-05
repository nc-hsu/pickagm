"""
Module for the different correlation models that can be employed for record selection
"""
import re
import numpy as np
from   numpy.typing import NDArray
from openquake.hazardlib.imt import IMT
from openquake.hazardlib.gsim.base import GMPE
import pandas as pd
import importlib.resources as res
from scipy import stats
from scipy.optimize import brentq
from scipy.stats._distn_infrastructure import rv_frozen

from pickagm.typing import GmmMap, CorrMap, IMTName 
from pickagm.corrmodels import conditional_correlation_matrix
from pickagm.contexts import (
    SiteRupContextBuilder, 
    get_rup_ctxs, 
    get_site_rup_ctxs
    )
from pickagm.dfops import interpolate_for_sa

NormalRV = rv_frozen  # alias for readability

SCALING_EXPONENTS = {
    "SA": 1.0,
    "PGA": 1.0,
    "PGV": 1.0,
    "RSD595": 0,
    "CAV": 1.0,
    "IA": 2.0,
    "AvgSA": 1.0
    }


def gcim_distributions(disagg_dst: pd.DataFrame,
                       imtl: float,
                       gmm_map: GmmMap, 
                       selection_imts: list[IMT],
                       conditioning_imt: IMT, 
                       corr_model_map: CorrMap,
                       ctx_builder: SiteRupContextBuilder,
                       occurence: bool,
                       percentiles=[0.16, 0.5, 0.84],
                       n_pts=200):
    # Disagg data must have columns "TRT", "Mag", "Dist"
    # Only a full disaggregation including TRT is implemented and the results
    disagg_type = "TRT_Mag_Dist_Eps"

    mu_cond, sig_cond, ordered_imts, site_rup_ctxs, sorted_weights, eps = \
        _conditional_im_distributions(
            disagg_type, disagg_dst, imtl, gmm_map, selection_imts, 
            conditioning_imt, corr_model_map, ctx_builder, occurence,
            n_samples=None, rng=None)

    df_stats, pdfs, cdfs = distribution_stats(
        ordered_imts, mu_cond, sig_cond, sorted_weights, percentiles, n_pts)

    return df_stats, pdfs, cdfs 


def gcim_simulation(n_samples: int,
                    disagg_dst: pd.DataFrame,
                    imtl: float,
                    gmm_map: GmmMap, 
                    selection_imts: list[IMT],
                    conditioning_imt: IMT, 
                    corr_model_map: CorrMap,
                    ctx_builder: SiteRupContextBuilder,
                    occurence: bool,
                    rup_rng: np.random.Generator,
                    sim_rng: np.random.Generator):
    # Disagg data must have columns "TRT", "Mag", "Dist"
    # Only a full disaggregation including TRT is implemented and the results
    disagg_type = "TRT_Mag_Dist_Eps"

    mu_cond, sig_cond, ordered_imts, site_rup_ctxs, sorted_weights, eps = \
        _conditional_im_distributions(
            disagg_type, disagg_dst, imtl, gmm_map, selection_imts, 
            conditioning_imt, corr_model_map, ctx_builder, occurence,
            n_samples=n_samples, rng=rup_rng)

    # loop through the columns (each column is one simulation)
    sims = generate_simulations(
        n_samples, mu_cond, sig_cond, ordered_imts, conditioning_imt, 
        site_rup_ctxs, corr_model_map, sim_rng)

    return sims, mu_cond, sig_cond, site_rup_ctxs, eps 


def generate_simulations(
        n_samples, mu_cond, sig_cond, imts, conditioning_imt, site_rup_ctxs, 
        corr_model_map, rng):
    # loop through the columns (each column is one simulation)
    sims = []
    for ii in range(mu_cond.shape[1]): # loop through the columns
        trt = site_rup_ctxs[ii].trt
        corr_df = corr_model_map[trt]
        corr_df = _sort_corr_df(
            corr_df, conditioning_imt.string, [imt.string for imt in imts])

        cond_rho_mat = conditional_correlation_matrix(
            corr_df, conditioning_imt.string)
        sims.append(generate_simulation(
            mu_cond[:, ii], sig_cond[:, ii], cond_rho_mat, rng))

    sims = pd.DataFrame(np.hstack(sims), 
                        index=[imt.string for imt in imts],
                        columns=list(range(n_samples)))
    return sims 


def generate_simulation(mu_cond, sig_cond, cond_rho_mat, 
                        rng:np.random.Generator):
    # mu_cond and sig_cond are vectors
    # rho_mat the conditional correlation matrices that correspond to 
    # mu_cond and sig_cond
    
    mu_cond = mu_cond.reshape((-1, 1))                      # shape (M,1)
    sig_cond = sig_cond.reshape((-1, 1))                    # shape (M,1)
    u = rng.standard_normal(len(mu_cond)).reshape(-1, 1)    # shape (M,1)
    L = np.linalg.cholesky(cond_rho_mat)                    # shape (M,M)
    v = L @ u

    # create a simulation of the the imt values for this M-R-TRT sample
    sim = np.exp(mu_cond + sig_cond * v)
    return sim


def _sort_corr_df(corr_df, conditioning_imt: str, imts: list[str]):
    im_order = imts + [conditioning_imt]
    corr_df = corr_df.loc[im_order, im_order]
    return corr_df


def _get_weights(disagg_dst, occurence):
    if occurence: # weights are from P(m|X=x) occurence probability
        weights = disagg_dst["P(m|X=x)"].to_numpy()
    else: # weights are from P(m|X>x) exceedance probability
        weights = disagg_dst["P(m|X>x)"].to_numpy()
    return weights


def _sort_ctxs_weights_by_trt(site_rup_ctxs, weights):
    trts = np.unique(site_rup_ctxs.trt)
    ctx_trt_subsets = {}
    weights_trt_subsets = {}
    for trt in trts:
        # Filter context for this specific TRT
        trt_mask = site_rup_ctxs.trt == trt
        ctx_trt_subsets[trt] = site_rup_ctxs[trt_mask]
        weights_trt_subsets[trt] = weights[trt_mask]

    sorted_weights = np.concatenate([w for w in weights_trt_subsets.values()])
    return ctx_trt_subsets, sorted_weights


def _group_imts(selection_imts):
    groups = dict.fromkeys([imt.name for imt in selection_imts])
    for k in groups.keys():
        groups[k] = [imt for imt in selection_imts if imt.name == k]
    return groups


def _calculate_eps_by_trt(ctx_trt_subsets, gmm_map, conditioning_imt, imtl):
    # calculate the conditioning epsilon for all ctxs
    epss = {}
    for trt, ctxs in ctx_trt_subsets.items():
        # calculate conditioning epsilon values
        gmm_im_cond = gmm_map[trt][conditioning_imt.name]
        m, s, *_ = _compute_gmm(gmm_im_cond, ctxs, [conditioning_imt])
        eps = (np.log(imtl) - m) / s
        epss[trt] = eps
    return epss


def conditional_im_distributions(
        site_rup_ctx: np.recarray,
        gmm: GMPE,
        imts: IMT|list[IMT],
        rhos: list[float]|np.ndarray,
        eps: list[float]|np.ndarray
        ) -> tuple[float, float, float, float]:
    
    # ctx  : recarray of length N
    # imts : list of length M
    # rho  : list of length N
    # eps  : list of length M

    rhos = np.array(rhos)[:, np.newaxis]  # shape (M, 1)
    eps = np.array(eps)[np.newaxis, :]    # shape (1, N)

    if isinstance(imts, IMT):
        imts = [imts]

    mu, sig, *_ = _compute_gmm(gmm, site_rup_ctx, imts)
    mu_cond = mu + sig * rhos * eps
    sig_cond = sig * np.sqrt(1 - rhos ** 2)
    
    return mu, sig, mu_cond, sig_cond


def _conditional_im_distributions(
        disagg_type: str,
        disagg_dst: pd.DataFrame,
        imtl: float,
        gmm_map: GmmMap, 
        selection_imts: list[IMT],
        conditioning_imt: IMT, 
        corr_model_map: CorrMap,
        ctx_builder: SiteRupContextBuilder,
        occurence: bool, 
        n_samples: float|None=None,
        rng: np.random.Generator|None=None):
    """ 
    Get the mu_cond and sig_cond for the selection_imts conditional on
    conditioning_imt and the imtl for n_samples from the disagg_dist.
    """
    # the disaggregation on epsilon doesn't affect the calculations so we can
    # collapse these rows to get the actual rupture scenarios.
    ddst = disagg_dst.groupby(["TRT", "Mag", "Dist"]).sum().reset_index()

    # filter out all the options that have a probability of zero as they
    # will contribute nothing to the gcim or samples and just add calculation time otherwise
    mask = ddst["P(m|X=x)"] != 0
    ddst = ddst.loc[mask, :]

    weights = _get_weights(ddst, occurence)

    # create the rupture contexts
    if n_samples is not None:
        rup_samples = ddst.sample(
        n=n_samples, replace=True, random_state=rng, weights=weights)
        # reassign weights so it has length n_samples
        weights = _get_weights(rup_samples, occurence)
    else:
        rup_samples = ddst
    site_rup_ctxs = get_site_rup_ctxs(disagg_type, rup_samples, ctx_builder)

    # sort ctxs and weights by trt
    ctx_trt_subsets, sorted_weights = _sort_ctxs_weights_by_trt(site_rup_ctxs, weights)

    # sort imts by imt.name
    imt_groups = _group_imts(selection_imts)

    # calculate the conditioning parameters first
    epss = _calculate_eps_by_trt(ctx_trt_subsets, gmm_map, conditioning_imt, imtl)
    
    # calculate the im distributions    
    mu_cond = []
    sig_cond = []
    for im_name, imts in imt_groups.items():
        mu_conds_for_im = [] 
        sig_conds_for_im = []

        for (trt, ctxs), eps in zip(ctx_trt_subsets.items(), epss.values()):
            # get rhos vector
            rhos = corr_model_map[trt].loc[conditioning_imt.string, 
                                           [imt.string for imt in imts]].to_numpy()
            
            # get appropriate gmm for im_name and trt
            gmm = gmm_map[trt][im_name]

            # calculate the conditional distributions
            _, _, mc, sc = conditional_im_distributions(
                ctxs, gmm, imts, rhos, eps.flatten())

            # append the results of the conditional distributions (hstack)
            mu_conds_for_im.append(mc)
            sig_conds_for_im.append(sc)

        mu_cond.append(np.hstack(mu_conds_for_im))
        sig_cond.append(np.hstack(sig_conds_for_im))
    mu_cond = np.vstack(mu_cond)
    sig_cond = np.vstack(sig_cond)
    
    ordered_imts = [imt for sublist in imt_groups.values() for imt in sublist] 
    ordered_ctxs = np.concatenate([ctx for ctx in ctx_trt_subsets.values()]).view(np.recarray)
    ordered_eps = [eps for sublist in epss.values() for eps in sublist]
    return mu_cond, sig_cond, ordered_imts, ordered_ctxs, sorted_weights, ordered_eps


def _compute_gmm(gmm, ctx, imts):
    mean = np.zeros((len(imts), len(ctx)))
    sig = np.zeros_like(mean)
    tau = np.zeros_like(mean)
    phi = np.zeros_like(mean)
    gmm.compute(ctx, imts, mean, sig, tau, phi)
    return mean, sig, tau, phi


def _distribution_mean(mu_cond, sig_cond, weights, axis):
    # expectation of a mixed distribution of log normal random variables
    mean_val = np.sum(weights * np.exp(mu_cond + 0.5*sig_cond**2), axis=axis)
    mean_val = mean_val.reshape((mu_cond.shape[0], -1))
    return mean_val


def _distribution_percentiles(percentiles, mu_cond, sig_cond, 
                             weights, im_low, im_high):
    # mu_cond: (M, N) matrix
    # im_low/im_high: (M,) vectors for each IMT
    M = mu_cond.shape[0]
    num_p = len(percentiles)
    results = np.zeros((M, num_p))

    for i in range(M):
        # Define the mixture CDF for the i-th IMT only
        # target_x is a scalar passed by the solver
        def mixture_cdf(target_x):
            # cdfs is shape (N,)
            cdfs = stats.lognorm.cdf(target_x, s=sig_cond[i, :], 
                                     scale=np.exp(mu_cond[i, :]))
            # Return scalar weighted sum (Equation 8)
            return np.sum(weights * cdfs)

        # Solve for each requested percentile
        for j, p in enumerate(percentiles):
            results[i, j] = brentq(lambda x: mixture_cdf(x) - p, 
                                   im_low[i], im_high[i])
    return results
    

def distribution_stats(imts, mu_cond, sig_cond, weights, percentiles, n_pts):
    # calculate the stats of the im distributions
    pdfs = {}
    cdfs = {}
    mean_vals = []
    percentile_vals = []
    
        # determine the pdf bounds and evaluate pdfs and cdfs
    im_low = np.exp(np.min(mu_cond - 3 * sig_cond, axis=1))  #(M,)
    im_high = np.exp(np.max(mu_cond + 3 * sig_cond, axis=1))
    # Create an (M, n_pts) matrix of x-values
    x_vals = np.array([np.geomspace(il, ih, n_pts) for il, ih in zip(im_low, im_high)])

    # mu/sig: (M, N) -> (M, N, 1)
    # x_vals: (M, n_pts,) -> (M, 1, n_pts)
    mu_3d = mu_cond[:, :, np.newaxis]
    sig_3d = sig_cond[:, :, np.newaxis]
    x_3d = x_vals[:, np.newaxis, :]

    # 3. Calculate 3D PDF Matrix
    # Resulting shape: (M, N, 500)
    pdf_matrix = stats.lognorm.pdf(x_3d, s=sig_3d, scale=np.exp(mu_3d))
    cdf_matrix = stats.lognorm.cdf(x_3d, s=sig_3d, scale=np.exp(mu_3d))

    # 4. Sum together weighted by rupture probabilities (Eq. 8)
    # weights (N,) broadcast to (1, N, 1)
    # Final shape: (M, n_pts) - One PDF/CDF with n_pts per IMT
    total_pdfs = np.sum(pdf_matrix * weights[np.newaxis, :, np.newaxis], axis=1)
    total_cdfs = np.sum(cdf_matrix * weights[np.newaxis, :, np.newaxis], axis=1)
    for imt, pdf, cdf, xs in zip(imts, total_pdfs, total_cdfs, x_vals):
        pdfs[imt.string] = np.hstack([xs.reshape(-1, 1), 
                                      pdf.reshape(-1, 1)])
        cdfs[imt.string] = np.hstack([xs.reshape(-1, 1), 
                                      cdf.reshape(-1, 1)])
            
    # calculate point stats
    mean_vals = _distribution_mean(mu_cond, sig_cond, weights, axis=1)
    percentile_vals = _distribution_percentiles(
        percentiles, mu_cond, sig_cond, weights, im_low, im_high)

    # organise the outputs:
    # vstack all of the results
    all_stats = np.hstack([mean_vals, percentile_vals])

    stats_labels = ["mean"] + [f"p{int(p*100)}" for p in percentiles]

    df_stats = pd.DataFrame(all_stats, 
                            index=[im.string for im in imts],
                            columns=stats_labels)

    return df_stats, pdfs, cdfs 


def ensemble_ks_test(record_ensemble_im_df: pd.DataFrame, target_cdfs: dict,
                     p_val=0.05):

    # get the ecdf of selected records for each im:   
    im_df = record_ensemble_im_df
    kst_results = {}

    for im, values in im_df.T.iterrows():
        # 2 sample ks-test with null hypothesis that the distributions are drawn
        # from teh same parent distribution

        if im not in target_cdfs.keys():
            # key error woulde raised for the conditioning im because it has no
            # target cdf. Therefore it is skipped.
            continue
        
        def cdf(x):
            xp = target_cdfs[im][:,0]
            fp = target_cdfs[im][:,1]
            return np.interp(x, xp, fp, left=0, right=1)

        kst_results[im] = stats.ks_1samp(values.to_numpy(), cdf)

    # check if all ims pass
    failed_ims = []
    for im, kst, in kst_results.items():
        # if the statistic is less than p_val we must reject the null hypothesis
        # and say that the distributions are statisitcally different
        if kst.pvalue < p_val or np.isnan(kst.pvalue):
            failed_ims.append(im) 

    passed = False if len(failed_ims) > 0 else True
    return passed, kst_results, failed_ims


def ensemble_ecdfs(record_ensemble_ims_df: pd.DataFrame):
    im_df = record_ensemble_ims_df
    ecdfs = {}
    for im, values in im_df.T.iterrows():
        ecdf = stats.ecdf(values.to_numpy()).cdf

        x_vals = np.repeat(ecdf.quantiles, 2)
        y_vals = np.repeat(ecdf.probabilities, 2)
        y_vals = np.concatenate(([0], y_vals[:-1]))

        ecdfs[im] = np.column_stack([x_vals, y_vals])
    return ecdfs


def ensemble_ks_bounds(target_cdfs: dict, n_sims: int, p_val:float=0.05) -> dict:
    # get upper and lower bounds for the specified p-val
    n_a = n_sims # number of samples in the sample distribution

    ks_bounds = {}
    for k, cdf in target_cdfs.items():
        x_vals = np.append([0], cdf[:, 0]).reshape((-1, 1))
        y_vals = np.append([0], cdf[:, 1]).reshape((-1, 1))
        n_b = len(cdf)  #  number of samples in the parent distribution (CDF) 
        n_eff = np.round((n_a * n_b) / (n_a + n_b))

        # kstwo.ppf is the inverse of the cdf of the two-sided Smirnov distribution.
        # By using we 1-p_val we are finding the value of D that is exceeded only
        # p-val of the time by random chance. i.e. "I think a p_val fluke is rare and
        # bet that the hypothesis that they are the same is wrong" 
        D_pval = stats.kstwo.ppf(1-p_val, n_eff)
        upper = np.clip(y_vals + D_pval, 0, 1).reshape((-1, 1))
        lower = np.clip(y_vals - D_pval, 0, 1).reshape((-1, 1))

        ks_bounds[k] = np.hstack([x_vals, lower, y_vals, upper])

    # dictionary: values = np.ndarray with columns for x_vals, lower_ks_bournd, 
    # cdf, upper_ks_bound
    return ks_bounds

    
def ensemble_r_score(ks_results: dict, weights: np.ndarray):
    # score of the misfit of the total ensemble, following Bradley (2012) 
    score = np.array([wi * ks.statistic ** 2 
                      for ks, wi in zip(ks_results.values(), weights)]).sum()
    return score
