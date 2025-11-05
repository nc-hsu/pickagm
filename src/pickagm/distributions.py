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
from scipy.stats._distn_infrastructure import rv_frozen

NormalRV = rv_frozen  # alias for readability

class IMDistribution:
    """
    Target distribution for a ground motion intensity measure.
    """
    
    def __init__(self, imt: IMT, gmm: GMPE, ctx: np.recarray, 
                 conditioning_imt: IMT|None,
                 conditioning_eps: float, rho: float) -> None:
        """_summary_

        Args:
            imt (IMT): _description_
            gmm (GMPE): _description_
            ctx (np.recarray): _description_
            corr_model (CorrelationModel): _description_

        Returns:
            _type_: _description_
        """
        self.imt = imt
        self.gmm = gmm
        self.ctx = ctx
        self.conditioning_imt = conditioning_imt
        self.conditioning_eps = conditioning_eps
        self.rho = rho
        self.uncond = self._unconditional_distribution()

        # conditioning values are provided
        if (conditioning_eps is not None) and (rho is not None):
            self.cond = self._conditional_distribution(conditioning_eps, rho)
        
        # insufficient condtioning values are provided
        elif (conditioning_eps is None) is not (rho is None):
            raise ValueError("conditioning_eps and rho must both "
                             "be None, or both have values")
        else:
            self.cond = None

        return None

    def _conditional_distribution(self, conditioning_eps:float, rho_ij: float
            ) -> NormalRV:
        
        mean, sigma = self.compute_gmm()
        conditional_mean = mean + rho_ij * conditioning_eps * sigma
        conditional_sigma = sigma * np.sqrt(1 - rho_ij**2)
        rv = stats.norm(loc=conditional_mean, scale=conditional_sigma)
        return rv

    def _unconditional_distribution(self) -> NormalRV:
        # just return the gmm result 
        mean, sigma = self.compute_gmm()
        rv = stats.norm(loc=mean, scale=sigma)
        return rv

    def compute_gmm(self) -> tuple[float, float]:
        # Compute mean and standard deviation using gmm and ctx
        N = len(self.ctx[self.ctx.dtype.names[0]]) # sites or rups in ctx
        M = 1  # only this class only handles one IMT at a time
        
        mean = np.zeros((M, N))
        sigma = np.zeros((M, N))
        tau = np.zeros((M, N))
        phi = np.zeros((M, N))

        self.gmm.compute(self.ctx, [self.imt], mean, sigma, tau, phi)
        return mean.flatten()[0], sigma.flatten()[0]
    
    def cdf(self, eps: float, cond: bool, npts:int=50) -> NDArray:
        # returns either the conditional or unconditional cdf between X values
        # corresponding to mu - eps*sigma and mu + eps* sigma

        if cond: rv = self.cond 
        else: rv = self.uncond

        x1 = rv.mean() - eps * rv.std()
        x2 = rv.mean() + eps * rv.std()
        xs = np.linspace(x1, x2, npts)

        return np.vstack([xs, rv.cdf(xs)]).T


    def __repr__(self) -> str:
        return f"TargetDistribution(imt={self.imt}, gmm={self.gmm}, " \
               f"gmm={self.ctx}, conditioning_imt={self.conditioning_imt}, "\
               f"conditioning_eps={self.conditioning_eps}, rho={self.rho})"


class CorrelationModel:
    """
    correlation model classes - used to construct a covariance matrix
    for an arbitrary number of ims as pandas data frame.
    
    Acts as a wrapper for the pandas data frame correlation matrix.
    """
    def __init__(self, model: str|None=None, tag:str="", 
                 only_imts:list[IMT]=[],
                 exclude_imts:list[IMT]=[], 
                 exclude_SA_T_gt: IMT|None=None,
                 exclude_SA_T_lt: IMT|None=None) -> None:
        # load the correlation model from data folder
        if model is not None:
            self.corr = load_correlation_model(
                model, only_imts, exclude_imts, 
                exclude_SA_T_gt, exclude_SA_T_lt)
        else:
            self.model = model
            self.corr = pd.DataFrame()
        
        if (model is not None) and tag == "":
            self.tag = model
        else:   
            self.tag = tag

    @property
    def correlations(self):
        return self.corr

    def get_correlation(self, im1: IMT, conditioning_im: IMT) -> float:
        """
        Get the correlation value for two imts.
        """
        # check if the IM exists in the correlation model
        if str(im1) not in self.corr.columns:
            if not str(im1).startswith("SA"):
                # raise error
                raise ValueError(f"IM {im1} not found in correlation model" 
                                 f"{self.model}")
            else:
                # try to interpolate column for the missing SA(im1)
                T1 = im1.period
                # Series indexed by rows (conditioning IMs)
                im_series = interpolate_correlations_for_sa(self.corr, T1)  

                # If the conditioning IM exists in the index, return directly
                if str(conditioning_im) in im_series.index:
                    return float(im_series.at[str(conditioning_im)])

                # If the conditioning IM is an SA and is not present in the 
                # rows, interpolate along the rows
                if str(conditioning_im).startswith("SA"):
                    T2 = conditioning_im.period
                    sa_pattern = re.compile(r"SA\(([\d.]+)s\)")

                    # Extract available SA rows (period, value) from the 
                    # interpolated column
                    sa_items = {}
                    for idx, val in im_series.items():
                        m = sa_pattern.match(idx)
                        if m:
                            sa_items[float(m.group(1))] = float(val)

                    if not sa_items:
                        raise ValueError(f"No SA rows available to interpolate " 
                                         f"conditioning IM {conditioning_im}")

                    # Sort by period and interpolate to T2
                    sa_items.sort(key=lambda x: x[0])
                    periods = np.array([p for p, _ in sa_items])
                    values = np.array([v for _, v in sa_items])

                    interp_val = np.interp(
                        T2, periods, values, left=np.nan, right=np.nan)
                    if np.isnan(interp_val):
                        raise ValueError(f"Could not interpolate correlation " 
                                         f"for {conditioning_im} vs {im1}")
                    return float(interp_val)

                # conditioning_im not there and not in SA -> cannot interpolate
                raise ValueError(f"Conditioning IM {conditioning_im} not "
                                 f"found in correlation model {self.model}")
            
        corr_value = self.corr.at[str(conditioning_im), str(im1)]
        return corr_value

    def get_SA_periods(self) -> list[float]:
        """
        Get the available spectral acceleration periods in the correlation 
        model.
        """
        sa_periods = []
        for imt_str in self.corr.columns:
            period = period_from_SA_IMT_str(imt_str)
            if period:
                sa_periods.append(period)
        return sa_periods
    
    def add_correlations_to_model(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Add rows and columns to the correlation corr dataframe.
        The DataFrame should have columns corresponding to the IMs in the 
        correlation model.
        """
        ... # todo:: be implemented

    def __repr__(self):
        return f"CorrelationModel(model={self.model})"

    def __str__(self):
        return f"{self.model}"


def period_from_SA_IMT_str(imt_str: str):
    sa_pattern = re.compile(r"SA\(([\d.]+)\)")
    m = sa_pattern.match(imt_str)
    if m:
        return float(m.group(1))
    

def load_correlation_model(correlation_model: str, 
                           only_imts: list[IMT]=[],
                           exclude_imts: list[IMT]=[], 
                           exclude_SA_T_gt: IMT|None=None,
                           exclude_SA_T_lt: IMT|None=None) -> pd.DataFrame:
    with res.files("pickagm.data").joinpath(f"{correlation_model}.csv") \
        .open("r", encoding="utf-8") as f:
        df = pd.read_csv(f, sep=",", index_col=[0], header=[0])

    if only_imts != []:
        imt_names = [str(imt) for imt in only_imts]
        df = df.loc[imt_names, imt_names] 

    if exclude_imts != []:
        imt_names = [str(imt) for imt in exclude_imts]
        df = df.drop(index=imt_names, columns=imt_names)

    if exclude_SA_T_gt is not None:
        imt_names = [imt for imt in df.columns 
                     if imt.startswith("SA") and 
                     (period_from_SA_IMT_str(imt) > 
                      period_from_SA_IMT_str(str(exclude_SA_T_gt))) 
                     ]
        df = df.drop(index=imt_names, columns=imt_names)

    if exclude_SA_T_lt is not None:
        imt_names = [imt for imt in df.columns 
                     if imt.startswith("SA") and 
                     (period_from_SA_IMT_str(imt) < 
                      period_from_SA_IMT_str(str(exclude_SA_T_lt))) 
                     ]
        df = df.drop(index=imt_names, columns=imt_names)

    return df


def interpolate_correlations_for_sa(
        df: pd.DataFrame | pd.Series, T, im_prefix="SA"):
    """
    Interpolate correlation values for spectral accelerations.

    - If `df` is a DataFrame: returns a pandas.Series (for scalar T) or 
      DataFrame (for multiple T)
      containing the interpolated column(s) indexed by the original rows.
    - If `df` is a Series: returns a single float (T must be scalar).

    Parameters
    ----------
    df : pandas.DataFrame or pandas.Series
        Data containing SA columns like "SA(0.02s)", "SA(0.05s)", ...
        If DataFrame, interpolation is performed for each row independently.
        If Series, interpolation is performed on that single series.
    T : float or array-like
        Target period(s) [s] for which to interpolate SA.
    im_prefix : str, default="SA"
        The prefix identifying the spectral acceleration columns.
    """
    pattern = re.compile(rf"{im_prefix}\(([\d.]+)\)")

    # Helper to extract mapping period -> column/key
    def _extract_keys_from_labels(labels):
        keys = {float(m.group(1)): lbl for lbl in labels 
                if (m := pattern.match(lbl))}
        if not keys:
            raise ValueError(f"No columns/labels found for IM type "
                             f"'{im_prefix}'")
        periods_sorted = np.array(sorted(keys.keys()))
        cols_sorted = [keys[p] for p in periods_sorted]
        return periods_sorted, cols_sorted

    T_arr = np.atleast_1d(T)

    # Case: Series input -> return a single float (require scalar T)
    if isinstance(df, pd.Series):
        if T_arr.size != 1:
            raise ValueError("When input is a Series, T must be a scalar value")
        periods, keys = _extract_keys_from_labels(df.index)
        values = np.array([float(df[k]) for k in keys])
        interp_val = np.interp(
            float(T_arr[0]), periods, values, left=np.nan, right=np.nan)
        if np.isnan(interp_val):
            raise ValueError(f"Could not interpolate SA at T={T_arr[0]}")
        return float(interp_val)

    # Case: DataFrame input -> return Series (for scalar T) or DataFrame 
    # (for multiple T)
    if isinstance(df, pd.DataFrame):
        periods, cols = _extract_keys_from_labels(df.columns)
        corr_data = df[cols].to_numpy()  # shape: (n_rows, n_periods)

        results = {}
        for t in T_arr:
            # interpolate each row along the SA-period axis
            interp_vals = np.apply_along_axis(
                lambda row: np.interp(
                    float(t), periods, row, left=np.nan, right=np.nan),
                axis=1,
                arr=corr_data
            )
            col_name = f"{im_prefix}({float(t):.3f})"
            results[col_name] = interp_vals

        res_df = pd.DataFrame(results, index=df.index).map(float)
        return res_df.iloc[:, 0] if res_df.shape[1] == 1 else res_df

    raise TypeError("df must be a pandas.DataFrame or pandas.Series")


def correlation_matrix(
        list_of_imts: list[IMT], correlation_model: CorrelationModel,
        conditioning_imt: IMT|None=None, incl_conditioning_imt: bool=False
        ) -> pd.DataFrame|tuple[pd.DataFrame, pd.Series]:
    
    """returns the correlation matrix of the chosen IMs depending on the
    correlation model. Periods of SA(T) are interpolated if they don't exist in
    the model already.

    Error will be raised if a non-SA imt is in list_of_imts but not in
    correlation_model

    If conditioning_imt != None then an additional vector is returned which is 
    a  (N,1) array of the correlations between each of the N dependent imts and 
    the conditioning imt. 
    
    If ncl_conditioning_imt is False and conditioning_imt is already in 
    list_of_imts the the correlation vector is removed from the correlation 
    matrix.

    if 

    Returns:
        _type_: _description_
    """
    
    if conditioning_imt and (not conditioning_imt in list_of_imts):
        imts_to_check = list_of_imts + [conditioning_imt]
    else:
        imts_to_check = list_of_imts

    SA_imts_not_in_model = [imt for imt in imts_to_check
                            if str(imt) not in \
                            correlation_model.correlations.index and 
                            str(imt).startswith("SA")]

    imt_correlations = correlation_model.correlations.copy()

    # interpolate new SA ims if needed and add to dataframe
    SA_imts_not_in_model = [imt for imt in SA_imts_not_in_model 
                            if str(imt).startswith("SA")]
    for imt in SA_imts_not_in_model:
        interp_values = interpolate_correlations_for_sa(
            correlation_model.correlations, period_from_SA_IMT_str(str(imt))
        )
        
        # add the column/row for this imt to the end. correlation is 1 with
        # self is 1
        imt_correlations[str(imt)] = interp_values
        imt_correlations.loc[str(imt), :] = interp_values.tolist() + [1.0]

    keys = [str(imt) for imt in imts_to_check]
    corr_matrix = imt_correlations.loc[keys, keys]

    if conditioning_imt:
        if str(conditioning_imt) in imt_correlations.columns:
            corr_vector = corr_matrix[str(conditioning_imt)] 

            if not incl_conditioning_imt:
                corr_matrix = corr_matrix.drop(index=str(conditioning_imt))
                corr_matrix = corr_matrix.drop(columns=str(conditioning_imt))
                corr_vector = corr_vector.drop(labels=str(conditioning_imt))

            return corr_matrix, corr_vector
        else:
            raise KeyError(f"imt {conditioning_imt} not in correlation matrix")
    
    return corr_matrix


def unconditional_covariance_matrix(
        imt_distributions: list[IMDistribution],
        correlation_matrix: pd.DataFrame) -> pd.DataFrame:
    
    sig_vec = np.array(
        [imtd.uncond.std() for imtd in imt_distributions]
        ).reshape(-1,1)
    sig_mat = np.matmul(sig_vec, sig_vec.T)

    labels = correlation_matrix.index
    cm = correlation_matrix.to_numpy()
    return pd.DataFrame( cm * sig_mat, index=labels, columns=labels) 


def covariance_vector(
        imt_distributions: list[IMDistribution],
        conditioning_imt_distribution: IMDistribution,
        conditioning_imt_correlation_vector: pd.Series, 
        ) -> pd.Series:
    
    sigmas = np.array([imtd.uncond.std() for imtd in imt_distributions]).flatten()
    Sigma_One = conditioning_imt_correlation_vector.to_numpy() * \
        sigmas * conditioning_imt_distribution.uncond.std()
    labels = conditioning_imt_correlation_vector.index
    return pd.Series(Sigma_One, index=labels)


def conditional_mean_vector(
        imt_distributions: list[IMDistribution]) -> pd.Series:
    labels = [str(imtd.imt) for imtd in imt_distributions]
    mean_values = []
    for imtd in imt_distributions:
        mv = imtd.cond.mean()
        if not np.isfinite(mv): 
            # catches the case where the contioninal IM has sigma == 0 and 
            # the rv can't calculate the mean and instead returns nan
            mv = imtd.cond.kwds["loc"]
        mean_values.append(mv)
    
    return pd.Series(mean_values, index=labels)


def conditional_sigma_vector(
        imt_distributions: list[IMDistribution]) -> pd.Series:
    labels = [str(imtd.imt) for imtd in imt_distributions]
    sigma_values = []
    for imtd in imt_distributions:
        sv = imtd.cond.std()
        if not np.isfinite(sv): 
            # catches the case where the contioninal IM has sigma == 0 and 
            # the rv can't calculate the std and instead returns nan
            sv = imtd.cond.kwds["scale"]
        sigma_values.append(sv)
    return pd.Series(sigma_values, index=labels)


def conditional_covariance_matrix(
        Sigma_0: pd.DataFrame, Sigma_1: pd.Series,
        conditioning_imt_distribution: IMDistribution) -> pd.DataFrame:
    
    labels = Sigma_0.index
    sigma_star = conditioning_imt_distribution.uncond.std()
    Sigma_1 = Sigma_1.to_numpy().reshape(-1,1)
    Sigma = Sigma_0.to_numpy() - (1 /  sigma_star ** 2) * \
          np.matmul(Sigma_1, Sigma_1.T)
    res = pd.DataFrame(Sigma, index=labels, columns=labels)
    return res


def create_im_distributions(
        imts: list[IMT], correlation_model: CorrelationModel,
        gmm_mapping: dict[str, GMPE], ctx: np.recarray, 
        conditioning_imt: IMT, conditioning_eps: float
        ) -> tuple[IMDistribution, list[IMDistribution]]:
    
    """
    create the im distributions used for record selection.

    imts: dependent IMTs

    Assumes:
    -   all imts are in the correlation model or can be interpolated from values
        that are in there

    """
       
    im_dists = []
    for imt in imts:
        if str(imt).startswith("SA"):
            key = "SA"
        else:
            key = str(imt)

        rho = correlation_model.get_correlation(imt, conditioning_imt)

        im_dists.append(IMDistribution(
            imt, gmm_mapping[key], ctx, conditioning_imt, conditioning_eps, rho 
        ))

    # distribution for the conditioning imt
    if str(conditioning_imt).startswith("SA"):
        key = "SA"
    else:
        key = str(conditioning_imt)

    imd_star = (IMDistribution(
        conditioning_imt, gmm_mapping[key], ctx, 
        conditioning_imt, conditioning_eps, 1.0 
        ))

    return imd_star, im_dists


def get_cs(imt_distributions: list[IMDistribution], 
           conditioning_imt_distribution: IMDistribution) -> NDArray:
    
    cond_imt_dist = conditioning_imt_distribution
    # filter out the imts that are SA
    if str(conditioning_imt_distribution.imt).startswith("SA"):
        sa_imtds = [imtd for imtd in imt_distributions + [cond_imt_dist]
                     if str(imtd.imt).startswith("SA")]
    else:
        sa_imtds = [imtd for imtd in imt_distributions
                     if str(imtd.imt).startswith("SA")]
        
    periods = [imtd.imt.period for imtd in sa_imtds]
    mus = [imtd.cond.mean() for imtd in sa_imtds]
    sigmas = [imtd.cond.std() for imtd in sa_imtds]

    points = [(p,m,s) for p, m, s in zip(periods, mus, sigmas)]
    sorted_points = sorted(points, key=lambda x: x[0])

    return np.array(sorted_points)


def sample_im_distribution(
        mean_vec: pd.Series, covariance_matrix: pd.DataFrame, n: int=1
        ) -> pd.DataFrame:
    """ use 'svd' method so that the semi-definite problem can be solved

    Args:
        mean_vec (pd.Series): _description_
        covariance_matrix (pd.DataFrame): _description_
        n (int, optional): _description_. Defaults to 1.

    Returns:
        pd.DataFrame: _description_
    """
    labels = mean_vec.index
    rng = np.random.default_rng()
    realisation = rng.multivariate_normal(
        mean_vec.to_numpy(), covariance_matrix.to_numpy(), 
        size=n, method="svd")

    return pd.DataFrame(realisation.T, index=labels)


def get_empirical_cdfs(df: pd.DataFrame) -> pd.Series:
    """_summary_

    Args:
        df (pd.DataFrame): a Dataframe where each row corresponds to an imt and
        each column a different simulation / record

    Returns:
        pd.Series: a series where each cell is scipy.stats.ecdf object
    """
    # Apply ECDF row-wise
    return df.apply(lambda row: stats.ecdf(row.to_numpy()), axis=1)

# todo:: create a correlation matrix for the BJ08 model
