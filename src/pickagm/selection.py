"""
    Module with algorithms for record selection

"""
import numpy as np
import pandas as pd
from numpy.typing import NDArray
from typing import Callable

from openquake.hazardlib.imt import IMT

from pickagm.distributions import SCALING_EXPONENTS

def sse_multi(
        eqdf: pd.DataFrame, im_sims: np.ndarray | pd.DataFrame) -> pd.DataFrame:
    """
    Compute (sum(ln(eqdb) - ln(im_sim))**2).sum(axis=1) for multiple im_sim 
    realizations.
    
    Parameters
    ----------
    eqdf : pd.DataFrame
        Data matrix (n_record x n_ims)
    im_sims : np.ndarray or pd.DataFrame
        Realizations matrix (n_realizations x n_ims)
    
    Returns
    -------
    pd.DataFrame
        Result of shape (n_rows, n_realizations)
    """
    df_array = np.log(eqdf.to_numpy())                 # (n_records, n_ims)
    vecs = np.log(np.asarray(im_sims))                 # (n_realizations, n_ims)
    
    # Compute all differences: 
    # broadcast df_log[:, None, :] against vecs_log[None, :, :]
    # resulting shape (n_records, n_realizations, n_ims)
    diff = df_array[:, None, :] - vecs[None, :, :]   
    
    # Square and sum over columns
    result = np.sum(diff ** 2, axis=2) # (n_records, n_realizations)                
    
    # Wrap back into DataFrame for nice labeling
    return pd.DataFrame(result, index=eqdf.index, columns=im_sims.index)


def scale_records_in_db(
        df: pd.DataFrame, conditioning_imt: IMT, conditioning_value: float
        ) -> tuple[pd.DataFrame, pd.Series]:
    """ 
    Assumes the conditioning_imt is in some way amplitude dependent. Trying to
    use RSD595 as a conditioning_imt will raise an error because it is not
    amplitude dependent

    Assumes the dataframe has a one level index which contains the im tags
    Assumes the dataframe contains only a single component. No imts are repeated
    
    """
    # (n_records x 1 im)
    sf_cond_imt = conditioning_value / df[str(conditioning_imt)] 
    gamma = get_sf_exponent(str(conditioning_imt))
    # (n_records x 1 im) amplitude sf
    alpha = sf_cond_imt.to_numpy() ** (1 / gamma)  
    
    sdf = scale_record_imts(df, alpha)
    
    return sdf, alpha


def scale_record_imts(df: pd.DataFrame, alpha: NDArray) -> pd.Series:
    # get the scale factors for the other ims
    
    sdf = df.copy()
    remaining_ims = list(sdf.columns)
    for imt_str in SCALING_EXPONENTS.keys():
        cols_to_scale = [c for c in remaining_ims if c.startswith(imt_str)] 
        if cols_to_scale:
            # remove ims from remaining ims
            remaining_ims = [c for c in remaining_ims if not c in cols_to_scale]
            
            # calculate the im scale factor    
            sf = alpha ** get_sf_exponent(imt_str) # (n_records x 1 im)

            # apply the scale factor to the columns   
            sdf.loc[:, cols_to_scale] = sdf.loc[:, cols_to_scale].mul(sf, axis=0)

    # if there are any tags in remaining_ims
    if remaining_ims:
        raise ValueError("Not all IMTs were scaled")
    return sdf


def get_sf_exponent(imt_str: str) -> float:
    if imt_str.startswith("SA"):
        return SCALING_EXPONENTS["SA"]
    else:
        return SCALING_EXPONENTS[imt_str]
    

def rank_by_score(
        score_df: pd.DataFrame, k: int | None = None) -> list[list[tuple]]:
    """

    """
    # number of records to return
    nvals = score_df.shape[0]
    m = score_df.shape[1]
    if k is None:
        k = min(m, nvals)

    # get integer positions of k smallest per column
    arr = score_df.to_numpy()
    idx = np.argsort(arr, axis=0)[:k, :]  # shape (k, m)

    # get the indices of the score_df corresponding to the ranks
    ranked_indices = []
    for ranked_records in idx.T:
        ranked_indices.append(score_df.index[ranked_records].to_list())

    return ranked_indices


def stack_components(
        eqdf: pd.DataFrame, component_tags: list[str]) -> list[pd.DataFrame]:
    
    dfs = [eqdf[tag].copy() for tag in component_tags]
    eqdf = pd.concat(dfs, axis=0, keys=component_tags)
    eqdf = eqdf.swaplevel(0, 1)
    eqdf = eqdf.sort_index(level=0)

    return eqdf
    

def scale_multiple_dbs(
        eqdfs: list[pd.DataFrame], 
        conditioning_imt: IMT, 
        conditioning_im_level: float
        ) -> dict[str, list[pd.DataFrame|pd.Series]]:

    scaled_eqdbs = [
        scale_records_in_db(db, conditioning_imt, conditioning_im_level) 
        for db in eqdfs]
    dfs = [i[0] for i in scaled_eqdbs]
    alphas = [i[1] for i in scaled_eqdbs]
    result = {"eqdfs" : dfs,
              "alphas": alphas}
    return  result  


# def get_ims_of_selected_records(
#         eqdf: pd.DataFrame, 
#         selected_indices: list,
#         selected_alphas: pd.Series,
#         scaled: bool=True
#         ) -> pd.DataFrame:
    
#     gm_ims = eqdf.loc[selected_indices]
#     if scaled:
#         gm_ims = scale_record_imts(gm_ims, selected_alphas.to_numpy())
    
#     return gm_ims


def best_individual_components(
        eqdf: pd.DataFrame,
        component_tags: list[str],
        conditioning_imt: IMT,
        conditioning_imt_level: float,
        realisations: pd.DataFrame,
        scoring_func: Callable[[pd.DataFrame, pd.DataFrame], pd.DataFrame],
        scaled: bool=True,
        max_sf: float|None=None,
        min_sf: float|None=None
        ) -> pd.DataFrame:
    """
    select the best records using best_choices and the finds the corresponding
    scaling factor using component tags and alphas and return the appropriately
    scaled ims

    The input DataFrame, df, must have a two-level MultiIndex on columns:
      level 0: integer rank
      level 1: ['dir', 'idx']

    'idx' is the index of the record in thr original databank
    'dir' is the component of the record (e.g. "U" or "V" for the ESM-database)

    component_tags is the list containing the possible values of 'dir'
    alphas is a list of pd.Series of n_records len where each value is the 
    amplitude of that record

    Assumes both the realisations and the ims in teh eq database are not in
    the form ln(IM)
    """
    eqdf_stacked = stack_components(eqdf, component_tags) 
    eqdf_scaled, alphas = scale_records_in_db(
    eqdf_stacked, conditioning_imt, conditioning_imt_level)
    
    eqdf_scaled, alphas = filter_based_on_scale_factor(
        eqdf_scaled, alphas, max_sf, min_sf)
    
    scores = scoring_func(eqdf_scaled, realisations)
    ranked_indices = rank_by_score(scores)
    selected_indices = select_indices(ranked_indices)

    # get ims of selected records:
    if scaled:
        ims_of_selected_records = eqdf_scaled.loc[selected_indices]
    else:
        ims_of_selected_records = eqdf_stacked.loc[selected_indices]

    return ims_of_selected_records


def filter_based_on_scale_factor(
        df: pd.DataFrame, alphas: pd.Series, 
        max_sf: float|None, min_sf: float|None
        ) -> tuple[pd.DataFrame, pd.Series]:
    
    alphas = pd.Series(alphas, index=df.index)

    if max_sf and min_sf:
        mask = (alphas < max_sf) & (alphas > min_sf)
    elif max_sf and (not min_sf):
        mask = alphas < max_sf
    elif (not max_sf) and min_sf:
        mask = alphas > min_sf
    else:
        return df, alphas
    return df[mask], alphas[mask]


def select_indices(index_lists: list[list]) -> list[tuple]:
    """ 
    select the highest ranking component for each column only if it has not
    already been selected for another record.

    one list for each simulated set of IMs
    each list contains the indices of the best fit records/components in order
    """
    selected_indices = []

    for index_list in index_lists:
        for index in index_list:
            if index not in selected_indices:
                selected_indices.append(index)
                break
    return selected_indices
    