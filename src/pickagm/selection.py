"""
    Module with algorithms for record selection

"""
import numpy as np
import pandas as pd
from numpy.typing import NDArray
from typing import Callable, Protocol

from openquake.hazardlib.imt import IMT
from openquake.hazardlib.gsim.base import GMPE

from pickagm.distributions import (
    SCALING_EXPONENTS, ensemble_ks_test, ensemble_r_score)
from pickagm.corrmodels import conditional_correlation_matrix
from pickagm.typing import GmmMap, CorrMap, IMTName
from pickagm.contexts import SiteRupContextBuilder, get_rup_ctxs


def _filter_on_bounds(df, col, bounds):
    lb, ub = bounds
    if lb and ub:
        df = df[(df[col] >= lb) & (df[col] <= ub)]
    elif lb:
        df = df[df[col] >= lb]
    elif ub:
        df = df[df[col] <= ub]
    return df


def preselection(df, magnitude_bounds, distance_bounds, vs30_bounds, usable_T):
    # preselection filtering

    if magnitude_bounds is not None:
        df = _filter_on_bounds(df, ("metadata", "mag"), magnitude_bounds)

    if distance_bounds is not None:
        df = _filter_on_bounds(df, ("metadata", "rjb"), distance_bounds)

    if vs30_bounds is not None:
        df = _filter_on_bounds(df, ("metadata", "vs30"), vs30_bounds)

    if usable_T is not None:
        df = df[df[("metadata", "max_usable_T")] >= usable_T]

    return df


def filter_on_sf(df, sf_bounds: tuple[float, float]|None):
    if sf_bounds is not None:
        df = _filter_on_bounds(df, ('metadata', 'alpha'), sf_bounds)
    return df


def filter_on_m(df, magnitude_bounds: tuple[float, float]|None):
    if magnitude_bounds is not None:
        df = _filter_on_bounds(df, ("metadata", "mag"), magnitude_bounds)
    return df


def filter_on_d(df, distance_bounds: tuple[float, float]|None):
    if distance_bounds is not None:
        df = _filter_on_bounds(df, ("metadata", "rjb"), distance_bounds)
    return df
                    

def filter_on_vs30(df, vs30_bounds: tuple[float, float]|None):
    if vs30_bounds is not None:
        df = _filter_on_bounds(df, ("metadata", "vs30"), vs30_bounds)
    return df


def filter_by_usable_T(df, usable_T: float|None):
    if usable_T is not None:
            df = df[df[("metadata", "max_usable_T")] >= usable_T]
    return df


def pick_records(df: pd.DataFrame, 
                 scores: np.ndarray,
                 max_count: int,
                 sim_trts: list[str]|None,
                 ok_trt_matches: dict[str, list[str]]|None):
    
    """ selects records from df based on the given scores. 
    Assumed that the lowest score is preferred
    df must have a a multiindex with level 0 headings of 'metadata' and 
    'ims_scaled'. """
    
    # the number of rows in scores should match the length of df
    assert scores.shape[0] == df.shape[0]
    
    empty_row = pd.Series(np.nan, index=df.columns)

    # loop across the columns and find the best fit records not already selected
    selected_records = []
    selected_sses = []
    for ii, c in enumerate(scores.columns): # loop through each sim
        if (sim_trts is not None) and (ok_trt_matches is not None):
            # get the scores only for the same trt type as the sim
            sim_trt = sim_trts[ii]
            allowable_trts = ok_trt_matches[sim_trt] 
            allowable_recs = df[df[("metadata", "trt")].isin(allowable_trts)].index
        else:
            # allow scores from all trt types
            allowable_recs = df.index
            
        col = scores.loc[allowable_recs, c].to_numpy()            # trt filtered, unsorted (original order)
        sorted_col = np.sort(col)     # sorted smallest to largest
        record_selected = False
        for v in sorted_col:    
            # loop through the scores for each record until we find the smallest 
            # score that hasn't been selected already and is a valid selection
            idxs = np.where(col == v)[0]
            for idx in idxs:        # if multiple have the same score take first
                record = df.loc[allowable_recs[idx], :]
                if valid_selection(record, pd.DataFrame(selected_records), max_count):
                    selected_records.append(record)
                    record_selected = True
                    selected_sses.append(v)
                    break 

            if record_selected: # moves onto the next sim
                break

        if not record_selected: # no record found
            selected_records.append(empty_row)

    selected_records = pd.DataFrame(selected_records)
    return selected_records 


def set_weights(weights, n):
    """ the weights should have the same length as the number of selection
    imts """
    if weights is None:
        weights = np.ones(n)
    else:
        assert len(weights) == n
        weights = np.array(weights)                             # shape (M,)
    weights = (weights / weights.sum())                     #  normalise to one
    return weights


def calculate_sse_scores(
        df: pd.DataFrame, 
        simulations: pd.DataFrame, 
        sig_cond: pd.DataFrame,
        weights: list|np.ndarray|None,
        conditioning_imt: IMT):
    """ convenience wrapper for record_sses that takes the whole record df"""
    scaled_im_df = df["ims_scaled"].drop(columns=conditioning_imt.string)
    sse = record_sses(scaled_im_df, simulations, sig_cond, weights)
    return sse 


def record_sses(scaled_rec_df: pd.DataFrame, simulation_df: pd.DataFrame, 
               sig_cond: pd.DataFrame|np.ndarray, weights: np.ndarray
               ) -> pd.DataFrame:
    rdf = scaled_rec_df
    sdf = simulation_df

    recs = np.log(rdf.to_numpy())                              # shape (R,M)
    sims = np.log(sdf.loc[rdf.columns, :].to_numpy())          # shape (M,N)
    if isinstance(sig_cond, pd.DataFrame):
        sigs =  sig_cond.loc[rdf.columns, :].to_numpy()            # shape (M,N)
    else:
        sigs = sig_cond

    # broadcast to three dimensions
    recs = recs[:, :, np.newaxis]              # shape (R,M) -> (R,M,1)
    sims = sims[np.newaxis, :, :]              # shape (R,M) -> (1,M,N)
    sigs = sigs[np.newaxis, :, :]              # shape (R,M) -> (1,M,N)
    weights = weights[np.newaxis, :, np.newaxis]   # shape (M,) -> (1,M,1) 

    # calculate the sse for all records and sims at the same time
    zs = (sims - recs) / sigs
    sse = np.sum(weights * zs ** 2, axis=1)         # shape (R,N)

    sse = pd.DataFrame(sse, index=rdf.index, columns=sdf.columns)

    return sse


def valid_selection(rec_id, selected_records: pd.DataFrame,
                     max_count):
    # get the data of the selected record
    sr = selected_records
    # Checks:
    if len(sr) == 0:    # no records selected yet
        return True
    if (_record_or_other_component_already_selected(rec_id, sr) or      
        _adding_exceeds_max_event_count(rec_id, sr, max_count)):
        return False    # not a valid selection
    return True


def _record_or_other_component_already_selected(rec_id, selected_records):
    # if the index, event_id and event time match then
    # Convert Series to dict for easy column-wise comparison
    cols_to_match = [("metadata", "event_id"), ("metadata", "station_code"), ("metadata", "location_code")] 
    rec = rec_id[cols_to_match]  # get key columns of trial record
    sr = selected_records[cols_to_match] # get key columns of selected records
    # Check if any row in 'sr' matches the values of 'rec'
    is_in_sr = (sr == rec).all(axis=1).any()
    return is_in_sr


def _adding_exceeds_max_event_count(
        rec_id: pd.Series, selected_records: pd.DataFrame, max_count: int):
    event_counts = selected_records.groupby(("metadata", "event_id")).size()

    if rec_id[("metadata","event_id")] not in event_counts.index:
        return False

    if event_counts[rec_id[("metadata","event_id")]] + 1 > max_count:
        return True
    return False


def scale_records_in_db(df:pd.DataFrame, 
                        conditioning_imt: IMT, 
                        conditioning_iml: float):
    """ takes a record database, conditioning_im and imtl and scales all the 
    ims for each record appropriately. 
    the record database must have a multilevel index with level 0 headings of
    'ims' and level 1 headings corresponding to the imt strings. The scaled
    ims are added as new columns under the heading 'ims_scaled'."""

    # df is assumed to be a df containing only columns of ims
    ims = df["ims"]
    
    # get scale factor for every record and add the new data to df
    srdf, alpha = _scale_records_in_db(ims, conditioning_imt, conditioning_iml)
    df = df.copy()
    df.loc[:, ("metadata", "alpha")] = alpha.copy()
    new_cols = pd.MultiIndex.from_product([['ims_scaled'], srdf.columns])
    df = pd.concat([df, srdf.set_axis(new_cols, axis=1)], axis=1)

    return df


def _scale_records_in_db(
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
    sf_cond_imt = conditioning_value / df[conditioning_imt.string] 
    gamma = get_sf_exponent(conditioning_imt.name)
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
       




if __name__ == "__main__":

    ...