""" 
Module with functions that perform df operations common to several modules
"""

import re
import numpy as np
import pandas as pd
from openquake.hazardlib.imt import SA
from scipy.interpolate import make_interp_spline

def interpolate_for_sa(df: pd.DataFrame | pd.Series, T: float | list[float],
                       pga_to_small_T: bool=False):
    
    is_series = isinstance(df, pd.Series)
    working_df = df.to_frame().T if is_series else df
    
    # Parse Periods from column headers
    pattern = re.compile(r"SA\(([\d.]+)\)")
    cols = [c for c in working_df.columns if pattern.search(c)]
    periods = np.array([float(pattern.search(c).group(1)) for c in cols])
    
    # Sort for interpolation safety
    idx = np.argsort(periods)
    periods, cols = periods[idx], np.array(cols)[idx]
    
    T_target = np.atleast_1d(T).astype(float)
    
    # Handle Bounds and PGA logic
    p_min, p_max = periods.min(), periods.max()
    
    if pga_to_small_T:
        # Map 0.0 to 0.001 for both source and target if desired
        periods = np.where(periods == 0, 1e-3, periods)
        T_targets_interp = np.where(T_target == 0, 1e-3, T_target)
        # Update bounds check for the mapped value
        p_min_check = 0.001 if p_min == 0 else p_min
    else:
        T_targets_interp = T_target
        p_min_check = p_min

    # Explicit check for both lower and upper bounds
    out_of_bounds = (T_targets_interp < p_min_check) | (T_targets_interp > p_max)
    if np.any(out_of_bounds):
        invalid_ts = T_targets_interp[out_of_bounds]
        raise ValueError(
            f"Target T {invalid_ts} out of bounds [{p_min}, {p_max}]."
        )

    # Vectorized Interpolation
    log_data = np.log(working_df[cols].to_numpy().T)    # interpolation axis is axis=0
    log_periods = np.log(periods)
    spline = make_interp_spline(log_periods, log_data, k=1)
    
    log_targets = np.log(T_targets_interp)
    res_values = np.exp(spline(log_targets)).T

    # Format Output
    col_names = [f"SA({t})" for t in T_target]
    res_df = pd.DataFrame(res_values, index=working_df.index, columns=col_names)

    if is_series:
        return float(res_df.iloc[0, 0])
    return res_df.iloc[:, 0] if len(T_target) == 1 else res_df


def label_to_str(label) -> str:
    """
    Flatten a MultiIndex column label into a single string like 'U_SA(0.02)'.
    - For tuples, finds the component (U/V/W) and the 'SA(...)' part if present.
    - Joins them with '_'.
    - Otherwise, joins all non-null string parts with '_'.
    """
    if isinstance(label, str):
        return label

    if isinstance(label, tuple):
        # Clean out None/NaN parts
        parts = [p for p in label if isinstance(p, str) and 
                 p.strip() not in ("", "nan")]

        # Otherwise, join all non-empty parts
        if parts:
            return "_".join(parts)

    return str(label)

# Helper to extract mapping period -> column/key
def extract_keys_from_labels(pattern, labels, return_keys=False):
    keys = {}
    for lbl in labels:
        lbl_str = label_to_str(lbl)
        m = pattern.search(lbl_str)
        if m:
            keys[float(m.group(1))] = lbl
    if not keys:
        raise ValueError("No SA columns/labels found in DataFrame")
    periods_sorted = np.array(sorted(keys.keys()))
    cols_sorted = [keys[p] for p in periods_sorted]
    
    if return_keys:
        return keys
    return periods_sorted, cols_sorted
    

def safe_multiindex_rename(df: pd.DataFrame, rename_dict: dict) -> pd.DataFrame:
    """
    Safely rename MultiIndex columns when rename keys have fewer levels than df.columns.
    
    Any rename key shorter than df.columns.nlevels is padded with None until the lengths match.
    """
    nlevels = df.columns.nlevels
    new_mapper = {}

    for key, new_name in rename_dict.items():
        # Ensure tuple
        if not isinstance(key, tuple):
            key = (key,)

        # Pad with None to match MultiIndex depth
        if len(key) < nlevels:
            key = key + (np.nan,) * (nlevels - len(key))

        new_mapper[key] = new_name

    # Perform rename
    return df.rename(columns=new_mapper)



