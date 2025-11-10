""" 
Module with functions that perform df operations common to several modules
"""

import re
import numpy as np
import pandas as pd
from openquake.hazardlib.imt import SA

def interpolate_for_sa(
        df: pd.DataFrame | pd.Series, T: list[float]
        ) -> pd.DataFrame|pd.Series|float:
    """
    Interpolate values for spectral accelerations.

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
    pattern = re.compile(rf"SA\(([\d.]+)\)")

    # Helper to extract mapping period -> column/key
    def _extract_keys_from_labels(labels):
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
        return periods_sorted, cols_sorted

    T_arr = np.atleast_1d(T)

    # Case: Series input -> return a single float (require scalar T)
    if isinstance(df, pd.Series):
        if T_arr.size != 1:
            raise ValueError("When input is a Series, T must be a scalar value")
        periods, keys = _extract_keys_from_labels(df.index)
        values = np.array([float(df[k]) for k in keys])
        log_Ts = np.log(periods)
        log_data = np.log(values)

        interp_val = np.exp(np.interp(
            np.log(float(T_arr[0])), log_Ts, log_data, 
            left=np.nan, right=np.nan))
        if np.isnan(interp_val):
            raise ValueError(f"Could not interpolate SA at T={T_arr[0]}")
        return float(interp_val)

    # Case: DataFrame input -> return Series (for scalar T) or DataFrame 
    # (for multiple T)
    if isinstance(df, pd.DataFrame):
        periods, cols = _extract_keys_from_labels(df.columns)
        data = df[cols].to_numpy()  # shape: (n_rows, n_periods)
        log_data = np.log(data)
        log_Ts = np.log(periods)
        
        results = {}
        for t in T_arr:
            # interpolate each row along the SA-period axis. Log interpolation
            interp_vals = np.exp(np.apply_along_axis(
                lambda row: np.interp(
                    np.log(float(t)), log_Ts, row, left=np.nan, right=np.nan),
                axis=1,
                arr=log_data
            ))
            col_name = str(SA(float(t)))
            if any(np.isnan(interp_vals)):
                raise ValueError(f"Could not interpolate SA at T={float(t)}")
        
            results[col_name] = interp_vals
            
        res_df = pd.DataFrame(results, index=df.index).map(float)
        return res_df.iloc[:, 0] if res_df.shape[1] == 1 else res_df

    raise TypeError("df must be a pandas.DataFrame or pandas.Series")


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

        # Extract component (U, V, W) and SA(...) if they exist
        comp = next((p for p in parts if p in ("U", "V", "W")), None)
        sa_part = next((p for p in parts if "SA(" in p), None)

        # If both found, return e.g. "U_SA(0.02)"
        if comp and sa_part:
            return f"{comp}_{sa_part}"

        # Otherwise, join all non-empty parts
        if parts:
            return "_".join(parts)

    return str(label)


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
