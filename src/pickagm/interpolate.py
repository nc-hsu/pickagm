import numpy as np
import pandas as pd
from typing import Iterable

Num = int | float
Nums = Num | Iterable[Num]

def _to_arr(x: Nums) -> np.ndarray:
    return np.atleast_1d(np.asarray(x, float))


def _prep(data: np.ndarray,
          axis: int) -> tuple[np.ndarray, tuple[int, ...]]:
    d = np.moveaxis(data, axis, 0)
    return d


def _restore(out: np.ndarray,
             axis: int) -> np.ndarray:
    return np.moveaxis(out, 0, axis)


def _check_dims(d: np.ndarray,
                x: np.ndarray,
                axis: int) -> bool:

    if axis == 0:
        s = "columns"
    elif axis == 1:
        s = "rows"

    if len(x) != d.shape[1]:
        raise ValueError(f"For axis {axis} "
                         f"the number of {s} in data ({d.shape[1]}) " 
                         f"must equal the number of x values ({len(x)})")
    return


def _check_extrapolation(x: np.ndarray,
                         nx: np.ndarray):
    if any(nx > np.max(x)) or any(nx < np.min(x)):
        raise ValueError("Extrapolation is not supported. " \
                         "All new x values must be: min(x) <= new x <= max(x)")

    
def _interp_vec(d: np.ndarray,
                x: np.ndarray,
                nx: np.ndarray) -> np.ndarray:
    idx = np.searchsorted(x, nx, side="right") - 1
    idx = np.clip(idx, 0, len(x) - 2)
    x0 = x[idx]
    x1 = x[idx + 1]
    w = (nx - x0) / (x1 - x0)
    y0 = d[:, idx]
    y1 = d[:, idx + 1]
    return y0 + w * (y1 - y0)


def _merge_arrays(orig_x: np.ndarray,
                  new_x: np.ndarray,
                  orig_y: np.ndarray,
                  new_y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:

    # Combine x-values
    x_all = np.concatenate([orig_x, new_x])
    # True where value is from the original grid
    # is_orig = np.concatenate([
    #     np.ones_like(orig_x, bool),
    #     np.zeros_like(new_x, bool)
    # ])

    # Sort x-values
    order = np.argsort(x_all)
    x_sorted = x_all[order]
    # is_orig_sorted = is_orig[order]

    # Build merged y-values
    # orig_y: (Nrows, n_orig)
    # new_y:  (Nrows, n_new)
    y_all = np.concatenate([orig_y, new_y], axis=1)
    y_sorted = np.take_along_axis(y_all, order[None, :], axis=1)

    # For duplicates: original overwrites interpolated
    _, idx_last = np.unique(x_sorted, return_index=True)
    keep = np.zeros_like(x_sorted, bool)
    keep[idx_last] = True

    # Ensure original data is retained when duplicated
    # original entries override interpolated ones
    out_x = x_sorted[keep]
    out_y = y_sorted[:, keep]

    return out_x, out_y


def interp_linear_xy(d: np.ndarray,
                     x: Iterable[Num],
                     nx: Nums,
                     axis: int = 1,
                     return_mode: str = "interp"):

    y_new = _interp_vec(d, x, nx)

    if return_mode == "interp":
        return nx, _restore(y_new, axis)

    # merged mode
    x_m, y_m = _merge_arrays(x, nx, d, y_new)
    return x_m, _restore(y_m, axis)


def interp_linear_x_logy(d: np.ndarray,
                         x: Iterable[Num],
                         nx: Nums,
                         axis: int = 1,
                         return_mode: str = "interp"):
    if np.any(d <= 0):
        raise ValueError("log(y) needs y>0.")
    logy = np.log(d)
    y_new = np.exp(_interp_vec(logy, x, nx))

    if return_mode == "interp":
        return nx, _restore(y_new, axis)

    x_m, y_m = _merge_arrays(x, nx, d, y_new)
    return x_m, _restore(y_m, axis)


def interp_logx_logy(d: np.ndarray,
                     x: np.array,
                     nx: np.array,
                     axis: int = 1,
                     return_mode: str = "interp"):
    if np.any(d <= 0):
        raise ValueError("log(y) needs y>0.")
    if np.any(x <= 0) or np.any(nx <= 0):
        raise ValueError("log(x) needs x>0.")
    lx = np.log(x)
    lnx = np.log(nx)
    logy = np.log(d)
    y_new = np.exp(_interp_vec(logy, lx, lnx))

    if return_mode == "interp":
        return nx, _restore(y_new, axis)

    x_m, y_m = _merge_arrays(x, nx, d, y_new)
    return x_m, _restore(y_m, axis)


def interpolate(data: np.ndarray,
                x_vals: Iterable[Num],
                new_x_vals: Nums,
                scheme: str = "linear_xy",
                axis: int = 0,
                return_mode: str = "interp") -> np.ndarray:
    """Interpolates columns of an m x n 2D numpy array

    Args:
        data (np.ndarray): m x n array of data to interpolate
        x_vals (Iterable[Num]): x values corresponding to the columns
        new_x_vals (Nums): new x values that we want to find the columns for
        scheme (str, optional): The type of interpolation to perform. One of:
                - "linear_xy"
                - "linear_x_logy"
                - "logx_logy"
                Defaults to "linear_xy".
        axis (int): the axis along which the interpolation is performed.
            axis = 0 is interpolation along a row (i.e a row is the y-values)
            axis = 1 is interpolation down a column (i.e. c column is the y-vals)
            Defaults to 0
        return_mode (str): either "interp" or "merged". "interp just reutrns 
            the interpolated data. "merged" returns the interpolated data merged
            with the original data and sorted accoring to the x values.

    Raises:
        ValueError: Unknown interpolation scheme is given

    Returns:
        np.ndarray: the interpolated columns in the order of the x_values 
        provided
    """
    x = _to_arr(x_vals)
    nx = _to_arr(new_x_vals)
    d = _prep(data, axis)
    
    _check_dims(d, x, axis)
    _check_extrapolation(x, nx)

    if scheme == "linear_xy":
        return interp_linear_xy(d, x, nx, axis, return_mode)
    if scheme == "linear_x_logy":
        return interp_linear_x_logy(d, x, nx, axis, return_mode)
    if scheme == "logx_logy":
        return interp_logx_logy(d, x, nx, axis, return_mode)
    raise ValueError(f"Unknown scheme: {scheme}")


def interpolate_SA_df(
        df: pd.DataFrame, 
        new_periods: list[float], 
        scheme: str, 
        return_mode: str) -> pd.DataFrame:
    # df has a 3-level multiindex for the columns. 
    # (component, IM (e.g. PGA, SA, RSD595), period)
    # assumes that df corresponds to a single component e.g. "rotD50" or "U"
    # interpolation between PGA and the lowest SA(T) is supported if PGA is 
    # present in the data frame
    # returns the original dataframe with the interpolated SA values inserted
    # the dataframe is then sorted lexicographically
    
    idx = pd.IndexSlice
    component = df.columns.get_level_values(0)[0]

    # check if PGA is included
    im_types = df.columns.get_level_values(1)
    if "PGA" in df.columns.get_level_values(1):
        data_array = df.loc[:, idx[:, ["PGA", "SA"], :]].to_numpy()
        periods = list(df.loc[:, idx[:, ["PGA", "SA"], :]].columns.get_level_values(2))
        periods[0] = 0
        mask = ~im_types.isin(["PGA", "SA"])
        rest_df = df.loc[:, mask]

    else:
        data_array = df.loc[:, idx[:, "SA", :]].to_numpy()
        periods = df.loc[:, idx[:, "SA", :]].columns.get_level_values(2)
        mask = ~im_types.isin(["SA"])
        rest_df = df.loc[:, mask]
    
    new_periods, new_data = interpolate(
        data_array, periods, 
        new_periods, scheme=scheme, 
        axis=0, return_mode=return_mode)
    
    new_column_tuples = [(component, "SA", p) if p > 0 
                         else (component, "PGA", 0) for p in new_periods]
    new_columns = pd.MultiIndex.from_tuples(
        new_column_tuples, names=df.columns.names)
    new_df = pd.DataFrame(new_data, columns=new_columns, index=df.index)

    # combine new_df and rest_df
    combined_df = pd.concat([new_df, rest_df], axis=1)
    combined_df.reindex(columns=sorted(combined_df.columns))

    return combined_df
