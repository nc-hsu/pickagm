"""
Module for stuff related to accessing and working with the different earthquake
databases that are available
"""
import re
from typing import Callable
import importlib.resources as res

import numpy as np
import pandas as pd
from openquake.hazardlib.imt import IMT, SA, PGA, PGV, RSD595, CAV, IA

from pickagm.dfops import (
    interpolate_for_sa, safe_multiindex_rename, label_to_str)


esm_nonSA_im_translations = {
    "pga": str(PGA()),
    "pgv": str(PGV()),
    "T90": str(RSD595()),
    "CAV": str(CAV()),
    "ia": str(IA()),
    "pgd": "PGD"
}

esm_unit_conversions = {
    str(PGA()): 1/981.0,
    "SA": 1/981.0,
    "AvgSA": 1/981.0
}

def esm_multilevel_im_translations() -> dict:
    components =  ["U", "V", "W"]
    translation_dict = {}
    for component in components:
        for im, translation  in esm_nonSA_im_translations.items():
            translation_dict[(component, im)] = (component, translation)
    return translation_dict


def load_flatfile(filename: str, n_rows: int) -> pd.DataFrame:
    with res.files("pickagm.data").joinpath(f"{filename}.csv") \
        .open("r", encoding="utf-8") as f:
        df = pd.read_csv(f, sep=",", header=[0], nrows=n_rows)
    return df


def load_esm_flatfile(n_rows: int) -> pd.DataFrame:
    df = load_flatfile("ESM_flatfile_SA", n_rows)
    return df


def filter_imt_data_esm(flatfile: pd.DataFrame, imts: list[IMT]):
    """ 
    returns a dataframe with only the im data need for the selection based
    on the imts in imts. If the database does not contain values for certain
    imt (and they cannot be interpolated i.e. SA) then they are currently 
    ignored and no warning is thrown. 
    """
    df = flatfile    
    df = rename_im_columns(
        df, rename_sa_columns_esm, esm_multilevel_im_translations)
    df = get_imt_df(df, ["U", "V", "W"], imts)
    df = scale_columns_for_units(df, esm_unit_conversions)
    return df
 

def rename_im_columns(
        df: pd.DataFrame, 
        sa_func: Callable[[pd.DataFrame], pd.DataFrame], 
        imt_translation_func: Callable[[], dict]
        ) -> pd.DataFrame:
    
    df = sa_func(df)
    df.columns = df.columns.str.split('_', expand=True)
    df = safe_multiindex_rename(df, imt_translation_func())
    return df
    

def rename_sa_columns_esm(df: pd.DataFrame) -> pd.DataFrame:
    """
    Rename columns matching pattern X_TYY_YYY -> X_SA(YY.YYY)
    where X ∈ {U, V, W} and YY_YYY represents a float value.
    The numeric part is parsed as a float for consistent formatting.
    """
    def rename_col(col: str) -> str:
        match = re.fullmatch(r'(\w+)_T(\d{1,2})_(\d+)', col)
        if match:
            comp, int_part, dec_part = match.groups()
            period = float(f"{int_part}.{dec_part}")
            return f"{comp}_{SA(period)}"
        return col  # leave unchanged if no match

    return df.rename(columns=rename_col)
  

def get_imt_df(
        df: pd.DataFrame, gm_component_tags: list[str], imts: list[IMT]
        ) -> pd.DataFrame:
    """ 
    returns a dataframe only containing the imt values for each record that
    we need for the selection
    """
    dfs = []
    
    for component in gm_component_tags:
        sub_df = df[component]
        ims_in_df = [str(imt) for imt in imts if str(imt) in sub_df.columns]
        sa_periods_to_interp = [imt.period for imt in imts if 
                         str(imt).startswith("SA") and 
                         not (str(imt) in sub_df.columns)]
        
        interped_ims = interpolate_for_sa(sub_df, sa_periods_to_interp)
        sub_df = sub_df[ims_in_df]
        sub_df.columns = [label_to_str(label) for label in sub_df.columns]
        new_sub_df = pd.merge(sub_df, interped_ims,
                              left_index=True, right_index=True)
        new_sub_df = new_sub_df[sorted(new_sub_df.columns)]
        columns = [(component, c) for c in new_sub_df.columns]
        new_sub_df.columns = pd.MultiIndex.from_tuples(columns)
        dfs.append(new_sub_df)

    # Source - https://stackoverflow.com/questions/38089010/merge-a-list-of-pandas-dataframes
    # Posted by cottontail
    # Retrieved 5.11.2025, License - CC-BY-SA 4.0
    im_df = dfs[0]
    for d in dfs[1:]:
        im_df = im_df.merge(d, left_index=True, right_index=True)

    return im_df


def scale_columns_for_units(
        df: pd.DataFrame, conversion_map: dict) -> pd.DataFrame:
    for tag, sf in conversion_map.items():
        mask = df.columns.get_level_values(1).str.startswith(tag)
        df.loc[:,mask] *= sf
    return df


def extract_im_data(df: pd.DataFrame, idx_im_data: int):
    im_df = df[df.columns[idx_im_data:]]
    return im_df.copy()


def extract_metadata(flatfile_df: pd.DataFrame, idx_im_data: int):
    return flatfile_df[flatfile_df.columns[:idx_im_data]].copy()


def expand_SA_column_names(col_names: list[str]):
    pattern = re.compile(rf"SA\(([\d.]+)\)")
    new_names = []
    for c in col_names:
        cs = c.split("_")
        m = pattern.search(c[1])
        if m: t = float(m.group(1))
        else: t = "None"
        new_names.append((cs[0], "SA", t))
    return new_names
