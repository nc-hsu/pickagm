from typing import Protocol

import numpy as np
import pandas as pd


class SiteRupContextBuilder(Protocol):

    def __init__(self, n_ctxs: int, site_params: dict, *args, **kwargs) -> None:
        ...

    def ctx(self, rup_params: dict) -> np.recarray:
        ...
    ...

def get_rup_ctxs(disagg_type: str, rup_samples: pd.DataFrame) -> list[dict]:
    
    disagg_mappping = {
        "Mag": "mag",
        "Dist": "rjb",
        "Lon": "lon",
        "Lat": "lat",
        "TRT": "trt",
        "Eps": "eps"
    }
    
    disagg_params = disagg_type.split("_")
    for param in disagg_params:
        if param not in disagg_mappping:
            raise ValueError(f"Disagg type {disagg_type} is not supported")
    
    rup_ctxs = []
    for s in rup_samples[disagg_params].itertuples(index=False):
        rup_ctx = {}
        for param in disagg_params:
            rup_ctx[disagg_mappping[param]] = getattr(s, param)
        
        rup_ctxs.append(rup_ctx)
    
    return rup_ctxs


def get_site_rup_ctxs(disagg_type: str, 
                      disagg_dist: pd.DataFrame, 
                      ctx_builder: SiteRupContextBuilder):
    
    rup_ctxs = get_rup_ctxs(disagg_type, disagg_dist) # do them all!
    site_rup_ctxs = [ctx_builder.ctx(rup_ctx) for rup_ctx in rup_ctxs]
    site_rup_ctxs = np.concatenate(site_rup_ctxs).view(np.recarray)
    return site_rup_ctxs

