import json
from pathlib import Path
from typing import Protocol

import h5py
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

CORR_COEFFS_FOLDER = Path(__file__).parent / "im_corr_data"

## This code is based on the correlation models in OpenQuake, but it is adapted
## to allow consider all IMs instead of just SA.
class CorrelationModel(Protocol):

    rho: dict[str, pd.DataFrame] # dict of correlation matrices for each residual type
    model: str

    def get_correlation(self, imt1: tuple, imt2: tuple, residual_type) -> float:
        ...

    def build_correlation_matrices(self, 
                                   incl_nonSA_imts: list[str]|None=None,
                                   incl_SA_Ts: list[float]|None=None,
                                   excl_SA_T_gt: float|None=None,
                                   excl_SA_T_lt: float|None=None) -> dict[str, pd.DataFrame]:
        ...


class EmpiricalCorrelationModel:
    """
    Intermediate class to handle the interpolation of empirical correlation 
    data and the construction of the correlation matrices.
    """
    def __init__(self, rho_arrays, rho_imts, 
                 incl_nonSA_imts: list[str]|None=None,
                 incl_SA_Ts: list[float]|None=None):
        """ 
        # TODO: add more details to docstring      
        :param np.ndarray avg_periods : 
            The periods used for calculating the spectral acceleration
        
        :param dict[str, np.ndarray] rho_arrays: 
            The correlation data. The keys rho_arrays are strings that represent
            the type of residual that the associated array corresponds to. 
        
        :param np.ndarray rho_periods: 
            A 1D array of the spectral periods that correspond to the rows and 
            columns of the correlation data

        Notes:
        1)  It is assumed that the same periods apply to all correlation 
            matrices in the same correlation model.
        
        2)  The keys for rho_arrays follow the notation of Al Atik et al (2010).
            Each key maps to a method that returns the appropriate
            correlation coefficient. Currently the this is:
            "total"  -> get_correlation()
            "dB_e"   -> get_between_event_correlation()
            "dS2S_s" -> get_between_site_correlation()
            "dWS_es" -> get_within_event_correlation()
            The class can be instantiated using a dict that doesn't contain all
            the keys, however if the mapped function is called then a ValueError
            will be raised. The provided dict must have as a minimum the key
            "total".
        """
        # rho_arrays is dictionary of the arrays of corr coeffs that are part
        # of the correlation model. The keys of this dictionary indicate what
        # correlation functions can be called without raising an error. 

        if not isinstance(rho_arrays, dict):
            raise TypeError("rho_arrays must be of type dict")
        if "total" not in rho_arrays.keys():
            raise ValueError("rho_arrays must contain the key 'total'")
        
        # store the raw correlation data
        self.raw_rhos = rho_arrays
        self.raw_imts = rho_imts
        self.raw_Ts = np.array([imt[1] for ii, imt in enumerate(self.raw_imts) 
                                if imt[0] == "SA"]) 
        self.valid_residual_types = rho_arrays.keys() 

        MIN_PERIOD = min(self.raw_Ts)
        MAX_PERIOD = max(self.raw_Ts)

        # set some attributes for the requested IMTs and periods.
        # self.included_nonSA_imts = incl_nonSA_imts
        # self.desired_SA_Ts = incl_SA_Ts
        # self.exclude_SA_T_gt = exclude_SA_T_gt
        # self.exclude_SA_T_lt = exclude_SA_T_lt
        self._nonsa_idxs = [ii for ii, imt in enumerate(self.raw_imts) 
                           if not imt[0] == "SA"]
        self._sa_idxs = [ii for ii, imt in enumerate(self.raw_imts) 
                        if imt[0] == "SA"]

        # create the interpolation functions for each residual type / array
        self.sa_interpers = {}

        # Do some checks to make sure the requested IMTs are in the model or 
        # can be interpolated from the available data.
        
        # check that non-SA IMTs are in the model
        if incl_nonSA_imts is not None:
            invalid_imts = [imt for imt in incl_nonSA_imts 
                            if imt not in [ii[0] for ii in self.raw_imts]]
            if invalid_imts != []: 
                    raise ValueError(f"Requested IMT(s) {invalid_imts} "
                                     f"are not in the model")
            

        # check that SA IMTs are in the model or can be interpolated from the
        # available data
        if incl_SA_Ts is not None: 
            if any(incl_SA_Ts < MIN_PERIOD) and ["PGA", "None"] not in rho_imts: # then we can't interpolate -> raise error
                raise ValueError(f"Period ({min(incl_SA_Ts):.3f}) is less than "
                                 f"the minimum allowable period for the "
                                 f"correlation model ({MIN_PERIOD:.3f}).")
            if any(incl_SA_Ts > MAX_PERIOD):
                raise ValueError(f"Period ({max(incl_SA_Ts):.3f}) is greater than "
                                 f"the maximum allowable period for the "
                                 f"correlation model ({MAX_PERIOD:.3f}).")

        # if we get to this point then the model is valid and we can build the 
        # desired correlation matrices for the requested IMTs.  

        self.rho = self.build_correlation_matrices(
            incl_nonSA_imts, incl_SA_Ts)


    def build_correlation_matrices(self, 
                                   incl_nonSA_imts: list[str]|None=None,
                                   incl_SA_Ts: list[float]|None=None):
        """
        Interpolate the raw correlation data for required periods (avg_periods).
        Returns a dict with keys indicating the type of residual and the values 
        being 2D numpy arrays of interpolated corr coeffs.
        """
        # interpolates the raw coeff data for all the provided matrices and 
        # returns a dictionary with the keys indicating the type of residual
        # and the values being 2D numpy arrays of interpolated corr coeffs.

        corr_mats = {}
        for residual_type in self.valid_residual_types:
            corr_mats[residual_type] = self._build_correlation_matrix(
                residual_type, incl_nonSA_imts, incl_SA_Ts)
        return corr_mats
    

    def _build_correlation_matrix(self, 
                                  residual_type,
                                  incl_nonSA_imts,
                                  incl_SA_Ts):
        
        # get the raw data for the specified residual type
        raw_rhos = self.raw_rhos[residual_type]
        
        #? if the period 0.0s is included we need to get the PGA
        # get the desired nonSA imts
        nonsa_idx = []
        if incl_nonSA_imts is not None:
            nonsa_idx += [ii for ii, imt in enumerate(self.raw_imts) 
                        if imt[0] in incl_nonSA_imts]
            
        elif ((any(incl_SA_Ts < min(self.raw_Ts)) and (incl_nonSA_imts is None)) or
              (any(incl_SA_Ts < min(self.raw_Ts)) and ("PGA" not in incl_nonSA_imts))): 
            # we need to include PGA for interpolation
            incl_nonSA_imts = ["PGA"]
            nonsa_idx += [ii for ii, imt in enumerate(self.raw_imts) 
                          if imt[0] == "PGA"]
        else:
            nonsa_idx = self._nonsa_idxs # the non SA imts will be removed later

        # get the SA-SA correlations for the desired periods
        if incl_SA_Ts is None:
            periods_to_interp = self.raw_Ts 
        
        elif all(incl_SA_Ts >= min(self.raw_Ts)):
            periods_to_interp = incl_SA_Ts
            idxs_for_sa_data = self._sa_idxs
            original_periods = self.raw_Ts
        
        else: 
            # Remove the 0 period to interp because we already have the PGA
            periods_to_interp = np.array(incl_SA_Ts)
            periods_to_interp = periods_to_interp[periods_to_interp != 0]
            idx_pga = [ii for ii, imt in enumerate(self.raw_imts) 
                                  if imt[0] == "PGA"]
            # Add the PGA to the orignal data
            idxs_for_sa_data = idx_pga + self._sa_idxs
            original_periods = np.append(np.array(0), self.raw_Ts)
    
        interper = self._create_SA_interper(
            #? PGA needs to be included here if a period less than the smallest 
            #? SA(T) in the correlation data
            raw_rhos[np.ix_(idxs_for_sa_data, idxs_for_sa_data)], 
            original_periods, 
            residual_type)
        rho_SA_SA = self._interpolate_matrix(interper, periods_to_interp)

        # get the nonSA-SA correlations for the desired periods
        rho_nonSA_SA = []
        for row in raw_rhos[nonsa_idx,:]:
            rho_nonSA_SA.append(
                np.interp(periods_to_interp, original_periods, row[idxs_for_sa_data]))
        rho_nonSA_SA = np.array(rho_nonSA_SA) # shape (n_nonSA, n_SA)
        
        # get the nonSA-nonSA correlations (no interpolation necessary)
        rho_nonSA_nonSA = raw_rhos[np.ix_(nonsa_idx, nonsa_idx)]

        # combine the nonSA-nonSA, nonSA-SA and SA-SA correlations into a single matrix
        n_nonSA = len(nonsa_idx)
        n_SA = len(periods_to_interp)
        rho_mat = np.zeros((n_nonSA + n_SA, n_nonSA + n_SA))
        rho_mat[:n_nonSA, :n_nonSA] = rho_nonSA_nonSA
        rho_mat[:n_nonSA, n_nonSA:] = rho_nonSA_SA
        rho_mat[n_nonSA:, :n_nonSA] = rho_nonSA_SA.T
        rho_mat[n_nonSA:, n_nonSA:] = rho_SA_SA

        # create a dataframe with the appropriate row and column labels
        imt_labels = [imt[0] for ii, imt in enumerate(self.raw_imts) 
                      if ii in nonsa_idx] +\
                     [f"SA({np.round(T,3)})" for T in periods_to_interp]
        corr_df = pd.DataFrame(rho_mat, 
                               index=imt_labels, 
                               columns=imt_labels)
        
        if incl_nonSA_imts is None: # remove nonSAimts
            drop_imts = [imt for ii, imt in enumerate(corr_df.index) if ii in nonsa_idx]
            corr_df = corr_df.drop(columns=drop_imts, index=drop_imts)

        return corr_df


    def _create_SA_interper(self, corr_data, periods, residual_type):
        interper = RegularGridInterpolator(
                points=(periods, periods),
                values=corr_data,
                method="linear"
                )
        self.sa_interpers[residual_type] = interper
        return interper


    def _interpolate_matrix(self, interper, new_Ts):
        """
        2D interpolation of the correlation matrix.
        Returns a symmetric nxn matrix with 1s on the diagonal
        
        :param RegularGridInterpolator interper: 
            interpolation function to be used
        
        :param np.ndarray new_Ts:
            array of periods at which the matrix should be interpolated
        """
        # create an array of interp points with shape (N*N, 2)
        t1, t2 = np.meshgrid(new_Ts, new_Ts, indexing="ij") # t1 and t2 have shape (N, N)
        interp_pts = np.column_stack([t1.ravel(), t2.ravel()])
        interped_rhos = interper(interp_pts).reshape(len(new_Ts), -1)

        # make sure the diagonal terms are equal to 1.0
        np.fill_diagonal(interped_rhos, 1.0)

        return interped_rhos

    
    def _check_valid_residual_type(self, residual_type):
        if not residual_type in self.valid_residual_types:
            raise ValueError(f"Invalid type of residual ({residual_type}) for "
                             f"correlation model.")
    
    def get_correlation(self, imt1: tuple, imt2: tuple, residual_type) -> float:
        """
        gets the correlation coefficient for the specified imts and residual type.
        Doesn't allow for interpolation, so the imts must be in the calculated
        correlation matrix already.
        # TODO: update docstring to be more specific about the imt format
        
        :param float t1:
            First period of interest.

        :param float t2:
            Second period of interest.

        :return float rho:
            The estimated correlation coefficient.
        """
        # check if the residual type is valid for this model
        self._check_valid_residual_type(residual_type)

        rho_df = self.corr_mats[residual_type]

        # check that the imts are valid
        self._check_valid_imt(imt1, rho_df)
        self._check_valid_imt(imt2, rho_df)

        return rho_df.loc[imt1, imt2]
    
    def _check_valid_imt(self, imt, corr_df):
        if imt not in corr_df.index:
            raise ValueError(f"IMT {imt} is not in the correlation matrix.")


    def __str__(self):
        return f"{self.model}"
    

def load_corr_coeff_data_from_hdf5(fp, group_name):
    """
    Loads an HDF5 group into a dictionary of numpy arrays.
    """
    with h5py.File(fp, "r") as f:
        if group_name not in f:
            raise KeyError(f"Group '{group_name}' not found.")
        
        group = f[group_name]
        
        # Dictionary comprehension: 
        # item[()] loads the dataset into a numpy array
        rho_dict = {
            key: item[()] 
            for key, item in group.items() 
            if isinstance(item, h5py.Dataset)
        }

        imt_dict = {
            key: json.loads(item.attrs["imts"])
            for key, item in group.items()
            if isinstance(item, h5py.Dataset)
        }

    return rho_dict, list(imt_dict.values())[0] # assuming the same imts for all datasets in the group


class ClemettCorrelationModelAsc(EmpiricalCorrelationModel):
    """
    Correlation model for the active shallow crust regions in Europe
    Clemett and Gündel (2026) - Paper in preparation
    """
    def __init__(self, 
                 incl_nonSA_imts: list[str]|None=None,
                 incl_SA_Ts: list[float]|None=None):
        
        self.model = "Clemett26ASC"

        fp = CORR_COEFFS_FOLDER / "clemett_corr_coeffs.hdf5"
        rho_arrays, rho_periods = load_corr_coeff_data_from_hdf5(fp, "asc")
        super().__init__(rho_arrays, rho_periods, 
                         incl_nonSA_imts, incl_SA_Ts) 
        

class ClemettCorrelationModelSInter(EmpiricalCorrelationModel):
    """
    Correlation model for subduction interface regions in Europe
    Clemett and Gündel (2026) - Paper in preparation
    """
    def __init__(self, 
                 incl_nonSA_imts: list[str]|None=None,
                 incl_SA_Ts: list[float]|None=None):
        
        self.model = "Clemett26SInter"

        fp = CORR_COEFFS_FOLDER / "clemett_corr_coeffs.hdf5"
        rho_arrays, rho_periods = load_corr_coeff_data_from_hdf5(fp, "sinter")
        super().__init__(rho_arrays, rho_periods, 
                         incl_nonSA_imts, incl_SA_Ts)


class ClemettCorrelationModelSSlab(EmpiricalCorrelationModel):
    """
    Correlation model for the subduction slab regions in Europe
    Clemett and Gündel (2026) - Paper in preparation
    """
    def __init__(self, 
                 incl_nonSA_imts: list[str]|None=None,
                 incl_SA_Ts: list[float]|None=None):
        
        self.model = "Clemett26SSlab"

        fp = CORR_COEFFS_FOLDER / "clemett_corr_coeffs.hdf5"
        rho_arrays, rho_periods = load_corr_coeff_data_from_hdf5(fp, "sslab")
        super().__init__(rho_arrays, rho_periods, 
                         incl_nonSA_imts, incl_SA_Ts)


class ClemettCorrelationModelVrancea(EmpiricalCorrelationModel):
    """
    Correlation model for the Vrancea deep seismicity region in Europe
    Clemett and Gündel (2026) - Paper in preparation
    """
    def __init__(self, 
                 incl_nonSA_imts: list[str]|None=None,
                 incl_SA_Ts: list[float]|None=None):
        
        self.model = "Clemett26Vrancea"

        fp = CORR_COEFFS_FOLDER / "clemett_corr_coeffs.hdf5"
        rho_arrays, rho_periods = load_corr_coeff_data_from_hdf5(fp, "vrancea")
        super().__init__(rho_arrays, rho_periods, 
                         incl_nonSA_imts, incl_SA_Ts)


def conditional_correlation_matrix(
        corr_df: pd.DataFrame, conditioning_imt: str) -> pd.DataFrame:
    
    # check if conditioning imt is in the correlation matrix
    if conditioning_imt not in corr_df.index:
        raise ValueError(f"Conditioning IMT {conditioning_imt} is not in the "
                         f"correlation matrix.")

    # convert to numpy array for calculations / remove conditionin imt
    rho = corr_df.drop(index=conditioning_imt, columns=conditioning_imt) \
                 .to_numpy() 
    
    # get vector of conditional correlations between imt and conditioning imt
    rho_v = corr_df.loc[conditioning_imt, :].drop(conditioning_imt) \
                                            .to_numpy() \
                                            .reshape(-1,1)

    numerator = rho - rho_v @ rho_v.T
    denominator = np.sqrt((1 - rho_v ** 2)) @ np.sqrt((1 - rho_v ** 2)).T
    rho_cond = numerator / denominator

    # force 1.0 on diagonal (should be the case already but just to be safe)
    np.fill_diagonal(rho_cond, 1.0)

    # convert back to dataframe with appropriate labels
    new_index = corr_df.index.drop(conditioning_imt)
    corr_cond_df = pd.DataFrame(rho_cond, index=new_index, columns=new_index)

    return corr_cond_df


# register the correlation models in a dictionary for easy access
CORR_MODELS = {
    "clemett_asc": ClemettCorrelationModelAsc,
    "clemett_sinter": ClemettCorrelationModelSInter,
    "clemett_sslab": ClemettCorrelationModelSSlab,
    "clemett_vrancea": ClemettCorrelationModelVrancea,
}


if __name__ == "__main__":
    from openquake.hazardlib.imt import PGA, RSD595, SA, AvgSA
    sa_periods = np.round(np.logspace(np.log10(0.025), np.log10(8), num=30), 3)
    nonSA_imts = [AvgSA([0,3]).string, RSD595().string, PGA().string]  # strings match the correlation matrix
    # cond_imt = AvgSA([0,3])  # strings match the correlation matrix
    # selection_imts = [RSD595(), PGA()] + [SA(period) for period in sa_periods]  # these need to be imt instances to work with GMMs
    
    corr_model = ClemettCorrelationModelAsc(
        incl_nonSA_imts=nonSA_imts, 
        incl_SA_Ts=sa_periods,
        )
    
    # rho_cond = conditional_correlation_matrix(corr_model.rho["total"], cond_imt.string)
    pass