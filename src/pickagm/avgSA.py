import numpy as np
from openquake.hazardlib.imt import IMT, SA, PGA
from openquake.hazardlib.gsim.base import GMPE

def compute_ln_AvgSA(
        ln_SA_data: np.ndarray, 
        weights: np.ndarray|None=None
        ) -> float|np.ndarray:
    # ln_SA_data: is a 1D or 2D numpy array containing the natural log of the 
    # SA values required to calculate the AvgSA. If 2D the average is
    # calculated along a row. i.e. each column represents a different period. If
    # it is 1D then the average is performed on all values in the vector
    
    # weights: is a 1D array containing the weights associated with each SA(T1)
    # returns ln_AvgSA
    
    if len(ln_SA_data.shape) == 1:
        n_values = len(ln_SA_data)
    elif len(ln_SA_data.shape) == 2:
        n_values = ln_SA_data.shape[1] # number of columns is the number of periods

    if not weights is None:
        if len(weights) != n_values:
            raise ValueError(f"length of 'weights' ({len(weights)}) does not " 
                             f"match 'n' ({n_values})")
        
        raise NotImplementedError("weighting of periods is not yet implemented")

    ln_avgSA = (1 / n_values) * ln_SA_data.sum(axis=1)
    return ln_avgSA


def compute_std_dev_ln_AvgSA(
        std_dev_ln_SA: np.ndarray,
        correlation_matrix: np.ndarray,
        weights: np.ndarray|None=None) -> np.ndarray:
    # std_dev_ln_SA:  is a 1D or 2D numpy array containing the std. dev of
    # the natural log of the SA values used to calculate the AvgSA. 
    # std_dev_ln_SA: can be either the total, between-event or within-event
    # std. deviation.
    # If 2D each column represents a different period and each row a different
    # rupture. If 1D then each values represent the values associated with a 
    # single rupture
    
    # Calculation according to Kohrangi et al. (2018) "Ground-motion models for 
    # average spectral acceleration in a period range: direct and indirect 
    # methods", in Bulletin of Earthquake Engineering, Eq. 3

    if len(std_dev_ln_SA.shape) == 1:
        n_values = len(std_dev_ln_SA)
    elif len(std_dev_ln_SA.shape) == 2:
        n_values = std_dev_ln_SA.shape[1] # number of cols is the no. of periods
    else:
        raise ValueError(f"Dimension = {len(std_dev_ln_SA.shape)}. Must be 1 or 2")

    if not weights is None:
        if len(weights) != n_values:
            raise ValueError(f"length of 'weights' ({len(weights)}) does not " 
                             f"match 'n' ({n_values})")
        
        raise NotImplementedError("weighting of periods is not yet implemented")

    # create an nxmxm data cube where each layer contains the products of 
    # sigma_i, sigma_j, and rho for that rupture. Uses numpy broadcasting
    # sigma_ln_AvgSA[:, :, np.newaxis] has shape (n, m, 1)
    # sigma_ln_AvgSA[:, np.newaxis, :] has shape (n, 1, m)
    ssr = std_dev_ln_SA[:, :, np.newaxis] * \
          std_dev_ln_SA[:, np.newaxis, :] * \
          correlation_matrix
    
    std_dev_ln_SA = np.sqrt((1 / n_values) ** 2 * np.sum(ssr, axis=(1, 2)))
    return std_dev_ln_SA


class indirect_AvgSA_GMPE():
    def __init__(
            self,
            gmm: GMPE|list[GMPE],
            rho_total: np.ndarray|None = None,
            rho_deltaB: np.ndarray|None = None,
            rho_deltaW: np.ndarray|None = None,
            gmm_weights: np.ndarray|None = None,
            mean_only: bool = False,
            avg_periods: list[float]|None=None,):
        
        
        # either rho_total needs to be specified or 
        # both rho_deltaB and rho_deltaW need to be specified
        if not mean_only:
            if (rho_deltaB is not None and rho_deltaW is not None):
                self.use_dB_dW = True
            else:
                self.use_dB_dW = False
            
            if rho_total is None and not self.use_dB_dW: 
                    raise ValueError("Either 'rho_total' or both 'rho_deltaB' and "
                                    "'rho_deltaW must not be None")

        if not isinstance(gmm, list):
            gmm = [gmm]

        if gmm_weights:
            raise NotImplementedError("weighting of GMMs is not implemented")
        
        self.gmm = gmm
        self.rho = rho_total
        self.rho_dB = rho_deltaB
        self.rho_dW = rho_deltaW
        self.gmm_weights = gmm_weights
        self.mean_only = mean_only
        
        if avg_periods is not None:
            self.imts = [SA(t) if t > 0 else PGA() for t in avg_periods]
        else:
            self.imts = None
        

    def compute(
            self,
            ctx: np.recarray,
            imts: IMT|list[IMT]|None, 
            mean: np.ndarray,
            sigma: np.ndarray|None=None,
            tau: np.ndarray|None=None,
            phi: np.ndarray|None=None) -> tuple[np.ndarray, np.ndarray]:
        
        # mean sigma tau and phi are all 1 x len(ctx) arrays instead of the
        # standard len(imts) x len(ctx) because the compute function should
        # always return None and the imts are averaged to produce a single value 
        
        if self.imts is not None:
            # the stored imts are used
            imts = self.imts

        # the standard arrays of shape len(imts) x len(ctx) still need to be 
        # to the openquake GMPE to calculate all the im values for the ctx
        sa_mean = np.zeros((len(imts), len(ctx)))
        sa_sigma = np.zeros_like(sa_mean)
        sa_tau = np.zeros_like(sa_mean)
        sa_phi = np.zeros_like(sa_mean)

        # for each GMM get the mean and std. dev for all ruptures
        for g in self.gmm:
            g.compute(ctx, imts, sa_mean, sa_sigma, sa_tau, sa_phi)

        # calculate the mean and std. dev of the ln_AvgSA
        mu_ln_AvgSA = compute_ln_AvgSA(sa_mean.T) # transpose so columns are for each T
        mean[:,:] = mu_ln_AvgSA.T

        if not self.mean_only:
            if self.use_dB_dW:
                tau_ln_AvgSA = compute_std_dev_ln_AvgSA(sa_tau.T, self.rho_dB)
                phi_ln_AvgSA = compute_std_dev_ln_AvgSA(sa_phi.T, self.rho_dW)
                sigma_ln_AvgSA = np.sqrt(tau_ln_AvgSA ** 2 + phi_ln_AvgSA ** 2)
                tau[:,:] = tau_ln_AvgSA.T
                phi[:,:] = phi_ln_AvgSA.T
            else:
                sigma_ln_AvgSA = compute_std_dev_ln_AvgSA(sa_sigma.T, self.rho)
            
            sigma[:,:] = sigma_ln_AvgSA.T
        return
    
    def __str__(self):
        return f"[indirect_AvgSA_GMPE with {self.gmm}]"
