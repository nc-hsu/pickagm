import pandas as pd
import numpy as np
from typing import Protocol
from openquake.hazardlib.gsim.base import GMPE

TRT = str  # TectonicRegionType: alias for readability
IMTName = str  # Name of imt. from IMT.name: alias for readability
GmmMap = dict[TRT, dict[IMTName, GMPE]]  # alias for readability
CorrMap = dict[TRT, pd.DataFrame]  # alias for readability

