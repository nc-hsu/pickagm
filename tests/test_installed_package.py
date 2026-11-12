"""
Smoke test for an *installed* copy of pickagm.

Why this test exists
--------------------
pickagm reads two kinds of bundled data files at run time:

* the ESM flatfile (``data/ESM_flatfile_SA.csv``), found through
  ``importlib.resources.files("pickagm.data")`` in ``eqdbases.py``;
* the correlation coefficients (``im_corr_data/clemett_corr_coeffs.hdf5``),
  found through ``Path(__file__).parent / "im_corr_data"`` in
  ``corrmodels.py``.

With an editable install (``pip install -e .``) Python imports the package
straight from ``src/pickagm/``, where every file is always present. So
packaging mistakes, such as a data file that was never declared as package
data, stay hidden. This test is meant to be run after installing the built
wheel into a fresh virtual environment. It then checks the files that a
normal user actually receives.

How to run it
-------------
Install the wheel into a fresh environment, then run pytest from the
repository root::

    pip install dist/pickagm-0.1.0-py3-none-any.whl pytest
    python -m pytest tests/test_installed_package.py

By default the test fails if ``pickagm`` was imported from a ``src`` folder
(an editable install, or the source tree on ``sys.path``), because then it
would prove nothing about the wheel. Set the environment variable
``PICKAGM_ALLOW_EDITABLE=1`` to run the remaining checks against an
editable install anyway; the check that ``BB17.csv`` is absent is then
skipped, since that file is present in the source tree by design.
"""

import importlib
import importlib.resources as res
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import pickagm

# Every module in the package. Importing each one also checks that the
# runtime dependencies (numpy, pandas, scipy, h5py, openquake.hazardlib) are
# installed, because the modules import them at the top.
MODULES = [
    "avgSA",
    "contexts",
    "corrmodels",
    "dfops",
    "distributions",
    "eqdbases",
    "interpolate",
    "selection",
    "typing",
]

# The four registered correlation models, by their key in
# corrmodels.CORR_MODELS. Listed here as well so that the test notices if
# the registry ever changes.
EXPECTED_CORR_MODELS = {
    "clemett_asc",
    "clemett_sinter",
    "clemett_sslab",
    "clemett_vrancea",
}

# True when the user has asked to allow an editable/source-tree import.
ALLOW_EDITABLE = os.environ.get("PICKAGM_ALLOW_EDITABLE") == "1"

# The folder that pickagm was imported from, e.g.
# .../site-packages/pickagm (installed) or .../src/pickagm (source tree).
PACKAGE_DIR = Path(pickagm.__file__).resolve().parent


def _imported_from_src() -> bool:
    """True if pickagm was imported from a folder called ``src``."""
    return PACKAGE_DIR.parent.name == "src"


def test_imported_from_site_packages():
    """pickagm must come from site-packages, not from the source tree.

    If it comes from ``src/``, the rest of this file would test the
    repository instead of the wheel, so fail loudly unless the user opted
    in with PICKAGM_ALLOW_EDITABLE=1.
    """
    if ALLOW_EDITABLE:
        pytest.skip("PICKAGM_ALLOW_EDITABLE=1: source-tree import allowed")
    assert not _imported_from_src(), (
        f"pickagm was imported from {PACKAGE_DIR}, a 'src' folder. Install "
        "the built wheel into a fresh environment, or set "
        "PICKAGM_ALLOW_EDITABLE=1 to run against an editable install."
    )
    assert PACKAGE_DIR.parent.name == "site-packages", (
        f"pickagm was imported from {PACKAGE_DIR}, not from site-packages."
    )


@pytest.mark.parametrize("name", MODULES)
def test_module_imports(name):
    """Each module of the package can be imported."""
    module = importlib.import_module(f"pickagm.{name}")
    assert module.__name__ == f"pickagm.{name}"


def test_corr_model_registry():
    """CORR_MODELS holds exactly the four Clemett models."""
    from pickagm.corrmodels import CORR_MODELS

    assert set(CORR_MODELS) == EXPECTED_CORR_MODELS


@pytest.mark.parametrize("key", sorted(EXPECTED_CORR_MODELS))
def test_corr_model_loads_hdf5(key):
    """Each correlation model can be built from the bundled HDF5 file.

    Constructing a model opens ``im_corr_data/clemett_corr_coeffs.hdf5``
    next to ``corrmodels.py``, so this fails if the HDF5 was not installed.
    The constructor needs a list of non-SA IM names and a numpy array of SA
    periods; the periods used here lie inside the range stored in the file
    (0.025 s to 8 s).
    """
    from pickagm.corrmodels import CORR_MODELS

    periods = np.array([0.1, 0.35, 1.0])
    model = CORR_MODELS[key](incl_nonSA_imts=["PGA"], incl_SA_Ts=periods)

    # "total" is the one residual type every model must provide.
    assert "total" in model.rho
    rho = model.rho["total"]
    assert isinstance(rho, pd.DataFrame)

    # One row/column for PGA plus one per requested period.
    expected_labels = ["PGA", "SA(0.1)", "SA(0.35)", "SA(1.0)"]
    assert list(rho.index) == expected_labels
    assert list(rho.columns) == expected_labels

    # A correlation matrix has ones on its diagonal and values in [-1, 1].
    assert np.allclose(np.diag(rho.to_numpy()), 1.0)
    assert np.all(np.abs(rho.to_numpy()) <= 1.0)


def test_load_esm_flatfile():
    """The bundled ESM flatfile can be read through importlib.resources."""
    from pickagm import eqdbases

    df = eqdbases.load_esm_flatfile(5)
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 5
    # A few columns that the rest of the package relies on.
    for column in ["event_id", "station_code", "U_pga", "U_T0_350"]:
        assert column in df.columns


def test_bb17_not_installed():
    """BB17.csv stays in the repository but must not be in the wheel."""
    if _imported_from_src():
        pytest.skip("source-tree import: BB17.csv is present there by design")
    from pickagm.corrmodels import CORR_COEFFS_FOLDER

    assert not (CORR_COEFFS_FOLDER / "BB17.csv").exists()
    assert not res.files("pickagm.im_corr_data").joinpath("BB17.csv").is_file()
    # The file that *should* be there, as a sanity check of the path above.
    assert (CORR_COEFFS_FOLDER / "clemett_corr_coeffs.hdf5").is_file()
