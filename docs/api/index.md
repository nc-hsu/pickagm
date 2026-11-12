# API reference

These pages are generated from the code of the installed package. In
version 0.1.0 most functions and classes do not have docstrings yet, so many
entries show only the name and signature (the arguments and, where the code
gives them, their types). Each entry has a **[source]** link to its code.

The modules are listed in the order in which a record selection uses them
(see {doc}`../workflow`), followed by the helper modules.

```{eval-rst}
.. autosummary::
   :toctree: generated

   pickagm.contexts
   pickagm.corrmodels
   pickagm.avgSA
   pickagm.distributions
   pickagm.eqdbases
   pickagm.selection
   pickagm.dfops
   pickagm.interpolate
   pickagm.typing
```

## Constants and type aliases

Autodoc only lists module-level values that have a docstring or doc
comment, so the ones used throughout the guide are listed here by hand,
with their values.

```{eval-rst}
.. autodata:: pickagm.typing.TRT
.. autodata:: pickagm.typing.IMTName
.. autodata:: pickagm.typing.GmmMap
.. autodata:: pickagm.typing.CorrMap
.. autodata:: pickagm.corrmodels.CORR_MODELS
.. autodata:: pickagm.corrmodels.CORR_COEFFS_FOLDER
   :no-value:
.. autodata:: pickagm.distributions.SCALING_EXPONENTS
.. autodata:: pickagm.distributions.NormalRV
.. autodata:: pickagm.eqdbases.esm_nonSA_im_translations
.. autodata:: pickagm.eqdbases.esm_unit_conversions
```

`CORR_COEFFS_FOLDER` is the `im_corr_data` folder next to `corrmodels.py`
in the installed package; its full path depends on where pickagm is
installed.
