"""
Sphinx configuration for the pickagm documentation.

Build locally (from the repository root) with:

    pip install ".[docs]"
    sphinx-build -W -b html docs docs/_build/html

`-W` turns every warning into an error, so a clean build means no warnings.
Read the Docs runs the same build using .readthedocs.yaml.

The pages are written in Markdown and read by MyST. The API reference is
generated from the installed package by autodoc/autosummary, so pickagm and
its real dependencies (including OpenQuake) must be importable. Nothing is
mocked: if a module cannot be imported, the build fails, exactly as it would
for a user.
"""

import inspect
from importlib.metadata import version as _pkg_version

# -- Project information -----------------------------------------------------

project = "pickagm"
author = "Nicholas Clemett"
copyright = "Nicholas Clemett"

# The version comes from the installed package metadata (pyproject.toml).
# The package has no __version__ attribute, so this is the one place that
# knows it. `release` is the full version string shown in the docs; `version`
# is the short X.Y form.
release = _pkg_version("pickagm")
version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------------

extensions = [
    # Markdown pages, with MyST extras enabled below.
    "myst_parser",
    # Pull signatures and docstrings out of the installed package.
    "sphinx.ext.autodoc",
    # Build the API summary tables and one page per module (api/index.md).
    "sphinx.ext.autosummary",
    # Understand Google/NumPy style docstrings as well as reST field lists.
    "sphinx.ext.napoleon",
    # Add "[source]" links that show the highlighted source of each object.
    "sphinx.ext.viewcode",
    # Render maths with MathJax in the browser (used by $...$ in MyST).
    "sphinx.ext.mathjax",
    # Add a "copy" button to every code block.
    "sphinx_copybutton",
]

# Only Markdown and reST sources; Markdown is read by MyST.
source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}

# The landing page (docs/index.md).
root_doc = "index"

# Files and folders that are not part of the documentation.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- MyST (Markdown) -----------------------------------------------------------

myst_enable_extensions = [
    # $inline$ and $$display$$ maths in Markdown.
    "dollarmath",
    # LaTeX environments such as \begin{aligned} ... \end{aligned}.
    "amsmath",
    # ::: fenced blocks for notes and warnings (admonitions).
    "colon_fence",
]
# Give every heading up to level 3 an anchor, so pages can link to sections.
myst_heading_anchors = 3

# -- autodoc / autosummary ----------------------------------------------------

# Write the per-module API pages (into docs/api/generated/) at build time.
autosummary_generate = True
# Document only objects defined in each module, not names it imports from
# elsewhere (otherwise, for example, SCALING_EXPONENTS would appear in both
# `distributions` and `selection`).
autosummary_imported_members = False

# Options applied to every automodule/autoclass directive:
#   members          document the public functions and classes of a module;
#   undoc-members    include them even when they have no docstring (most do
#                    not yet, and their signatures are still useful);
#   show-inheritance show the base class of each class.
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}
# Keep the order of the source file, which follows the logic of each module,
# instead of sorting alphabetically.
autodoc_member_order = "bysource"
# Show type hints in the signature, where the code already has them.
autodoc_typehints = "signature"
# Show the class docstring only. Several __init__ docstrings use free-form
# notes that are not valid reST, so they are left out of the build.
autoclass_content = "class"


def _skip_type_docstring_of_data(app, what, name, obj, options, lines):
    """Drop the docstring that a module-level value only borrows from its type.

    For a value such as ``CORR_MODELS = {...}`` (listed with ``autodata`` in
    api/index.md), autodoc falls back to ``obj.__doc__``, which is the
    docstring of the built-in type (here ``dict``). That text is not meant
    for these pages and is not valid reST, so it would fail the ``-W``
    build. A doc comment written for the value itself (``#: ...`` above it
    in the code) comes from the source instead, differs from ``obj.__doc__``
    and is kept.
    """
    if what != "data":
        return
    own_doc = getattr(obj, "__doc__", None)
    if own_doc and "\n".join(lines).strip() == inspect.cleandoc(own_doc).strip():
        lines.clear()


def setup(app):
    """Register the docstring hook above with Sphinx."""
    app.connect("autodoc-process-docstring", _skip_type_docstring_of_data)


# -- Napoleon ------------------------------------------------------------------

# Most existing docstrings use reST field lists (":param x: ..."), which
# Sphinx reads natively. Napoleon additionally reads the Google-style
# "Args:/Returns:/Raises:" sections used in a few places.
napoleon_google_docstring = True
napoleon_numpy_docstring = True

# -- HTML output ---------------------------------------------------------------

html_theme = "furo"
html_title = f"pickagm {release}"
# Furo shows an "edit/view source" link; point it at the GitHub repository.
html_theme_options = {
    "source_repository": "https://github.com/nc-hsu/pickagm/",
    "source_branch": "main",
    "source_directory": "docs/",
}

# -- sphinx-copybutton -----------------------------------------------------------

# Strip shell and Python prompts when a reader copies a code block.
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
