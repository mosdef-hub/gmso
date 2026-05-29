# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------

project = "gmso"
copyright = "2024, mosdef-hub, Vanderbilt University"
author = (
    "Matt Thompson, Alex Yang, Ray Matsumoto, Parashara Shamaprasad, "
    "Umesh Timalsina, Co D. Quach, Ryan S. DeFever, Justin Gilmer, "
    "Nicholas C. Craven, Christopher R. Iacovella, Brad Crawford, and Chris Jones"
)

version = "0.16.1"
release = "0.16.1"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx_copybutton",
    "sphinx_autodoc_typehints",
    "numpydoc",
]

# -- Autodoc settings --------------------------------------------------------

autodoc_default_options = {
    "members": True,
    "inherited-members": False,
    "undoc-members": False,
    "show-inheritance": True,
    "member-order": "bysource",
    "exclude-members": "model_config, model_fields",
}

autodoc_mock_imports = ["foyer"]

autodoc_typehints = "description"
autodoc_typehints_format = "short"
autodoc_typehints_description_target = "documented"
autodoc_preserve_defaults = True

# -- Autosummary settings ----------------------------------------------------

autosummary_generate = True
autosummary_imported_members = False

# -- Napoleon settings -------------------------------------------------------

napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_references = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_use_ivar = False
napoleon_preprocess_types = True
napoleon_attr_annotations = True

# -- Numpydoc settings -------------------------------------------------------

numpydoc_show_class_members = False
numpydoc_class_members_toctree = False

# -- Intersphinx mapping -----------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3.11", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "sympy": ("https://docs.sympy.org/latest", None),
    "unyt": ("https://unyt.readthedocs.io/en/stable", None),
    "networkx": ("https://networkx.org/documentation/stable", None),
    "pydantic": ("https://docs.pydantic.dev/latest", None),
}

# -- sphinx-autodoc-typehints settings ---------------------------------------

always_document_param_types = False
typehints_fully_qualified = False
simplify_optional_unions = True

# -- General settings --------------------------------------------------------

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

master_doc = "index"

# -- HTML output settings ----------------------------------------------------

html_theme = "sphinx_rtd_theme"

html_theme_options = {
    "navigation_depth": 4,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "includehidden": True,
    "titles_only": False,
}

html_static_path = ["_static"]

html_context = {
    "display_github": True,
    "github_user": "mosdef-hub",
    "github_repo": "gmso",
    "github_version": "main",
    "conf_py_path": "/docs/",
}

# -- Copy button settings ----------------------------------------------------

copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
