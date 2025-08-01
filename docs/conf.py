# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath(".."))


# -- Project information -----------------------------------------------------

project = "gmso"
copyright = "2024, mosdef-hub, Vanderbilt University"
author = "Matt Thompson, Alex Yang, Ray Matsumoto, Parashara Shamaprasad, Umesh Timalsina, Co D. Quach, Ryan S. DeFever, Justin Gilmer, Nicholas C. Craven, Christopher R. Iacovella, Brad Crawford, and Chris Jones"

# The full version, including alpha/beta/rc tags
version = "0.14.0"
release = "0.14.0"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosummary",
]

autosummary_generate = True
autodoc_default_flags = ["members", "inherited-members"]

napoleon_numpy_docstring = True
napoleon_use_admonition_for_notes = True
napoleon_use_param = False
napoleon_use_ivar = True

intersphinx_mapping = {"python": ("https://docs.python.org/3.11", None)}
# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]

_python_doc_base = "http://docs.python.org/3.11"

source_suffix = ".rst"

master_doc = "index"
