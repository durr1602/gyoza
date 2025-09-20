# Configuration file for the Sphinx documentation builder.

import sys, os

from setuptools_scm import get_version

sys.path.insert(0, os.path.abspath("../../workflow/scripts/"))

# -- Project information

project = 'gyōza'
copyright = '2025, Durand'
author = 'Durand'
version = get_version(root="../..", relative_to=__file__)

if os.environ.get("READTHEDOCS") == "True":
    version = ".".join(version.split(".")[:3])

release = version

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx_tabs.tabs',
    'myst_parser',
    'autodoc2'
]

autodoc2_packages = [
    "../../workflow/scripts",
]

autodoc2_output_dir = "api"
autodoc2_render_plugin = "myst"


# allow md docstring
myst_enable_extensions = [
    "fieldlist",
]

# fix duplicate object description issue for class attributes (e.g.: dataclass)
napoleon_use_ivar = True

# Napoleon settings for NumPy-style
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_admonition_for_notes = False

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'