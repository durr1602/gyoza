# Configuration file for the Sphinx documentation builder.

import sys, os

sys.path.insert(0, os.path.abspath("../../workflow/scripts/"))

# -- Project information

project = 'durr106-gyoza'
copyright = '2025, Durand'
author = 'Durand'
version = '1.1.6'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx_tabs.tabs',
    "autodoc2"
]

autodoc2_packages = [
    "../../workflow/scripts",
]

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