# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "SpaceVecAlg"
copyright = "2025, CNRS-AIST JRL, LIRMM"
author = "CNRS-AIST JRL, LIRMM"
release = "2.12.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",  # For Google/NumPy docstring support
]

templates_path = ["_templates"]
exclude_patterns = []

# autodoc options
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'

html_theme = "sphinx_book_theme"
html_theme_options = {
    "repository_url": "https://github.com/jrl-umi3218/SpaceVecAlg",
    "use_repository_button": True,
    "use_download_button": False,
}
html_title = ""
# html_logo = "_static/images/mc_rtc_logo.png"

html_static_path = ["_static"]
