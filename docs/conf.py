# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Imports -----------------------------------------------------------------

import configparser
import datetime
import os
import re
import semantic_version
import shutil
import sphinx_bootstrap_theme
import sys

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

docssrc_dir = os.path.abspath(os.path.join(__file__, ".."))
project_dir = os.path.dirname(docssrc_dir)
sys.path.insert(0, project_dir)

# -- Sphinx Setup ------------------------------------------------------------

def setup(app):
    # Add custom stylesheet
    app.add_css_file("css/main.css")
    app.add_js_file("js/apitoc.js")
    app.add_js_file("js/example-admonition.js")

# -- Project information -----------------------------------------------------

import gecco

project = 'GECCO'
copyright = '2020-{}, Zeller group, EMBL'.format(datetime.datetime.now().year)
author = gecco.__author__

# The parsed semantic version
semver = semantic_version.Version.coerce(gecco.__version__)
# The short X.Y version
version = "v{v.major}.{v.minor}.{v.patch}".format(v=semver)
# The full version, including alpha/beta/rc tags
release = str(semver)

# Project URLs
_parser = configparser.ConfigParser()
_parser.read(os.path.join(project_dir, "setup.cfg"))
project_urls = dict(
    map(str.strip, line.split(" = ", 1))
    for line in _parser.get("metadata", "project_urls", fallback="").splitlines()
    if line.strip()
)

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "sphinx_bootstrap_theme",
    "recommonmark"
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "monokailight"

# The name of the default role for inline references
default_role = "py:obj"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "bootstrap"

# Add any paths that contain custom themes here, relative to this directory.
html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
    # Bootswatch (http://bootswatch.com/) theme.
    "bootswatch_theme": "sandstone",
    # Choose Bootstrap version.
    "bootstrap_version": "3",
    # Tab name for entire site. (Default: "Site")
    "navbar_site_name": "Documentation",
    # HTML navbar class (Default: "navbar") to attach to <div> element.
    # For black navbar, do "navbar navbar-inverse"
    "navbar_class": "navbar",
    # Render the next and previous page links in navbar. (Default: true)
    "navbar_sidebarrel": True,
    # Render the current pages TOC in the navbar. (Default: true)
    "navbar_pagenav": False,
    # A list of tuples containing pages or urls to link to.
    "navbar_links": [
        ("Repository", _parser.get("metadata", "home-page").strip(), True),
    ],
    #  + [
    #     (k, v, True)
    #     for k, v in project_urls.items()
    #     if k in {"Builds"}
    # ],
    # Render admonition using panels instead of alerts (this is a PR under
    # review there: https://github.com/ryan-roemer/sphinx-bootstrap-theme/pull/199)
    # but hopefully it will get merged one day.
    "admonition_use_panel": True,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
html_sidebars = {
    "*": ["localtoc.html"],
    "api/*": ["localtoc.html"],
}

# If given, this must be the name of an image file (path relative to the
# configuration directory) that is the logo of the docs. It is placed at the
# top of the sidebar; its width should therefore not exceed 200 pixels.
html_logo = os.path.join("_static", "img", "gecco.png")

# If given, this must be the name of an image file (path relative to the
# configuration directory) that is the favicon of the docs. Modern browsers
# use this as the icon for tabs, windows and bookmarks. It should be a
# Windows-style icon file (.ico), which is 16x16 or 32x32 pixels large.
html_favicon = os.path.join("_static", "img", "gecco.ico")

# Hide the `source` button in the navbar.
html_show_sourcelink = False

# -- Options for imgmath extension -------------------------------------------

imgmath_image_format = "svg"

# -- Options for napoleon extension ------------------------------------------

napoleon_include_init_with_doc = True
napoleon_include_special_with_doc = False
napoleon_include_private_with_doc = False
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_rtype = False

# -- Options for autodoc extension -------------------------------------------

autoclass_content = "class"
autodoc_member_order = 'bysource'
autosummary_generate = []

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    'sklearn': ('https://scikit-learn.org/stable', None),
    'sklearn-crfsuite': ("https://sklearn-crfsuite.readthedocs.io/en/latest/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "statsmodels": ("https://tedboy.github.io/statsmodels_doc/", None),
    "biopython": ("https://biopython.org/docs/1.77/api/", None),
}

# -- Options for todo extension ----------------------------------------------

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

# -- Options for recommonmark extension --------------------------------------

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# -- Options for nbsphinx extension ------------------------------------------

nbsphinx_execute = 'auto'
